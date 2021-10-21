import argparse
import re
import pickle
import os.path
from os import path
import pysam


def parse_vcf_to_dict(fn_vcf):
    """ load the vcf file to a dict{list_info} according to chromosome name """
    print("Start parsing the vcf file", fn_vcf)
    in_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    count_het = 0
    chr_SNP_vcf = {}
    chr_gap_vcf = {}
    for segment in in_vcf_file:
        ref_name = segment.contig
        pos      = segment.pos - 1
        ref      = segment.ref
        list_alt = segment.alts
        hap_info = str(segment).split()[9] # "0|0", "1|0", "0|1" tag
        if len(ref) > 1 or len(list_alt[-1]) > 1:  # separate SNP and gap sites
            if ref_name not in chr_gap_vcf.keys():
                chr_gap_vcf[ref_name] = [[pos, ref, list_alt, hap_info]]
            else:
                chr_gap_vcf[ref_name].append([pos, ref, list_alt, hap_info])
        else:
            if ref_name not in chr_SNP_vcf.keys():
                chr_SNP_vcf[ref_name] = [[pos, ref, list_alt, hap_info]]
            else:
                chr_SNP_vcf[ref_name].append([pos, ref, list_alt, hap_info])
    in_vcf_file.close()
    return chr_SNP_vcf, chr_gap_vcf
   

def add_gap(dict_var_gap, idx, pattern):
    """
    dict_var_gap{}
    - keys: position
    - values: dict_pattern{}
              - keys: pattern
              - values: occurance
    """
    if dict_var_gap.get(idx):
        if dict_var_gap[idx].get(pattern):
            dict_var_gap[idx][pattern] += 1
        else:
            dict_var_gap[idx][pattern] = 1
    else:
        dict_var_gap[idx] = {pattern: 1}


def warp_by_cigar(cigar_tuples, sequence, query_name, dict_var_gap, start_pos):
    mod_sequence = ""
    mod_id = 0
    for pair_info in cigar_tuples:
        code, runs = pair_info
        if code == 0 or code == 7 or code == 8: # M or = or X
            mod_sequence += sequence[mod_id:mod_id + runs]
            mod_id += runs
        elif code == 1: # I
            add_gap(dict_var_gap, start_pos+mod_id-1, ('I', sequence[mod_id:mod_id+runs]))
            mod_id += runs
        elif code == 2: # D
            add_gap(dict_var_gap, start_pos+mod_id-1, ('D', runs))
            mod_sequence += '-'*runs
        elif code == 4 or code == 5: # S or H, pysam already parsed
            pass
        else:
            print ("ERROR: unexpected cigar code in sequence", query_name)
    return mod_sequence


def parse_sam_to_dict(fn_sam):
    """ load the sam file to a dict{list_info} according to chromosome name """
    print("Start parsing the sam file", fn_sam)
    in_sam_file = pysam.AlignmentFile(fn_sam, "r")
    chr_sam = {}
    dict_chr_var_gap = {}
    for segment in in_sam_file:
        flag = segment.flag
        if (flag & 4): # bitwise AND 4, segment unmapped
            continue
        ref_name = segment.reference_name
        if ref_name.endswith('A') or ref_name.endswith('B'): # ignores haplotype suffixes
            ref_name = ref_name[:-1]
            
        seq_name  = segment.query_name
        cigar_tuples = segment.cigartuples
        start_pos = segment.reference_start # start position in genome coordiante
        mapq      = segment.mapping_quality
        rg_tag    = segment.get_tag("RG")
        
        sequence  = segment.query_alignment_sequence # aligned sequence without SoftClip part
        if ref_name not in chr_sam.keys(): # initialize the dictionary according to chromosome name
            chr_sam[ref_name] = [[], [], [], []]
            dict_chr_var_gap[ref_name] = {}
        # warp the read sequence according to CIGAR
        mod_sequence = warp_by_cigar(cigar_tuples, sequence, seq_name, dict_chr_var_gap[ref_name], start_pos)
        try:
            assert len(mod_sequence) == (segment.reference_end - start_pos)
        except:
            print("WARNING! Warping fail at:", segment.query_name)
        
        chr_sam[ref_name][0].append(start_pos) #position
        chr_sam[ref_name][1].append(mod_sequence) #sequence
        chr_sam[ref_name][2].append(mapq)
        chr_sam[ref_name][3].append(rg_tag)
    in_sam_file.close()
    return chr_sam, dict_chr_var_gap



def compare_vcf_sam(fn_vcf, fn_sam, fn_output):
    chr_SNP_vcf, chr_gap_vcf  = parse_vcf_to_dict(fn_vcf)
    chr_sam, dict_chr_var_gap = parse_sam_to_dict(fn_sam)

    # COMPARE VCF AND SAM FILE FOR SNPs
    f = open(fn_output, 'w')
    f.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tGAP_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\n")
    
    for ref_name, list_het_info in sorted(chr_SNP_vcf.items()):
        list_sam_pos   = chr_sam[ref_name][0]
        list_sam_reads = chr_sam[ref_name][1]
        list_mapq      = chr_sam[ref_name][2]
        list_rg_tag    = chr_sam[ref_name][3]
        
        start_flag = False
        start_low_bd = 0  # last variant site's starting point

        total_ref_count = 0
        total_alt_count = 0
        total_gap_count = 0
        total_other_count = 0
        for het_info in list_het_info:
            target_pos, ref_base, list_alt, hap_info = het_info

            ref_count = 0
            alt_count = 0
            gap_count = 0
            other_count = 0
            count_reads = 0
            sum_mapq = 0
            
            count_a = 0.0
            count_b = 0.0
            read_dis = 0.0
            if target_pos == 'N/A':
                f.write('\n')
                continue

            for idx in range(start_low_bd, len(list_sam_pos)):
                align_st_pos = list_sam_pos[idx]
                ran = range(align_st_pos, align_st_pos + len(list_sam_reads[idx]))

                if target_pos in ran:
                    if not start_flag:
                        start_flag = True
                        start_low_bd = idx
                    sum_mapq += list_mapq[idx]
                    count_reads += 1
                    if 'hapA' in list_rg_tag[idx]:
                        count_a += 1.0
                    elif 'hapB' in list_rg_tag[idx]:
                        count_b += 1.0
                    try:  # scan the read base crossing the target position
                        allele = list_sam_reads[idx][target_pos - align_st_pos]
                        if allele == ref_base:
                            ref_count += 1
                            total_ref_count += 1
                        elif allele in list_alt:
                            alt_count += 1
                            total_alt_count += 1
                        elif allele == '-':
                            gap_count += 1
                            total_gap_count += 1
                        else:
                            other_count += 1
                            total_other_count += 1
                    except:
                        print("in here. pos is: ", pos, "   and align_st_pos is: ", align_st_pos)
                        print("sam[align_st_pos]: ", list_sam_reads[idx])
                elif align_st_pos > target_pos: # the aligned sequence exceed the variant site
                    break
            start_flag = False

            # write report for each target_position
            f.write(ref_name + "\t" + str(target_pos+1) + "\t") # report is in 1-based genome format the same as vcf
            if ref_count + alt_count == 0:
                f.write("N/A")
            else:
                f.write(format(ref_count / float(ref_count + alt_count), '.8f')) 
            f.write("\t" + str(ref_count) + "\t" + str(alt_count) + "\t" + str(gap_count) + "\t" + str(other_count) + "\t" + str(count_reads) + "\t" + str(sum_mapq) + "\t")
            if count_a + count_b == 0: # read distribution section
                f.write("N/A")
            else:
                if hap_info[0] =='0': # hapA
                    f.write(format((count_a)/(count_a+count_b), '.8f'))
                else: # hapB
                    f.write(format((count_b)/(count_a+count_b), '.8f'))
            f.write("\n")
        
        # report the total count
        print()
        print(ref_name)
        print("total ref base:", total_ref_count)
        print("total alt base:", total_alt_count)
        print("total gap:", total_gap_count)
        print("total other base:", total_other_count)
    f.close()

    # Prepare for the gapped reference bias analysis
    """ 
    for ref_name, list_variant in sorted(chr_gap_vcf.items()):
        for list_info_pairs in list_variant:
            pos, ref, list_alt, hap_info = list_info_pairs
            if dict_chr_var_gap[ref_name].get(pos):
                print(ref, list_alt, hap_info, dict_chr_var_gap[ref_name][pos])
    """


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    
    fn_vcf = args.vcf
    fn_sam = args.sam
    fn_fasta = args.fasta #TODO not used
    fn_output = args.out
    
    compare_vcf_sam(fn_vcf, fn_sam, fn_output)
