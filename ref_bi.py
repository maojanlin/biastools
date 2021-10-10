import argparse
import re
import pickle
import os.path
from os import path
import pysam

def main(fn_vcf, fn_sam, fn_fasta, fn_output):
    # read reference file
    ref_file = open(fn_fasta, 'r')  #TODO need to fix this
    reference = ''
    for line in ref_file:
        if not line.startswith('>'):
            reference += line.strip()
    ref_file.close()
    



    vcf_file_name = fn_output+'.chr_vcf.pickle'
    sam_file_name = fn_output+'.chr_sam.pickle'
    
    ######## LOAD THE VCF FILE ########
    # read the fn_vcf file, build the pickle file it is not exist
    in_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    count_het = 0
    if not path.exists(vcf_file_name):
        chr_vcf = {}
        for segment in in_vcf_file:
            ref_name = segment.contig
            pos      = segment.pos - 1
            ref      = segment.ref
            list_alt = segment.alts
            if len(ref) > 1 or len(list_alt[0]) > 1:  # TODO skip all the multiple bases variant site
                continue
            
            if ref_name not in chr_vcf.keys():
                chr_vcf[ref_name] = [[pos], [[ref, ''.join(list_alt)]]]
            else:
                chr_vcf[ref_name][0].append(pos)
                chr_vcf[ref_name][1].append([ref, ''.join(list_alt)])
        
        #pickle.dump(chr_vcf, open(vcf_file_name, 'wb'))
        print ('Dump to {}'.format(vcf_file_name))
    else:
        #print ('Load from {}'.format(vcf_file_name))
        chr_vcf = pickle.load(open(vcf_file_name, "rb"))
    in_vcf_file.close()
    


    ######## LOAD THE SAM/BAM FILE ########
    in_sam_file = pysam.AlignmentFile(fn_sam, "r")
    chr_sam = {}
    if not path.exists(sam_file_name):
        for segment in in_sam_file:
            flag = segment.flag
            if (flag & 4): # bitwise AND 4, segment unmapped
                continue
            ref_name = segment.reference_name
            #: ignores haplotype suffixes
            if ref_name.endswith('A') or ref_name.endswith('B'):
                ref_name = ref_name[:-1]
            cigar = segment.cigarstring #### maybe cigartuple would be enough
            parsed_cigar = segment.cigartuples
            mapq  = segment.mapping_quality
            start_pos = segment.reference_start # start position in genome coordiante
            sequence  = segment.query_alignment_sequence
            rg_tag    = segment.get_tag("RG")
            mod_sequence = ""
            # warp the read sequence according to CIGAR
            mod_id = 0
            start = start_pos
            num_del = sum([p[1] for p in parsed_cigar if p[0] == 2]) # collect the Ds
            num_ins = sum([p[1] for p in parsed_cigar if (p[0] == 1 or p[0]==4 or p[0] == 5)])
            for pair_info in parsed_cigar:
                code = pair_info[0]
                runs = pair_info[1]
                if code == 0: # M
                    mod_sequence += sequence[mod_id:mod_id + runs]
                    mod_id += runs
                elif code == 1: # I
                    mod_id += runs
                elif code == 2: # D
                    mod_sequence += '-'*runs
                elif code == 4: # S
                    mod_id += runs
                elif code == 5: # H
                    mod_id += runs
                else:
                    print ("ERROR: unexpected cigar code", segment.cigarstring)
            try:
                assert len(mod_sequence) == (segment.reference_end - start_pos)
            except:
                print("WARNING! Warping fail at:", segment.query_name)

            
            if ref_name not in chr_sam.keys():
                chr_sam[ref_name] = [[], [], [], [], [], []]
            chr_sam[ref_name][0].append(start_pos) #position
            chr_sam[ref_name][1].append(mod_sequence) #sequence
            chr_sam[ref_name][2].append(mapq)
            chr_sam[ref_name][3].append(rg_tag)
        #pickle.dump(chr_sam, open(sam_file_name, 'wb'))
        print ('Dump to {}'.format(sam_file_name))
    else:
        #print ('Load from {}'.format(sam_file_name))
        chr_sam = pickle.load(open(sam_file_name, 'rb'))
    in_sam_file.close()


    ######## OUTPUT THE BIAS REPORT ########
    f = open(fn_output, 'w')
    f.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tGAP_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\n")
    
    chr_list = sorted(chr_vcf.keys())
    for ref_name in chr_list:
        list_het_site = chr_vcf[ref_name][0]
        list_het_info = chr_vcf[ref_name][1]
        
        list_sam_pos   = chr_sam[ref_name][0]
        list_sam_reads = chr_sam[ref_name][1]
        list_mapq      = chr_sam[ref_name][2]
        list_rg_tag    = chr_sam[ref_name][3]
        
        start_flag = False
        start_low_bd = 0
        count_pos = 0

        total_ref_count = 0
        total_alt_count = 0
        total_gap_count = 0
        total_other_count = 0
        hap_map = ref_hap_map(fn_vcf)

        for target_pos in list_het_site:
            ref_count = 0
            alt_count = 0
            gap_count = 0
            other_count = 0
            reads_at_het = 0
            sum_mapq = 0
            # reference_hap = find_ref_hap(pos, fn_vcf)#added this
            reference_hap = hap_map[target_pos] # replace reading every time to a map
            count_a = 0.0#added this
            count_b = 0.0#added this
            read_dis = 0.0
            if target_pos == 'N/A':
                f.write('\n')
                count_pos += 1
                continue

            for idx in range(start_low_bd, len(list_sam_pos)):
                align_st_pos = list_sam_pos[idx]
                ran = range(align_st_pos, align_st_pos + len(list_sam_reads[idx]))

                if target_pos in ran:
                    if not start_flag:
                        start_flag = True
                        start_low_bd = idx
                    reads_at_het += 1
                    sum_mapq += list_mapq[idx]
                    if 'hapA' in list_rg_tag[idx]:
                        count_a += 1.0
                    elif 'hapB' in list_rg_tag[idx]:
                        count_b += 1.0
                    try:  # scan the read base crossing the target position
                        allele = list_sam_reads[idx][target_pos - align_st_pos]
                        if allele in list_het_info[count_pos][0]:
                            ref_count += 1
                            total_ref_count += 1
                        elif allele in list_het_info[count_pos][1]:
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
                else:
                    if align_st_pos > target_pos:
                        break
            start_flag = False

            # write report for each target_position
            f.write(ref_name + "\t" + str(target_pos+1) + "\t") # report is in 1-based genome format the same as vcf
            if ref_count + alt_count == 0:
                f.write("N/A")
            else:
                f.write(format(ref_count / float(ref_count + alt_count), '.8f')) 
            f.write("\t" + str(ref_count) + "\t" + str(alt_count) + "\t" + str(gap_count) + "\t" + str(other_count) + "\t" + str(reads_at_het) + "\t" + str(sum_mapq) + "\t")
            if count_a + count_b == 0: # read distribution section
                f.write("N/A")
            else:
                if reference_hap =='hapA':
                    f.write(format((count_a)/(count_a+count_b), '.8f'))
                else:
                    f.write(format((count_b)/(count_a+count_b), '.8f'))
            f.write("\n")
            count_pos += 1
        
        # report the total count
        print()
        print(ref_name)
        print("total ref base:", total_ref_count)
        print("total alt base:", total_alt_count)
        print("total gap:", total_gap_count)
        print("total other base:", total_other_count)
    f.close()



def ref_hap_map(fn_vcf):
    """
    dict_hap_map{}
    - key: het_site
    - value: hapA or hapB
    """
    dict_hap_map = {}
    file_in = open(fn_vcf, 'r')
    for line in file_in:
        if line.startswith("#"):
            continue
        else:
            spl = line.split()
            site = int(spl[1]) - 1 # transform the 1-base into 0-base
            if dict_hap_map.get(site):
                continue
            if spl[9][0] == '0':
                dict_hap_map[site] = 'hapA'
            else:
                dict_hap_map[site] = 'hapB'
    file_in.close()
    return dict_hap_map



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_sam = args.sam
    fn_fasta = args.fasta
    fn_output = args.out
    print('Input vcf:', fn_vcf)
    print('Input sam:', fn_sam)
    print('Input fasta:', fn_fasta)
    print('Output file:', fn_output)
    main(fn_vcf, fn_sam, fn_fasta, fn_output)
