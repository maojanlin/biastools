import argparse
import re
from indelpost import Variant, VariantAlignment
import pickle
import os.path
from os import path
import pysam
import numpy as np


def parse_vcf_to_dict(fn_vcf):
    """ load the vcf file to a dict{list_info} according to chromosome name """
    print("Start parsing the vcf file", fn_vcf)
    in_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    count_het = 0
    dict_chr_vcf = {}
    for segment in in_vcf_file:
        ref_name = segment.contig
        pos      = segment.pos
        ref      = segment.ref
        list_alt = segment.alts
        hap_info = str(segment).split()[9] # "0|0", "1|0", "0|1" tag
        SNP_flag = True
        if len(ref) > 1 or len(list_alt[-1]) > 1:  # separate SNP and gap sites
            SNP_flag = False
        if ref_name not in dict_chr_vcf.keys():
            dict_chr_vcf[ref_name] = [[pos, ref, list_alt, hap_info, SNP_flag]]
        else:
            dict_chr_vcf[ref_name].append([pos, ref, list_alt, hap_info, SNP_flag])
    in_vcf_file.close()
    return dict_chr_vcf
   

def parse_sam_to_dict(fn_sam, reference): # utilizing indelpost decomposition
    """ load the sam file to a dict{list_info} according to chromosome name """
    print("Start parsing the sam file", fn_sam)
    in_sam_file = pysam.AlignmentFile(fn_sam, "r")
    dict_chr_sam = {}
    for segment in in_sam_file:
        flag = segment.flag
        if (flag & 4): # bitwise AND 4, segment unmapped
            continue
        ref_name = segment.reference_name
        if ref_name.endswith('A') or ref_name.endswith('B'): # ignores haplotype suffixes
            ref_name = ref_name[:-1]
            
        seq_name  = segment.query_name
        #cigar_tuples = segment.cigartuples
        ref_seq = segment.get_reference_sequence()
        start_pos = segment.reference_start # start position in genome coordiante, need +1 for vcf coordinate
        end_pos   = start_pos + len(ref_seq)
        mapq      = segment.mapping_quality
        rg_tag    = segment.get_tag("RG")
        
        sequence  = segment.query_alignment_sequence # aligned sequence without SoftClip part
        if ref_name not in dict_chr_sam.keys(): # initialize the dictionary according to chromosome name
            dict_chr_sam[ref_name] = [] 

        list_var = []
        if ref_seq != sequence:
            v = Variant(ref_name, start_pos+1, ref_seq, sequence, reference)
            for var in v.decompose_complex_variant():
                list_var.append((var.pos, var.ref, var.alt))
        dict_chr_sam[ref_name].append((start_pos+1, end_pos+1, mapq, rg_tag, list_var))
    in_sam_file.close()
    return dict_chr_sam




def fetch_nearby_cohort(
    var: pysam.VariantRecord,
    f_query_vcf: pysam.VariantFile,
    f_fasta: pysam.FastaFile,
    update_info: dict,
    query_info: dict=None,
    padding: int=0,
    debug: bool=False
    ) -> pysam.VariantRecord:
    ''' Fetch nearby cohorts and local REF haplotype for a variant.
    Inputs:
        - var: Target variant.
        - f_query_vcf: Queried VCF file (e.g. a cohort VCF).
        - f_fasta: REF FASTA file.
        - update_info: VCF INFO field to update. 'ID', 'Number', 'Type' are required.
            E.g.
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency estimate for each alternate allele'}
        - query_info: VCF INFO field to query. If not set, use `update_info`.
    Raises:
        - ValueError: If fetched variants don't share the same contig.
    '''
    # var.start: 0-based; var.pos: 1-based
    # Pysam uses 0-based
    # var_region = (var.contig, var.start, var.start + max(var.alleles))
    # Fetch cohort variants

    
    # COMPARE VCF AND SAM FILE FOR VARIANTs
    print("Start comparing the vcf and sam file")
    dict_chr_bias = {}
    for ref_name, list_sam_info in sorted(dict_chr_sam.items()):
        list_het_info = dict_chr_vcf[ref_name]
        dict_chr_bias[ref_name] = {}
        for het_info in list_het_info:
            # the list store the vcf information for ref_count, alt_count, gap_count, other_count, num_read, sum_mapq, hap_alt_count, hap_oth_count, SNP_flag
            dict_chr_bias[ref_name][het_info[0]] = [0,0,0,0,0,0,0,0,het_info[4]] 

        vcf_flag = False
        vcf_low_bd = 0
        for read_info in list_sam_info:
            start_pos, end_pos, mapq, rg_tag, list_var = read_info
            for idx in range(vcf_low_bd, len(list_het_info)):
                var_st_pos = list_het_info[idx][0]
                if var_st_pos in range(start_pos, end_pos): # if the vcf site is within the range of the read:
                    if not vcf_flag:
                        vcf_flag = True
                        vcf_low_bd = idx
                elif var_st_pos > end_pos: # the vcf site exceed the read
                    break
            vcf_flag = False
            list_var_cover = list_het_info[vcf_low_bd:idx]
            list_var_flag  = np.zeros(len(list_var_cover)) # 0:ref, 1:alt, 2:gap, 3:others
            
            # All-to-all comparison with indelPost
            for var_sam_info in list_var:
                s_pos, s_ref, s_alt = var_sam_info
                '''
                v1 = Variant("chr21", 14238195, 'TAATCCT', 'CTGGGCC',  reference)
                v2 = Variant("chr21", 14238193, 'TCTAATCCT', 'TCCTGGGCC',  reference)
                print("equivalent?", v1 == v2)
                '''
                bk_flag = False
                # mismatches and insertions: simple comparison
                for id_vcf, var_vcf_info in enumerate(list_var_cover):
                    v_pos, v_ref, v_list_alt, v_hap_info, SNP_flag = var_vcf_info
                    if s_pos == v_pos:
                        if s_alt in v_list_alt:
                            list_var_flag[id_vcf] = 1
                        else:
                            list_var_flag[id_vcf] = 3
                        bk_flag = True
                        break
                if len(s_ref) <= len(s_alt) or bk_flag:
                    pass
                else: # deletions: need to consider gaps
                    for id_vcf, var_vcf_info in enumerate(list_var_cover):
                        if var_vcf_info in range(s_pos+1,s_pos+len(s_ref)-len(s_alt)+1):
                            if list_var_flag[id_vcf] == 0:
                                list_var_flag[id_vcf] =2
            
            for id_vcf, flag_var in enumerate(list_var_flag):
                v_pos, _, _, v_hap_info, _ = list_var_cover[id_vcf]
                dict_chr_bias[ref_name][v_pos][int(flag_var)] += 1  # variant information
                dict_chr_bias[ref_name][v_pos][4] += 1              # read count
                dict_chr_bias[ref_name][v_pos][5] += mapq           # mapq
                if 'hapA' in rg_tag:                                # hap info
                    if v_hap_info[0] =='0': # hapA
                        dict_chr_bias[ref_name][v_pos][6] += 1
                    else:
                        dict_chr_bias[ref_name][v_pos][7] += 1
                elif 'hapB' in rg_tag:
                    if v_hap_info[0] =='0': # hapA
                        dict_chr_bias[ref_name][v_pos][7] += 1
                    else:
                        dict_chr_bias[ref_name][v_pos][6] += 1
                else:
                    print("CAUTION, rg_tag haplotype information incorrect!")
    return dict_chr_bias



def output_report(dict_chr_bias, fn_output):
    f_all = open(fn_output, 'w')
    f_gap = open(fn_output + '.gap', 'w')
    f_SNP = open(fn_output + '.SNP', 'w')
    f_all.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tGAP_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\tGAP\n")
    f_gap.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tGAP_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\n")
    f_SNP.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tGAP_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\n")
    output_string = ""
    for ref_name, dict_var_info in sorted(dict_chr_bias.items()):
        for var_pos, var_info in sorted(dict_var_info.items()):
            ref_count, alt_count, gap_count, other_count, num_read, sum_mapq, hap_alt_count, hap_oth_count, SNP_flag = var_info
            output_string = (ref_name + '\t' + str(var_pos) + '\t')
            if ref_count + alt_count == 0:
                output_string += ("N/A")
            else:
                output_string += (format(ref_count / float(ref_count + alt_count), '.8f')) 
            output_string += ("\t" + str(ref_count) + "\t" + str(alt_count) + "\t" + str(gap_count) + "\t" + str(other_count) + "\t" + str(num_read) + "\t" + str(sum_mapq) + "\t")
            if hap_alt_count + hap_oth_count == 0: # read distribution section
                output_string += ("N/A")
            else:
                output_string += (format((hap_alt_count)/(hap_alt_count + hap_oth_count), '.8f'))
            if SNP_flag:
                f_all.write(output_string + '\t' + '\n')
                f_SNP.write(output_string + '\n')
            else:
                f_all.write(output_string + '\t' + '.' + '\n')
                f_gap.write(output_string + '\n')
    f_all.close()
    f_gap.close()
    f_SNP.close()



def cohort_seq(
    cohort_vars: list,
    ref: str,
    padding: int=0,
    debug: bool=False
    ) -> list:
    ''' 
    Modified from matchall.match_allele
    Only return the cohort reference/alternative sequences

    Inputs:
        - cohort_vars: list of pysam.VariantRecord. Fetched nearby cohorts.
        - ref: Local REF haplotype

    Returns:
        - list_cohort_seqs: list of seqs nearby the variants
    '''
    if debug:
        print('cohorts')
        for c in cohort_vars:
            print(f'{c}'.rstrip())
            print('   ', c.start, c.stop, c.alleles)
        print('ref =', ref)
    
    #start = min(var.start, min([v.start for v in cohort_vars]))
    start = min([v.start for v in cohort_vars])
    """
    dict_hap_info = {}
    for i, alt in enumerate(var.alts):
        var_seq = ref[:var.start-start] + alt + ref[var.stop-start:]
        dict_hap_info[var_seq] = 0
        # dict_hap_info[var_seq] = 
    print("dict_hap_info",dict_hap_info)
    """
    # Loop through matched cohort variants
    list_cohort_seqs = []
    for c_var in cohort_vars:
        # print(c_var.info[info['ID']])
        # Loop through each allele (index=`i`) in a cohort variant
        list_sub_seqs = []
        for i, alt in enumerate(c_var.alts):
            c_var_seq = ref[:c_var.start-start] + alt + ref[c_var.stop-start:]
            list_sub_seqs.append(c_var_seq)
        list_cohort_seqs.append(list_sub_seqs)
    return list_cohort_seqs


def fetch_alt_seqs(
        var:pysam.VariantRecord,
        ref:str
        )-> list:
    """
    Inputs: 
        - variant
        - reference sequence
    Ouput:
        - list of alternatives sequences
    """
    list_alt_seqs = []
    for idx, alt in enumerate(var.alts):
        var_seq = ref[:var.start] + alt + ref[var.stop:]
        list_alt_seqs.append(var_seq)
    return list_alt_seqs


def parse_MD(md_tag):
    list_string = re.split('(\d+)', md_tag)
    list_chunk = []
    for ele in list_string:
        if ele.isdigit():
            list_chunk.append(('M', int(ele)))
        elif ele == "":
            continue
        elif ele[0] == '^':
            list_chunk.append(('D', ele[1:]))
        else:
            list_chunk.append(('m', ele))
    return list_chunk


def variant_boundary(
        segment: pysam.AlignedSegment,
        ) -> tuple:
    """giving the closest distatnt to boundary of the aligned read"""
    cigar_tuples = segment.cigartuples
    md_tag       = segment.get_tag("MD")
    md_tuples    = parse_MD(md_tag)
    
    # determine the left boundary
    left_bd = 0
    code, run = cigar_tuples[0]
    if code == 0 or code == 7 or code == 8: # M or = or X
        left_bd = run
    m_code, m_run = md_tuples[0]
    if m_code == 'M':
        left_bd = min(m_run, left_bd)
    else:
        left_bd = 0

    # determine the right boundary
    right_bd = 0
    code, run = cigar_tuples[-1]
    if code == 0 or code == 7 or code == 8:
        right_bd = run
    m_code, m_run = md_tuples[-1]
    if m_code == 'M':
        right_bd = min(right_bd, m_run)
    else:
        right_bd = 0

    return (left_bd, right_bd)


def compare_vcf_sam(
    f_vcf:   pysam.VariantFile,
    f_bam:   pysam.AlignmentFile,
    f_fasta: pysam.FastaFile,
    padding: int=5
    ) -> dict:
    """
    take the reads and print the related cohort variants sequences
    """
    # initialize the dictionary according to chromosome name
    # TODO but I am not sure the usage of dict_chr_sam
    dict_chr_sam = {}
    for ref_name in f_fasta.references:
        dict_chr_sam[ref_name] = []

    for segment in f_bam:
        flag = segment.flag
        if (flag & 4): # bitwise AND 4, segment unmapped
            continue
        # aligned read information
        ref_name = segment.reference_name
        seq_name  = segment.query_name
        #ref_seq = segment.get_reference_sequence()
        pos_start = segment.reference_start # start position in genome coordiante, need +1 for vcf coordinate
        pos_end   = segment.reference_end
        mapq      = segment.mapping_quality
        rg_tag    = segment.get_tag("RG")
        sequence  = segment.query_alignment_sequence # aligned sequence without SoftClip part
        
        left_bd, right_bd = variant_boundary(segment)
        left_pad  = max(0, left_bd - padding)
        right_pad = max(0, right_bd - padding)

        if right_pad == 0:
            seq_var_segment = segment.query_alignment_sequence[left_pad:]
        else:
            seq_var_segment = segment.query_alignment_sequence[left_pad: -right_pad]
        # print(segment.query_alignment_sequence)
        # print(segment.get_reference_sequence())
        print(seq_var_segment)
        # print(left_bd, right_bd)
        # print(left_pad, right_pad)

        related_vars = list(f_vcf.fetch(ref_name, pos_start, pos_end)) # list of pysam.variant
        if len(seq_var_segment) == 0: # TODO can adjust to the paddings
            print("ALL CLEAR!")
            continue
        else:
            cohort_bg = pos_start + left_pad
            cohort_ed   = pos_end  - right_pad
            cohort_vars = list(f_vcf.fetch(ref_name, cohort_bg, cohort_ed)) # list of pysam.variant
            if len(cohort_vars) == 0:
                #TODO
                print("NOT RELATED!")
                continue

            # TODO cohort_start and cohort_maxstop needs to be fixed too
            cohort_start = min(cohort_bg, min([v.start for v in cohort_vars]))
            cohort_maxstop = cohort_ed
            for v in cohort_vars:
                cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) for a in v.alleles]))
            print(cohort_bg, cohort_start, cohort_ed, cohort_maxstop)
            
            #for var in cohort_vars:
            #    print(var.pos, var.ref, var.alts)
            try:
                ref_seq = f_fasta.fetch(
                    reference=ref_name, start=cohort_start, end=cohort_maxstop)
            except:
                update_info_empty(var, update_info)
                print(f'Warning: encounter the edge of a contig. Set "{update_info["ID"]}" as the init value.', file=sys.stderr)
                return var
            # print(cohort_seq(
            #     cohort_vars=cohort_vars,
            #     ref=ref_seq,
            #     padding=padding,
            #     debug=False))
            # print("=============")
            for var in cohort_vars:
                ref_seq = f_fasta.fetch(reference=ref_name, start=var.start-5, end=var.stop+5)
                print(var.start, fetch_alt_seqs(var, ref_seq))
            print("=====================")
    return dict_chr_sam


def switch_var_seq(
        var     :pysam.VariantRecord,
        ref     :str,
        start   :int,
        genotype:int
        )-> tuple :
    """
    Switch the ref sequence according to the haplotype information
    """
    if genotype == 0:
        return ref, 0
    else:
        alt = var.alts[genotype - 1]
        return ref[:var.start-start] + alt + ref[var.stop-start:], len(var.ref) - len(alt)


def variant_seq(
        f_vcf   :pysam.VariantFile,
        f_fasta :pysam.FastaFile,
        padding :int=5
        )-> set: # set_conflict, dict_var_haps
    """
    Output
        dictionary containing the sequences nearby the variants
        - keys: ref_name
        - values: dict {}
                    - keys: var.start
                    - values: (seq_hap0, seq_hap1)
        set containing the conflict variants
        - note: not include the variants within the padding distance to conflict variants
    """
    dict_ref_var_seqs = {}
    for ref_name in f_fasta.references:
        dict_ref_var_seqs[ref_name] = {}

    set_conflict_vars = set()
    list_f_vcf = list(f_vcf)
    for var in list_f_vcf:
        ref_name = var.contig
        var_start = var.start - padding
        var_stop  = var.stop  + padding
        
        cohort_vars = list(f_vcf.fetch(var.contig, var_start, var_stop))
        if len(cohort_vars) > 1: # the case where variants in the padding area
            # Expanding to the padding variants' padding area
            cohort_start = min(var_start, min([v.start-padding for v in cohort_vars]))
            cohort_maxstop = var_stop
            for v in cohort_vars:
                cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) + padding for a in v.alleles]))
            # Iterate until there are no variants in the padding area
            while cohort_vars != list(f_vcf.fetch(var.contig, cohort_start, cohort_maxstop)):
                cohort_vars = list(f_vcf.fetch(var.contig, cohort_start, cohort_maxstop))
                cohort_start = min(cohort_start, min([v.start-padding for v in cohort_vars]))
                for v in cohort_vars:
                    cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) + padding for a in v.alleles]))

            # Iterative parameters
            ref_seq = f_fasta.fetch(reference=var.contig, start= cohort_start, end = cohort_maxstop)
            seq_hap0, seq_hap1 = ref_seq, ref_seq
            adj_hap0, adj_hap1 = cohort_start, cohort_start
            diff_hap0, diff_hap1     =  0,  0
            prev_start0, prev_start1 = -1, -1
            l_diff0 = var_start - cohort_start
            r_diff0 = cohort_maxstop - var_stop
            l_diff1 = l_diff0
            r_diff1 = r_diff0
            for c_var in cohort_vars: # Modify the iterative parameters
                hap_0, hap_1 = c_var.samples[0]['GT']
                if c_var.start > prev_start0 + diff_hap0: # checking if there are overlaps
                    adj_hap0 += diff_hap0
                    seq_hap0, diff_hap0 = switch_var_seq(c_var, seq_hap0, adj_hap0, hap_0)
                    prev_start0 = c_var.start
                    if c_var.start < var.start:
                        l_diff0 -= diff_hap0
                    elif c_var.start != var.start:
                        r_diff0 -= diff_hap0
                else: # overlapping variants are consider conflicts
                    set_conflict_vars.add(prev_start0)
                    set_conflict_vars.add(c_var.start)
                if c_var.start > prev_start1 + diff_hap1:
                    adj_hap1 += diff_hap1
                    seq_hap1, diff_hap1 = switch_var_seq(c_var, seq_hap1, adj_hap1, hap_1)
                    prev_start1 = c_var.start
                    if c_var.start < var.start:
                        l_diff1 -= diff_hap1
                    elif c_var.start != var.start:
                        r_diff1 -= diff_hap1
                else:
                    set_conflict_vars.add(prev_start1)
                    set_conflict_vars.add(c_var.start)

            seq_hap0 = seq_hap0[max(0,l_diff0):len(seq_hap0) - max(0,r_diff0)]
            seq_hap1 = seq_hap1[max(0,l_diff1):len(seq_hap1) - max(0,r_diff1)]
        else: # single variant
            ref_seq = f_fasta.fetch(reference=var.contig, start= var_start, end = var_stop)
            hap_0, hap_1 = var.samples[0]['GT']
            seq_hap0,_ = switch_var_seq(var, ref_seq, var_start, hap_0)
            seq_hap1,_ = switch_var_seq(var, ref_seq, var_start, hap_1)

        if dict_ref_var_seqs[ref_name].get((var.start)):
            print("WARNNING! Duplicate variant at contig:", var.contig, ",pos:", var.start)
        dict_ref_var_seqs[ref_name][(var.start)] = (seq_hap0, seq_hap1)
    return set_conflict_vars, dict_ref_var_seqs




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
    
    f_vcf   = pysam.VariantFile(fn_vcf)
    f_bam   = pysam.AlignmentFile(fn_sam)
    f_fasta = pysam.FastaFile(fn_fasta)
    padding = 10
    set_conflict_vars, dict_ref_var_seqs = variant_seq(
            f_vcf=f_vcf,
            f_fasta=f_fasta,
            padding=padding)
    # extend conflict set
    for pos in list(set_conflict_vars):
        for extend in range(pos-padding, pos+padding):
            set_conflict_vars.add(extend)
    for key, value in dict_ref_var_seqs['chr21'].items():
        #print(key, value)
        if key in set_conflict_vars:
            continue
        print('>' + str(key))
        print(value[1])
    """
    dict_chr_bias = compare_vcf_sam(
            f_vcf=f_vcf, 
            f_bam=f_bam, 
            f_fasta=f_fasta)"""
    #output_report(dict_chr_bias, fn_output)


