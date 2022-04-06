import argparse
import re
import pickle
import os.path
from os import path
import pysam
import numpy as np
from scipy.stats import chisquare


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



def map_read_to_ref(
    read_start   :int,
    read_end     :int,
    cigar_tuples :list
    ) -> dict:
    """
    return the mapping of ref->read position
    """
    dict_read_map = {}
    ref_curser  = read_start
    read_curser = 0
    for pair_info in cigar_tuples:
        code, runs = pair_info
        if code == 0 or code == 7 or code == 8: # M or = or X
            for pos in range(ref_curser, ref_curser + runs):
                dict_read_map[pos] = read_curser
                read_curser += 1
            ref_curser += runs
        elif code == 1: # I
            dict_read_map[ref_curser] = read_curser
            ref_curser  += 1
            read_curser += runs
        elif code == 2: # D
            for pos in range(ref_curser, ref_curser + runs):
                dict_read_map[pos] = read_curser
            read_curser += 1
            ref_curser += runs
        elif code == 4 or code == 5: # S or H, pysam already parsed
            pass
        else:
            print ("ERROR: unexpected cigar code in sequence", query_name)
    return dict_read_map


def match_hap(
        var_start   :int,
        read_map    :dict,
        seq_read    :str,
        seq_hap     :str,
        padding     :int
        ) -> bool:
    """
    Adjust and compare the two sequences
    """
    r_start = read_map[var_start]
    l_bound = r_start - padding
    r_bound = l_bound + len(seq_hap)
    if l_bound < 0:
        seq_hap = seq_hap[-l_bound:]
        l_bound = 0
    if r_bound > len(seq_read):
        seq_hap = seq_hap[:len(seq_read)-r_bound]
        r_bound = len(seq_read)
    if seq_read[l_bound:r_bound] == seq_hap:
        return True
    else:
        return False


"""
functions above are obsolete
"""

def chi_square_test(
        var_start:      int,
        list_pos_start: list
        ) -> float:
    """
    transform the list pos_start into distribution
    """
    #list_pos_start, list_pos_end = list_start_end
    if len(list_pos_start) < 2:
        return 0 # None, cannot test
    #print(var_start, list_pos_start)
    bucket_num = 5
    bucket_len = int(100/ bucket_num)
    list_count = np.zeros(bucket_num)
    for ele in list_pos_start:
        input_idx = int((var_start - ele)/bucket_len)
        if input_idx >= bucket_num:
            #print("skip")
            continue
        list_count[input_idx] += 1
    #print(list_count)
    return chisquare(list_count)[1]
    

def interval_variance(
        var_start:      int,
        list_start_end: list
        ) -> float:
    list_pos_start, list_pos_end = list_start_end
    if len(list_pos_start) < 2:
        return None
    print(var_start)
    print(list_pos_start)
    list_interval = []
    old_pos = 0 #list_pos_start[0]
    for pos in list_pos_start:
        list_interval.append(pos - old_pos)
        old_pos = pos
    list_interval = list_interval[1:]
    mean_interval  = sum(list_interval)/len(list_interval)
    print(list_interval)
    var = 0
    for interval in list_interval:
        var += (interval - mean_interval)*(interval - mean_interval)
    var = var/len(list_interval)
    return var



def output_report(
        f_vcf                   :pysam.VariantFile,
        dict_ref_bias           :dict,
        dict_set_conflict_vars  :dict,
        fn_output               :str
        ) -> None:
    """
    Output the reference bias report to three different files:
        - f_all: containing all the variants
        - f_gap: contains only insertions and deletions
        - f_SNP: contains only SNPs
    """
    f_all = open(fn_output, 'w')
    f_gap = open(fn_output + '.gap', 'w')
    f_SNP = open(fn_output + '.SNP', 'w')
    f_all.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tBOTH_COUNT\tNEITHER_COUNT\tNUM_READS\tSUM_MAPQ\tEVEN_P_VALUE\tREAD_DISTRIBUTION\tGAP\n")
    f_gap.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tBOTH_COUNT\tNEITHER_COUNT\tNUM_READS\tSUM_MAPQ\tEVEN_P_VALUE\tREAD_DISTRIBUTION\n")
    f_SNP.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tBOTH_COUNT\tNEITHER_COUNT\tNUM_READS\tSUM_MAPQ\tEVEN_P_VALUE\tREAD_DISTRIBUTION\n")
    for var in f_vcf:
        ref_name = var.contig
        hap = var.samples[0]['GT']
        # Filtering all the homozygous alleles or the alleles without reference
        if (hap[0] != 0 and hap[1] != 0) or (hap[0] == 0 and hap[1] == 0):
            continue
        if hap[0] == 0:
            idx_ref, idx_alt = 0, 1
        else:
            idx_ref, idx_alt = 1, 0
        # Filtering the conflict vars
        if var.start in dict_set_conflict_vars[ref_name]:
            continue
        n_read = dict_ref_bias[ref_name][var.start]['n_read']
        n_var  = dict_ref_bias[ref_name][var.start]['n_var']
        map_q  = dict_ref_bias[ref_name][var.start]['map_q']
        #p_value = interval_variance(var.start, dict_ref_bias[ref_name][var.start]['distribute'])
        p_value = chi_square_test(var.start, dict_ref_bias[ref_name][var.start]['distribute'][idx_alt])
        if p_value: # p_value is not None
            p_value = min(p_value, chi_square_test(var.start, dict_ref_bias[ref_name][var.start]['distribute'][idx_ref]))
        output_string = (ref_name + '\t' + str(var.start+1) + '\t')
        # n_var[0,1,2,3] = hap0, hap1, both, others
        if sum(n_var[:2]) == 0:
            output_string += ("N/A")
        else:
            #output_string += (format((n_var[idx_ref]+0.5*n_var[2]) / float(sum(n_var[:3])), '.8f'))
            output_string += (format((n_var[idx_ref]) / float(sum(n_var[:2])), '.8f'))
        output_string += ("\t" + str(n_var[idx_ref]) + "\t" + str(n_var[idx_alt]) + "\t" + str(n_var[2]) +"\t" + str(n_var[3]) + "\t" + str(sum(n_read)) + "\t" + str(sum(map_q)) + "\t")
        if p_value == None:
            output_string += ("0") + '\t'
        else:
            output_string += (format(p_value, '.8f')) + '\t'
        if sum(n_read) == 0:
            output_string += ("N/A")
        else:
            output_string += (format(n_read[idx_ref] / float(sum(n_read)), '.8f'))

        if len(var.ref) ==  len(var.alts[ hap[idx_alt] - 1]): # length of ref is equal to length of 
            f_all.write(output_string + '\t' + '\n')
            f_SNP.write(output_string + '\n')
        else:
            f_all.write(output_string + '\t' + '.\n')
            f_gap.write(output_string + '\n')
    
    f_all.close()
    f_gap.close()
    f_SNP.close()


def hap_inside(
        seq_read    :str,
        seq_hap     :str,
        padding     :int
        ) -> bool:
    """
    Finding if the haplotype is in the read
    Also considering the boundary condition
    One padding side can be omitted
    """
    if seq_hap in seq_read:
        return True
    else:
        len_hap = len(seq_hap)
        for idx in range(1,padding):
            # checking read left side
            if seq_hap[idx:] == seq_read[:len_hap - idx]:
                return True
            # checking read right side
            if seq_hap[:-idx] == seq_read[idx - len_hap:]:
                return True
    return False


def return_locate_cigar(
        read_start  :int,
        target_pos  :int,
        cigar_tuples:tuple
        ) -> int:
    """
    return the cigar value of a location
    according to the CIGAR string
    """
    ref_curser  = read_start -1
    read_curser = 0
    for pair_info in cigar_tuples:
        code, runs = pair_info
        if code == 0 or code == 7 or code == 8: # M or = or X
            ref_curser += runs
            if ref_curser > target_pos:
                return 0
            else:
                read_curser += runs
        elif code == 1: # I
            ref_curser  += 1
            if ref_curser > target_pos:
                return -runs
            else:
                read_curser += runs
        elif code == 2: # D
            ref_curser += runs
            if ref_curser > target_pos:
                return runs
            else:
                read_curser += 1
        elif code == 4 or code == 5: # S or H, pysam already parsed
            pass
        else:
            print ("ERROR: unexpected cigar code in sequence")
    return 0


def locate_by_cigar(
        read_start  :int,
        target_pos  :int,
        cigar_tuples:tuple
        ) -> int:
    """
    return the location of a specific reference position in the read
    according to the CIGAR string
    """
    ref_curser  = read_start
    read_curser = 0
    for pair_info in cigar_tuples:
        code, runs = pair_info
        if code == 0 or code == 7 or code == 8: # M or = or X
            ref_curser += runs
            if ref_curser > target_pos:
                return read_curser + (runs - ref_curser + target_pos)
            else:
                read_curser += runs
        elif code == 1: # I
            #ref_curser  += 1
            if ref_curser > target_pos:
                return read_curser
            else:
                read_curser += runs
        elif code == 2: # D
            ref_curser += runs
            if ref_curser > target_pos:
                return read_curser
            #else:
            #    read_curser += 1
        elif code == 4 or code == 5: # S or H, pysam already parsed
            pass
        else:
            print ("ERROR: unexpected cigar code in sequence")
    return read_curser


def match_to_hap(
        seq_name    :str, # for debug
        read_start  :int,
        read_end    :int,
        var_start   :int,
        seq_read    :str,
        seq_hap     :str,
        cigar_tuples:tuple,
        padding     :int,
        l_min_req   :int,
        r_min_req   :int,
        start_flag  :bool=True
        ) -> int:
    """
    1. Find the matching point of the variant on the read
    2. Extend the padding on the read
    3. compare the read to haplotype sequences
    """
    if read_start > var_start: # Not cover
        return -1
    elif read_end < var_start: # Not cover
        return -1
    
    # locating the variant site on the read
    r_start = locate_by_cigar(
            read_start=read_start,
            target_pos=var_start,
            cigar_tuples=cigar_tuples
            )
    
    # Matching
    if start_flag:  # From var.start
        l_bound = r_start - padding
        r_bound = l_bound + len(seq_hap)
    else:           # From var.stop
        r_bound = r_start + padding
        l_bound = r_bound - len(seq_hap)

    min_match = 0 # minimum match length
    if l_bound < 0:
        seq_hap = seq_hap[-l_bound:]
        l_bound = 0
        min_match = r_min_req # minimum len to cover variant
    if r_bound > len(seq_read):
        seq_hap = seq_hap[:len(seq_read)-r_bound]
        r_bound = len(seq_read)
        if min_match != 0:
            print("WARNING! Both l_bound and r_bound exceed the read!!")
        min_match = l_min_req # minimum len to cover variant
    if r_bound - l_bound < min_match:
        return -1 # Not cover
    if seq_read[l_bound:r_bound].upper() == seq_hap.upper():
        return 1 # Match
    else:
        return 0 # Not match


def compare_sam_to_haps(
    f_vcf           :pysam.VariantFile,
    f_sam           :pysam.AlignmentFile,
    dict_ref_alts   :dict,
    dict_set_conflict_vars: dict
    ) -> dict:
    """
    Input:  f_sam file
    Output: ref bias dictionary according to variants
    """
    # build up the ref bias dictionary
    dict_ref_var_bias = {}
    for ref_name in dict_ref_alts.keys():
        dict_ref_var_bias[ref_name] = {}
        for start_pos in dict_ref_alts[ref_name]:
            # n_var has hap0, hap1, both, and others
            dict_ref_var_bias[ref_name][start_pos] = {'n_read':[0,0], 'n_var':[0,0,0,0], 'map_q':[0,0], 'distribute':[[],[],[],[]]}
    
    # parameters for pipeline design
    count_others  = [0,0]
    count_both    = [0,0]
    count_error   = [0,0]
    count_correct = [0,0]

    # scanning all the read alignments
    dict_errors = {}
    for segment in f_sam:
        flag = segment.flag
        if (flag & 4): # bitwise AND 4, segment unmapped
            continue
        # aligned read information
        ref_name     = segment.reference_name
        seq_name     = segment.query_name
        pos_start    = segment.reference_start # start position in genome coordiante, need +1 for vcf coordinate
        pos_end      = segment.reference_end
        cigar_tuples = segment.cigartuples
        mapq         = segment.mapping_quality
        rg_tag       = segment.get_tag("RG")
        read_seq     = segment.query_alignment_sequence # aligned sequence without SoftClip part
        
        related_vars = list(f_vcf.fetch(ref_name, pos_start, pos_end)) # list of pysam.variant
        #fetching the sequence in the read_seq regarding to the variant
        for var in related_vars:
            if var.start in dict_set_conflict_vars[ref_name]: # neglecting the conflict variant sites
                continue
            seq_hap0, seq_hap1, diff_hap0, diff_hap1 = dict_ref_alts[ref_name][var.start]
            if seq_hap0 == seq_hap1:
                continue

            if diff_hap0 !=0: # if hap0 is a gap:
                diff_read = return_locate_cigar(
                        read_start=pos_start, 
                        target_pos=var.start, 
                        cigar_tuples=cigar_tuples
                        )
                if diff_read == diff_hap0:
                    match_flag_0 = 1
                    match_flag_1 = 0
                else:
                    match_flag_0 = 0
                    match_flag_1 = match_to_hap(seq_name, pos_start, pos_end, var.start, read_seq, seq_hap1, cigar_tuples, 0, 1, 1, True)
            elif diff_hap1 !=0: # if hap1 is a gap:
                diff_read = return_locate_cigar(
                        read_start=pos_start, 
                        target_pos=var.start, 
                        cigar_tuples=cigar_tuples
                        )
                if diff_read == diff_hap1:
                    match_flag_0 = 0
                    match_flag_1 = 1
                else:
                    match_flag_0 = match_to_hap(seq_name, pos_start, pos_end, var.start, read_seq, seq_hap0, cigar_tuples, 0, 1, 1, True)
                    match_flag_1 = 0
            else:
                match_flag_0 = match_to_hap(seq_name, pos_start, pos_end, var.start, read_seq, seq_hap0, cigar_tuples, 0, 1, 1, True)
                match_flag_1 = match_to_hap(seq_name, pos_start, pos_end, var.start, read_seq, seq_hap1, cigar_tuples, 0, 1, 1, True)

            if match_flag_0 == 1 and match_flag_1 == 1:
                print("Both Trouble!", seq_name, var.start, seq_hap0, seq_hap1)

            # 5. Assign Values
            if match_flag_0 == -1 and match_flag_1 == -1:
                continue
            if match_flag_0 == 1 and match_flag_1 == 1:
                dict_ref_var_bias[ref_name][var.start]['n_var'][2] += 1
            elif match_flag_0 == 1:
                dict_ref_var_bias[ref_name][var.start]['n_var'][0] += 1
                # record the starting position of each read cover the variant
                dict_ref_var_bias[ref_name][var.start]['distribute'][0].append(pos_start)
                dict_ref_var_bias[ref_name][var.start]['distribute'][2].append(pos_end)
            elif match_flag_1 == 1:
                dict_ref_var_bias[ref_name][var.start]['n_var'][1] += 1
                # record the starting position of each read cover the variant
                dict_ref_var_bias[ref_name][var.start]['distribute'][1].append(pos_start)
                dict_ref_var_bias[ref_name][var.start]['distribute'][3].append(pos_end)
            else:
                dict_ref_var_bias[ref_name][var.start]['n_var'][3] += 1
            
            # standard updating of read number and mapping quality
            if 'hapA' == rg_tag:
                dict_ref_var_bias[ref_name][var.start]['n_read'][0] += 1
                dict_ref_var_bias[ref_name][var.start]['map_q'][0]  += mapq
            elif 'hapB' == rg_tag:
                dict_ref_var_bias[ref_name][var.start]['n_read'][1] += 1
                dict_ref_var_bias[ref_name][var.start]['map_q'][1]  += mapq
            else:
                print("WARNING, there is a read without haplotype information!!")

            # TODO DEBUG PURPOSE!
            if seq_hap0 != seq_hap1: # only count heterozygous site
                if (len(var.ref) == 1 and max([len(seq) for seq in var.alts]) == 1):
                    gap_flag = 0
                else:
                    gap_flag = 1
                if match_flag_0 and match_flag_1:
                    count_both[gap_flag] += 1
                elif match_flag_0 == False and match_flag_1 == False:
                    count_others[gap_flag] += 1
                elif ('hapA' == rg_tag) and match_flag_0:
                    count_correct[gap_flag] += 1
                elif ('hapB' == rg_tag) and match_flag_1:
                    count_correct[gap_flag] += 1
                else:
                    count_error[gap_flag] += 1
    
    print("count correct:", count_correct)
    print("count error:", count_error)
    print("count both:", count_both)
    print("count others:", count_others)
    return dict_ref_var_bias


def len_var_seq(
        var     :pysam.VariantRecord,
        genotype:int
        )-> tuple :
    """
    Switch the ref sequence according to the haplotype information
    """
    if genotype == 0:
        return 0, var.ref
    else:
        alt = var.alts[genotype - 1]
        return len(var.ref) - len(alt), alt


def variant_seq(
        f_vcf       :pysam.VariantFile,
        f_fasta     :pysam.FastaFile
        )-> tuple: # dict_set_conflict_vars, dict_var_haps, dict_cohort
    """
    Output
        - dict_set_conflict_vars: the dictionary marking the overlaping variants
        - dict_ref_alts:
            in each contig:
            - key: var.start
            - values: [varseq_hap0, varseq_hap1]
                # not only store the varseq but also indicating the variant length
    """
    dict_ref_alts = {}
    dict_set_conflict_vars = {}
    for ref_name in f_fasta.references:
        dict_ref_alts[ref_name] = {}
        dict_set_conflict_vars[ref_name] = set()

    old_ref_name = ""
    for var in f_vcf:
        ref_name = var.contig
        if old_ref_name != ref_name: # changing the contig
            # Reset the parameters
            overlap0,    overlap1    =  0,  0
            prev_start0, prev_start1 = -1, -1
            old_ref_name = ref_name
        
        hap_0, hap_1 = var.samples[0]['GT']
        diff_hap0, var_seq0 = len_var_seq(var, hap_0)
        diff_hap1, var_seq1 = len_var_seq(var, hap_1)
        if var.start > prev_start0 + overlap0 and var.start > prev_start1 + overlap1: # checking if there are overlaps
            dict_ref_alts[ref_name][var.start] = [var_seq0, var_seq1, diff_hap0, diff_hap1]
            # hap0
            prev_start0 = var.start
            overlap0 = len(var_seq0) - 1 if (diff_hap0 == 0) else diff_hap0
            # hap1
            prev_start1 = var.start
            overlap1 = len(var_seq1) - 1 if (diff_hap1 == 0) else diff_hap1
        else: # overlapping variants are consider conflicts
            dict_set_conflict_vars[ref_name].add(prev_start1)
            dict_set_conflict_vars[ref_name].add(var.start)
    return dict_set_conflict_vars, dict_ref_alts





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
    f_sam   = pysam.AlignmentFile(fn_sam)
    f_fasta = pysam.FastaFile(fn_fasta)
    var_chain = 25
    print("Start building the variant maps...")
    dict_set_conflict_vars, dict_ref_alts = variant_seq(
            f_vcf=f_vcf,
            f_fasta=f_fasta
            )
    # extend conflict set
    for ref_name in dict_set_conflict_vars.keys():
        for pos in list(dict_set_conflict_vars[ref_name]):
            for extend in range(pos-var_chain, pos+var_chain):
                dict_set_conflict_vars[ref_name].add(extend)
    
    print("Start comparing reads to the variant map...")
    dict_ref_bias = compare_sam_to_haps(
            f_vcf=f_vcf,
            f_sam=f_sam,
            dict_ref_alts=dict_ref_alts,
            dict_set_conflict_vars=dict_set_conflict_vars
            )
    f_vcf   = pysam.VariantFile(fn_vcf)
    print("Start output report...")
    output_report(
            f_vcf=f_vcf,
            dict_ref_bias=dict_ref_bias,
            dict_set_conflict_vars=dict_set_conflict_vars, 
            fn_output=fn_output
            )


