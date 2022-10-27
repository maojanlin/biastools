import argparse
import re
import pickle
import os.path
from os import path
import pysam
import numpy as np
import math
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
    p_value = chisquare(list_count)[1]
    if math.isnan(p_value):
        return 0
    else:
        return p_value
    

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


def get_division(num_1, num_2):
    if num_2 == 0:
        return 'nan'
    else:
        return format(num_1 / num_2, '.4f')


def output_report(
        f_vcf                   :pysam.VariantFile,
        dict_ref_bias           :dict,
        dict_set_conflict_vars  :dict,
        flag_real               :bool,
        fn_golden               :str,
        fn_output               :str
        ) -> None:
    """
    Output the reference bias report to three different files:
        - f_all: containing all the variants
        - f_gap: contains only insertions and deletions
        - f_SNP: contains only SNPs
    """
    if flag_real != True:
        with open(fn_golden, "rb") as f:
            dict_ref_var_name = pickle.load(f)

    f_all = open(fn_output, 'w')
    f_gap = open(fn_output + '.gap', 'w')
    f_SNP = open(fn_output + '.SNP', 'w')
    if flag_real:
        f_all.write("CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tOTHER\tGAP\n")
        f_gap.write("CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tOTHER\n")
        f_SNP.write("CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tOTHER\n")
    else:
        f_all.write("CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tOTHER\tMAP_BALANCE\tMAP_REF\tMAP_ALT\tMIS_MAP\tSIM_BALANCE\tSIM_REF\tSIM_ALT\tGAP\n")
        f_gap.write("CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tOTHER\tMAP_BALANCE\tMAP_REF\tMAP_ALT\tMIS_MAP\tSIM_BALANCE\tSIM_REF\tSIM_ALT\n")
        f_SNP.write("CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tOTHER\tMAP_BALANCE\tMAP_REF\tMAP_ALT\tMIS_MAP\tSIM_BALANCE\tSIM_REF\tSIM_ALT\n")
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
        p_value = min(p_value, chi_square_test(var.start, dict_ref_bias[ref_name][var.start]['distribute'][idx_ref]))

        output_string = (ref_name + '\t' + str(var.start+1) + '\t')
        output_string += (str(sum(n_read)) + "\t" + get_division(sum(map_q), sum(n_read)) + "\t" + format(p_value, '.4f') + '\t')
        # n_var[0,1,2,3] = hap0, hap1, both, others
        output_string += get_division(n_var[idx_ref], sum(n_var[:2])) + "\t" + str(n_var[idx_ref]) + "\t" + str(n_var[idx_alt]) +"\t" + str(n_var[3])
        if flag_real != True: # Golden Information
            # mapping balance information
            output_string += "\t" + get_division(n_read[idx_ref], sum(n_read[:2])) + '\t' + str(n_read[idx_ref]) + '\t' + str(n_read[idx_alt]) + '\t' + str(n_read[2])  
            read_info = dict_ref_var_name[ref_name][var.start]
            # simulation balance information
            output_string += '\t' + get_division(read_info[idx_ref+2], sum(read_info[2:])) + '\t' + str(read_info[idx_ref+2]) + '\t' + str(read_info[idx_alt+2])  

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
    #l_min_req = min(l_min_req, len(seq_read)) # make sure the min requirement is shorter than the read
    #r_min_req = min(r_min_req, len(seq_read))
    if l_bound < 0:
        seq_hap = seq_hap[-l_bound:]
        l_bound = 0
        min_match = r_min_req # minimum len to cover variant
    if r_bound > len(seq_read):
        seq_hap = seq_hap[:len(seq_read)-r_bound]
        r_bound = len(seq_read)
        if min_match != 0:
            print("WARNING! The read is totally contain in the variant!!", seq_name, "at", var_start)
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
    dict_ref_haps   :dict,
    dict_ref_gaps   :dict,
    dict_ref_cohorts:dict,
    dict_set_conflict_vars: dict, #For Debug only
    flag_real       :bool,
    fn_golden       :str,
    padding         :int=5
    ) -> dict:
    """
    Input:  f_sam file
    Output: ref bias dictionary according to variants
    """
    if flag_real != True:
        with open(fn_golden, "rb") as f:
            dict_ref_var_name = pickle.load(f)

    # build up the ref bias dictionary
    dict_ref_var_bias = {}
    for ref_name in dict_ref_haps.keys():
        dict_ref_var_bias[ref_name] = {}
        for start_pos in dict_ref_haps[ref_name]:
            # n_read is real answer: hapA, hapB, other_chr(mis-map)
            # n_var has hap0, hap1, both, and others
            # distribute is for even_p_value
            dict_ref_var_bias[ref_name][start_pos] = {'n_read':[0,0,0], 'n_var':[0,0,0,0], 'map_q':[0,0,0], 'distribute':[[],[],[],[]]}
    
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
        
        if flag_real == False:
            chr_tag, hap_tag = rg_tag.split('_')
        related_vars = list(f_vcf.fetch(ref_name, pos_start, pos_end)) # list of pysam.variant
        if len(related_vars) > 0: # extend the related_vars if there are cohort in the boundary
            new_start = pos_start
            new_end   = pos_end
            if dict_ref_cohorts[ref_name].get(related_vars[0].start):
                new_start = dict_ref_cohorts[ref_name][related_vars[0].start][0]  # cohort start
            if dict_ref_cohorts[ref_name].get(related_vars[-1].start):
                new_end   = dict_ref_cohorts[ref_name][related_vars[-1].start][1] # cohort end
            if new_start != pos_start or new_end != pos_end:
                related_vars = list(f_vcf.fetch(ref_name, new_start, new_end))
        #fetching the sequence in the read_seq regarding to the variant
        for var in related_vars:
            if var.start in dict_set_conflict_vars[ref_name]: # neglecting the conflict variant sites
                continue
            seq_hap0, seq_hap1 = dict_ref_haps[ref_name][var.start]

            # 1 for match, 0 for unmatch, -1 for not cover
            match_flag_0 = 0
            match_flag_1 = 0
            # 1. Cohort alignment
            if dict_ref_cohorts[ref_name].get(var.start): # Anchor Left
                cohort_start, cohort_stop, cohort_seq0, cohort_seq1, lpad_0, lpad_1, rpad_0, rpad_1 = dict_ref_cohorts[ref_name][var.start] 
                                                                                                                        # the left, right min-require covering at least 1bp in var
                                                                                                                        # performance better than covering all var
                match_flag_0 = match_to_hap(seq_name, pos_start, pos_end, cohort_start, read_seq, cohort_seq0, cigar_tuples, padding, lpad_0+1, rpad_0+1, True)
                match_flag_1 = match_to_hap(seq_name, pos_start, pos_end, cohort_start, read_seq, cohort_seq1, cigar_tuples, padding, lpad_1+1, rpad_1+1, True)
                if not ((match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_1 == 1 and  match_flag_0 !=1)): # Anchor Right
                    match_flag_0 = match_to_hap(seq_name, pos_start, pos_end, cohort_stop, read_seq, cohort_seq0, cigar_tuples, padding, lpad_0+1, rpad_0+1, False)
                    match_flag_1 = match_to_hap(seq_name, pos_start, pos_end, cohort_stop, read_seq, cohort_seq1, cigar_tuples, padding, lpad_1+1, rpad_1+1, False)
                if match_flag_0 == 1 and match_flag_1 == 1: # BOTH cases, determine by if there are indels in local alignment
                    if dict_ref_gaps[ref_name].get(var.start):
                        diff_hap0, diff_hap1 = dict_ref_gaps[ref_name][var.start]
                        diff_read = return_locate_cigar(
                                read_start=pos_start, 
                                target_pos=var.start, 
                                cigar_tuples=cigar_tuples
                                )
                        if diff_read == diff_hap0 and diff_read != diff_hap1:
                            match_flag_1 = 0
                        elif diff_read != diff_hap0 and diff_read == diff_hap1:
                            match_flag_0 = 0
            # 2. Local alignment
            flag_4 = False
            if match_flag_0 == match_flag_1: # both or others
                match_flag_0 = match_to_hap(seq_name, pos_start, pos_end, var.start, read_seq, seq_hap0, cigar_tuples, padding, padding+1, padding+1, True)
                match_flag_1 = match_to_hap(seq_name, pos_start, pos_end, var.start, read_seq, seq_hap1, cigar_tuples, padding, padding+1, padding+1, True)
                if not ((match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_1 == 1 and  match_flag_0 !=1)): # Anchor Right
                    match_flag_0 = match_to_hap(seq_name, pos_start, pos_end, var.stop, read_seq, seq_hap0, cigar_tuples, padding, padding+1, padding+1, False)
                    match_flag_1 = match_to_hap(seq_name, pos_start, pos_end, var.stop, read_seq, seq_hap1, cigar_tuples, padding, padding+1, padding+1, False)
                if match_flag_0 == 1 and match_flag_1 == 1: # BOTH cases, determine by if there are indels in local alignment
                    if  dict_ref_gaps[ref_name].get(var.start):
                        diff_hap0, diff_hap1 = dict_ref_gaps[ref_name][var.start]
                        diff_read = return_locate_cigar(
                            read_start=pos_start, 
                            target_pos=var.start, 
                            cigar_tuples=cigar_tuples
                            )
                        if diff_read == diff_hap0 and diff_read != diff_hap1:
                            match_flag_1 = 0
                        elif diff_read != diff_hap0 and diff_read == diff_hap1:
                            match_flag_0 = 0
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
            if flag_real: # no golden information
                dict_ref_var_bias[ref_name][var.start]['n_read'][0] += 1
                dict_ref_var_bias[ref_name][var.start]['map_q'][0]  += mapq
            else:
                if ref_name != chr_tag: # not the same chromosome
                    dict_ref_var_bias[ref_name][var.start]['n_read'][2] += 1
                    dict_ref_var_bias[ref_name][var.start]['map_q'][2] += 1
                elif dict_ref_var_name[ref_name].get(var.start) == None:
                    continue
                elif 'hapA' == hap_tag: # hapA
                    if seq_name in dict_ref_var_name[ref_name][var.start][0]: # check if the read name is in the golden set
                        dict_ref_var_bias[ref_name][var.start]['n_read'][0] += 1
                        dict_ref_var_bias[ref_name][var.start]['map_q'][0]  += mapq
                    else:
                        dict_ref_var_bias[ref_name][var.start]['n_read'][2] += 1
                        dict_ref_var_bias[ref_name][var.start]['map_q'][2] += 1
                elif 'hapB' == hap_tag: # hapB
                    if seq_name in dict_ref_var_name[ref_name][var.start][1]: # check if the read name is in the golden set
                        dict_ref_var_bias[ref_name][var.start]['n_read'][1] += 1
                        dict_ref_var_bias[ref_name][var.start]['map_q'][1]  += mapq
                    else:
                        dict_ref_var_bias[ref_name][var.start]['n_read'][2] += 1
                        dict_ref_var_bias[ref_name][var.start]['map_q'][2] += 1
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
                    elif ('hapA' == hap_tag) and match_flag_0:
                        count_correct[gap_flag] += 1
                    elif ('hapB' == hap_tag) and match_flag_1:
                        count_correct[gap_flag] += 1
                    else:
                        count_error[gap_flag] += 1
    print("count correct:", count_correct)
    print("count error:", count_error)
    print("count both:", count_both)
    print("count others:", count_others)
    return dict_ref_var_bias



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
        return ref, 0, len(var.ref)
    else:
        alt = var.alts[genotype - 1]
        return ref[:var.start-start] + alt + ref[var.stop-start:], len(var.ref) - len(alt), len(alt)


def variant_seq(
        f_vcf       :pysam.VariantFile,
        f_fasta     :pysam.FastaFile,
        var_chain   :int=15,
        padding     :int=5
        )-> tuple: # dict_set_conflict_vars, dict_var_haps, dict_cohort
    """
    Output
        dictionary containing the sequences nearby the variants
        - keys: ref_name
        - values: dict {}
                    - keys: var.start
                    - values: (seq_hap0, seq_hap1)
        set containing the conflict variants
        -values: dict_cohort {}
                    - keys: var.start
                    - values: (tuple)
                        - cohort start # anchor to the referene
                        - cohort seq 0 # chort seq still got paddings
                        - cohort seq 1
        - note: not include the variants within the padding distance to conflict variants
    """
    dict_ref_haps = {}
    dict_ref_gaps = {}
    dict_ref_cohorts = {}
    dict_set_conflict_vars = {}
    for ref_name in f_fasta.references:
        dict_ref_haps[ref_name] = {}
        dict_ref_gaps[ref_name] = {}
        dict_ref_cohorts[ref_name] = {}
        dict_set_conflict_vars[ref_name] = set()

    list_f_vcf = list(f_vcf)
    idx_vcf = 0 # While Loop Management
    while idx_vcf < len(list_f_vcf):
        var = list_f_vcf[idx_vcf]
        ref_name = var.contig
        
        cohort_vars = list(f_vcf.fetch(var.contig, var.start-var_chain, var.stop+var_chain))
        if len(cohort_vars) > 1: # the case where variants in the chaining area
            # Expanding to the chaining variants' chaining area
            cohort_start = min(var.start-var_chain, min([v.start-var_chain for v in cohort_vars]))
            cohort_maxstop = var.stop+var_chain
            for v in cohort_vars:
                cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) + var_chain for a in v.alleles]))
            # Iterate until there are no variants in the chaining area
            while cohort_vars != list(f_vcf.fetch(var.contig, cohort_start, cohort_maxstop)):
                cohort_vars = list(f_vcf.fetch(var.contig, cohort_start, cohort_maxstop))
                cohort_start = min(cohort_start, min([v.start-var_chain for v in cohort_vars]))
                for v in cohort_vars:
                    cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) + var_chain for a in v.alleles]))

            # Iterative parameters
            ref_seq = f_fasta.fetch(reference=var.contig, start= cohort_start, end = cohort_maxstop)
            seq_hap0, seq_hap1 = ref_seq, ref_seq
            adj_hap0, adj_hap1 = cohort_start, cohort_start
            diff_hap0, diff_hap1     =  0,  0
            overlap0,  overlap1      =  0,  0
            prev_start0, prev_start1 = -1, -1
            # parameters for cohort records
            conflict_flag = False
            # parameters keep track of the var positions
            list_start_hap = [[],[]]
            list_len_hap   = [[],[]]
            for c_var in cohort_vars: # Modify the iterative parameters
                hap_0, hap_1 = c_var.samples[0]['GT']
                if c_var.start > prev_start0 + overlap0: # checking if there are overlaps
                    adj_hap0 += diff_hap0
                    seq_hap0, diff_hap0, len_var= switch_var_seq(c_var, seq_hap0, adj_hap0, hap_0)
                    prev_start0 = c_var.start
                    overlap0 = len_var - 1 if (diff_hap0 == 0) else diff_hap0
                    list_start_hap[0].append(c_var.start - adj_hap0)
                    list_len_hap[0].append(len_var)
                else: # overlapping variants are consider conflicts
                    list_start_hap[0].append(-1)    # house keeping
                    list_len_hap[0].append(-1)      # house keeping
                    conflict_flag = True            # conflicts in the cohort
                    dict_set_conflict_vars[ref_name].add(prev_start0)
                    dict_set_conflict_vars[ref_name].add(c_var.start)
                if c_var.start > prev_start1 + overlap1:
                    adj_hap1 += diff_hap1
                    seq_hap1, diff_hap1, len_var = switch_var_seq(c_var, seq_hap1, adj_hap1, hap_1)
                    prev_start1 = c_var.start
                    overlap1 = len_var - 1 if (diff_hap1 == 0) else diff_hap1
                    list_start_hap[1].append(c_var.start - adj_hap1)
                    list_len_hap[1].append(len_var)
                else:
                    list_start_hap[1].append(-1)
                    list_len_hap[1].append(-1)
                    conflict_flag = True
                    dict_set_conflict_vars[ref_name].add(prev_start1)
                    dict_set_conflict_vars[ref_name].add(c_var.start)
                if diff_hap0 != 0 or diff_hap1 != 0:
                    dict_ref_gaps[ref_name][c_var.start] = (diff_hap0, diff_hap1)

            for idx, c_var in enumerate(cohort_vars):
                start0 = list_start_hap[0][idx]
                start1 = list_start_hap[1][idx]
                seq_0 = seq_hap0[start0 - padding:start0 + list_len_hap[0][idx] + padding]
                seq_1 = seq_hap1[start1 - padding:start1 + list_len_hap[1][idx] + padding]
                if dict_ref_haps[ref_name].get((c_var.start)):
                    print("WARNNING! Duplicate variant at contig:", var.contig, ",pos:", c_var.start)
                dict_ref_haps[ref_name][(c_var.start)] = (seq_0, seq_1)
            if not conflict_flag: # only generate the cohort if there are no conflict alleles
                seq_hap0 = seq_hap0[var_chain-padding:start0 + list_len_hap[0][idx] + padding]
                seq_hap1 = seq_hap1[var_chain-padding:start1 + list_len_hap[1][idx] + padding]
                max_cohort_stop = cohort_vars[-1].stop
                for idy, c_var in enumerate(cohort_vars): # the start and end position of each var in the cohort
                    lpad_0 = list_start_hap[0][idy] - (var_chain-padding)
                    lpad_1 = list_start_hap[1][idy] - (var_chain-padding)
                    rpad_0 = len(seq_hap0) - lpad_0 - list_len_hap[0][idy]
                    rpad_1 = len(seq_hap1) - lpad_1 - list_len_hap[1][idy]
                    dict_ref_cohorts[ref_name][(c_var.start)] = (var.start, max_cohort_stop, seq_hap0, seq_hap1, lpad_0, lpad_1, rpad_0, rpad_1) # c_var should be the last in cohort
            idx_vcf += len(cohort_vars) # While Loop Management
        else: # single variant
            var_start = var.start - padding
            var_stop  = var.stop  + padding
            ref_seq = f_fasta.fetch(reference=var.contig, start= var_start, end = var_stop)
            hap_0, hap_1 = var.samples[0]['GT']
            seq_hap0,diff_hap0,_ = switch_var_seq(var, ref_seq, var_start, hap_0)
            seq_hap1,diff_hap1,_ = switch_var_seq(var, ref_seq, var_start, hap_1)
            if dict_ref_haps[ref_name].get((var.start)):
                print("WARNNING! Duplicate variant at contig:", var.contig, ",pos:", var.start)
            dict_ref_haps[ref_name][(var.start)] = (seq_hap0, seq_hap1)
            if diff_hap0 != 0 or diff_hap1 != 0:
                dict_ref_gaps[ref_name][var.start] = (diff_hap0, diff_hap1)
            idx_vcf += 1 # While Loop Management
        
    return dict_set_conflict_vars, dict_ref_haps, dict_ref_cohorts, dict_ref_gaps




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-r', '--real_data', help='turn off hap_information warning for real data', action='store_true')
    parser.add_argument('-p', '--golden_pickle', help='the pickle file contain the golden information for report reference')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    
    fn_vcf = args.vcf
    fn_sam = args.sam
    fn_fasta = args.fasta
    flag_real = args.real_data
    fn_golden = args.golden_pickle
    fn_output = args.out
    
    f_vcf   = pysam.VariantFile(fn_vcf)
    f_sam   = pysam.AlignmentFile(fn_sam)
    f_fasta = pysam.FastaFile(fn_fasta)
    #var_chain = 15
    padding = 5
    var_chain = 25
    #padding   = 10
    print("Start building the variant maps...")
    dict_set_conflict_vars, dict_ref_haps, dict_ref_cohorts, dict_ref_gaps = variant_seq(
            f_vcf=f_vcf,
            f_fasta=f_fasta,
            var_chain=var_chain,
            padding=padding)
    # extend conflict set
    for ref_name in dict_set_conflict_vars.keys():
        for pos in list(dict_set_conflict_vars[ref_name]):
            for extend in range(pos-var_chain, pos+var_chain):
                dict_set_conflict_vars[ref_name].add(extend)
    
    print("Start comparing reads to the variant map...")
    dict_ref_bias = compare_sam_to_haps(
            f_vcf=f_vcf,
            f_sam=f_sam,
            dict_ref_haps=dict_ref_haps,
            dict_ref_gaps=dict_ref_gaps,
            dict_ref_cohorts=dict_ref_cohorts,
            dict_set_conflict_vars=dict_set_conflict_vars,
            flag_real=flag_real,
            fn_golden=fn_golden,
            padding=padding)
    
    f_vcf   = pysam.VariantFile(fn_vcf)
    print("Start output report...")
    output_report(
            f_vcf=f_vcf,
            dict_ref_bias=dict_ref_bias,
            dict_set_conflict_vars=dict_set_conflict_vars,
            flag_real=flag_real,
            fn_golden=fn_golden,
            fn_output=fn_output)


