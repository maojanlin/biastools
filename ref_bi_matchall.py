import argparse
import re
from indelpost import Variant, VariantAlignment
import pickle
import os.path
from os import path
import pysam
import numpy as np


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
    f_all.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\tGAP\n")
    f_gap.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\n")
    f_SNP.write("CHR\tHET_SITE\tREFERENCE_BIAS\tREF_COUNT\tALT_COUNT\tOTHER_COUNT\tNUM_READS\tSUM_MAPQ\tREAD_DISTRIBUTION\n")
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
        output_string = (ref_name + '\t' + str(var.start+1) + '\t')
        if sum(n_var[:2]) == 0:
            output_string += ("N/A")
        else:
            output_string += (format(n_var[idx_ref] / float(sum(n_var[:2])), '.8f'))
        output_string += ("\t" + str(n_var[idx_ref]) + "\t" + str(n_var[idx_alt]) + "\t" + str(n_var[2]) + "\t" + str(sum(n_read)) + "\t" + str(sum(map_q)) + "\t")
        if sum(n_read) == 0:
            output_string += ("N/A")
        else:
            output_string += (format(n_read[idx_ref] / float(sum(n_read[:2])), '.8f'))
        
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


def match_to_hap(
        read_start  :int,
        var_start   :int,
        seq_read    :str,
        seq_hap     :str,
        cigar_tuples:tuple,
        padding     :int
        ) -> bool:
    """
    Adjust and compare the two sequences
    """
    if read_start > var_start:
        return False

    ref_curser  = read_start
    read_curser = 0
    for pair_info in cigar_tuples:
        code, runs = pair_info
        if code == 0 or code == 7 or code == 8: # M or = or X
            ref_curser += runs
            if ref_curser > var_start:
                #TODO
                read_curser += (runs - ref_curser + var_start)
                break
            else:
                read_curser += runs
        elif code == 1: # I
            ref_curser  += 1
            if ref_curser > var_start:
                break
            else:
                read_curser += runs
        elif code == 2: # D
            ref_curser += runs
            if ref_curser > var_start:
                break
            else:
                read_curser += 1
        elif code == 4 or code == 5: # S or H, pysam already parsed
            pass
        else:
            print ("ERROR: unexpected cigar code in sequence")

    r_start = read_curser
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


def compare_sam_to_haps(
    f_vcf           :pysam.VariantFile,
    f_sam           :pysam.AlignmentFile,
    dict_ref_haps   :dict,
    padding         :int=5
    ) -> dict:
    """
    Input:  f_sam file
    Output: ref bias dictionary according to variants
    """
    # build up the ref bias dictionary
    dict_ref_var_bias = {}
    for ref_name in dict_ref_haps.keys():
        dict_ref_var_bias[ref_name] = {}
        for start_pos in dict_ref_haps[ref_name]:
            # n_var has hap0, hap1, and others
            dict_ref_var_bias[ref_name][start_pos] = {'n_read':[0,0], 'n_var':[0,0,0], 'map_q':[0,0]}
    
    # scanning all the read alignments
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
        """
        for idx, var in enumerate(related_vars): # make sure the variants are totally contained in the read
            if var.start < pos_start or var.stop > pos_end-1:
                related_vars.pop(idx)
        if related_vars == []:
            continue
        """
        #dict_read_map = map_read_to_ref(read_start=pos_start, read_end=pos_end,cigar_tuples=cigar_tuples)
        for var in related_vars:
            seq_hap0, seq_hap1 = dict_ref_haps[ref_name][var.start]

            #print("==============================")
            #fetching the sequence in the read_seq regarding to the variant
            match_flag = False
            if var.start >= pos_start:
                #if match_hap(var.start, dict_read_map, read_seq, seq_hap0, padding):
                if match_to_hap(pos_start, var.start, read_seq, seq_hap0, cigar_tuples, padding):
                    dict_ref_var_bias[ref_name][var.start]['n_var'][0] += 1
                    match_flag = True
                #if match_hap(var.start, dict_read_map, read_seq, seq_hap1, padding):
                if match_to_hap(pos_start, var.start, read_seq, seq_hap1, cigar_tuples, padding):
                    dict_ref_var_bias[ref_name][var.start]['n_var'][1] += 1
                    match_flag = True
            if not match_flag:
                if hap_inside(read_seq, seq_hap0, padding):
                    dict_ref_var_bias[ref_name][var.start]['n_var'][0] += 1
                    match_flag = True
                if hap_inside(read_seq, seq_hap1, padding):
                    dict_ref_var_bias[ref_name][var.start]['n_var'][1] += 1
                    match_flag = True
                if match_flag == False:
                    dict_ref_var_bias[ref_name][var.start]['n_var'][2] += 1
            
            # standard updating of read number and mapping quality
            if 'hapA' == rg_tag:
                dict_ref_var_bias[ref_name][var.start]['n_read'][0] += 1
                dict_ref_var_bias[ref_name][var.start]['map_q'][0]  += mapq
            elif 'hapB' == rg_tag:
                dict_ref_var_bias[ref_name][var.start]['n_read'][1] += 1
                dict_ref_var_bias[ref_name][var.start]['map_q'][1]  += mapq
            else:
                print("WARNING, there is a read without haplotype information!!")
        """
        print(dict_ref_var_bias['chr21'][14238188])
        print(dict_ref_var_bias['chr21'][14238200])
        print(dict_ref_var_bias['chr21'][14238206])
        """

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
        return ref, 0
    else:
        alt = var.alts[genotype - 1]
        return ref[:var.start-start] + alt + ref[var.stop-start:], len(var.ref) - len(alt)


def variant_seq(
        f_vcf   :pysam.VariantFile,
        f_fasta :pysam.FastaFile,
        padding :int=5
        )-> tuple: # dict_set_conflict_vars, dict_var_haps
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
    dict_ref_haps = {}
    dict_set_conflict_vars = {}
    for ref_name in f_fasta.references:
        dict_ref_haps[ref_name] = {}
        dict_set_conflict_vars[ref_name] = set()

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
                    dict_set_conflict_vars[ref_name].add(prev_start0)
                    dict_set_conflict_vars[ref_name].add(c_var.start)
                if c_var.start > prev_start1 + diff_hap1:
                    adj_hap1 += diff_hap1
                    seq_hap1, diff_hap1 = switch_var_seq(c_var, seq_hap1, adj_hap1, hap_1)
                    prev_start1 = c_var.start
                    if c_var.start < var.start:
                        l_diff1 -= diff_hap1
                    elif c_var.start != var.start:
                        r_diff1 -= diff_hap1
                else:
                    dict_set_conflict_vars[ref_name].add(prev_start1)
                    dict_set_conflict_vars[ref_name].add(c_var.start)

            seq_hap0 = seq_hap0[max(0,l_diff0):len(seq_hap0) - max(0,r_diff0)]
            seq_hap1 = seq_hap1[max(0,l_diff1):len(seq_hap1) - max(0,r_diff1)]
        else: # single variant
            ref_seq = f_fasta.fetch(reference=var.contig, start= var_start, end = var_stop)
            hap_0, hap_1 = var.samples[0]['GT']
            seq_hap0,_ = switch_var_seq(var, ref_seq, var_start, hap_0)
            seq_hap1,_ = switch_var_seq(var, ref_seq, var_start, hap_1)

        if dict_ref_haps[ref_name].get((var.start)):
            print("WARNNING! Duplicate variant at contig:", var.contig, ",pos:", var.start)
        
        ## FOR DEBUG!!!!!
        #if var.start > 14239000:
        #    break
        dict_ref_haps[ref_name][(var.start)] = (seq_hap0, seq_hap1)
    return dict_set_conflict_vars, dict_ref_haps




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
    padding = 5
    print("Start building the variant maps...")
    dict_set_conflict_vars, dict_ref_haps = variant_seq(
            f_vcf=f_vcf,
            f_fasta=f_fasta,
            padding=padding)
    # extend conflict set
    for ref_name in dict_set_conflict_vars.keys():
        for pos in list(dict_set_conflict_vars[ref_name]):
            for extend in range(pos-padding, pos+padding):
                dict_set_conflict_vars[ref_name].add(extend)
    """
    for key, value in dict_ref_haps['chr21'].items():
        #print(key, value)
        if key in dict_set_conflict_vars['chr21']:
            continue
        print('>' + str(key))
        print(value[0])
    """
    print("Start comparing reads to the variant map...")
    dict_ref_bias = compare_sam_to_haps(
            f_vcf=f_vcf,
            f_sam=f_sam,
            dict_ref_haps=dict_ref_haps,
            padding=padding)
    """
    for key, value in dict_ref_bias['chr21'].items():
        if value['n_read'][0] == 0:
            continue
        print(key, value)
    """
    f_vcf   = pysam.VariantFile(fn_vcf)
    print("Start output report...")
    output_report(
            f_vcf=f_vcf,
            dict_ref_bias=dict_ref_bias,
            dict_set_conflict_vars=dict_set_conflict_vars, 
            fn_output=fn_output)


