import argparse
import pickle
import pysam
import numpy as np



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
            dict_ref_alts[ref_name][var.start] = [var_seq0, var_seq1, hap_0, hap_1]
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


def hap_seq(
        var     :pysam.VariantRecord,
        genotype:int
        )-> str :
    """
    return variant sequence according to haplotype information
    """
    if genotype == 0:
        return var.ref
    else:
        return var.alts[genotype - 1]


def left_right_check(seq_hap0, seq_hap1):
    """
    Check the extension direction of the repetitiveness
    return:
        - 0: right side extension
        - 1: left side extension
        - 2: both sides are extensible
    """
    assert(seq_hap0 != seq_hap1)
    assert((seq_hap0 in seq_hap1) or (seq_hap1 in seq_hap0))
    len_0 = len(seq_hap0)
    len_1 = len(seq_hap1)
    if len_0 > len_1:
        if seq_hap0[:len_1] == seq_hap1:
            return 0 # right side repetitive
        elif seq_hap0[-len_1:] == seq_hap1:
            return 1 # left side repetitive
    else:
        if seq_hap1[:len_0] == seq_hap0:
            return 0 # right side repetitive
        elif seq_hap1[-len_0:] == seq_hap0:
            return 1 # left side repetitive
    return 2 # in the middle


def extend_ref_seq(
        seq_hap0,
        seq_hap1,
        ref_extend_0,
        ref_extend_1,
        flag_right=True
        )-> tuple:
    """
    Extend the seq_hap0 and seq_hap1 till they makes a difference
    """
    seq_hap0_extend = seq_hap0
    seq_hap1_extend = seq_hap1
    assert((seq_hap0_extend in seq_hap1_extend) or (seq_hap1_extend in seq_hap0_extend))
    len_iterate = min(len(ref_extend_0), len(ref_extend_1))
    if flag_right: # extend to the right
        for idx in range(len_iterate):
            seq_hap0_extend += ref_extend_0[idx]
            seq_hap1_extend += ref_extend_1[idx]
            if (seq_hap0_extend in seq_hap1_extend) or (seq_hap1_extend in seq_hap0_extend): # still indistinguishable
                continue
            else:
                return seq_hap0_extend, seq_hap1_extend, idx+1
    else: # extend to the left
        for idx in range(len_iterate):
            seq_hap0_extend = ref_extend_0[-idx-1] + seq_hap0_extend
            seq_hap1_extend = ref_extend_1[-idx-1] + seq_hap1_extend
            if (seq_hap0_extend in seq_hap1_extend) or (seq_hap1_extend in seq_hap0_extend): # still indistinguishable
                continue
            else:
                return seq_hap0_extend, seq_hap1_extend, idx+1
    return seq_hap0_extend, seq_hap1_extend, False


def extend_ref_seq_padding(
        seq_hap0,
        seq_hap1,
        ref_extend_0,
        ref_extend_1,
        flag_right=True,
        padding=5
        ):
    """
    Call the extend_ref_seq and add padding in the end
    """
    if flag_right:
        seq_hap0_extend, seq_hap1_extend, len_extend = extend_ref_seq(seq_hap0, seq_hap1, ref_extend_0[:-padding], ref_extend_1[:-padding], flag_right)
        if len_extend:
            return seq_hap0_extend + ref_extend_0[len_extend:len_extend+padding], seq_hap1_extend + ref_extend_1[len_extend:len_extend+padding], len_extend+padding
        else:
            return seq_hap0, seq_hap1, False
    else:
        seq_hap0_extend, seq_hap1_extend, len_extend = extend_ref_seq(seq_hap0, seq_hap1, ref_extend_0[padding:], ref_extend_1[padding:], flag_right)
        if len_extend:
            return ref_extend_0[-len_extend-padding:-len_extend] + seq_hap0_extend, ref_extend_1[-len_extend-padding:-len_extend] + seq_hap1_extend, len_extend+padding
        else:
            return seq_hap0, seq_hap1, False



def nearest_left_right_var(
        left_0, 
        right_0, 
        f_hap0_fasta,
        left_1,
        right_1,
        f_hap1_fasta,
        ref_name,
        left_extend=40,
        right_extend=40
        ) -> tuple:
    left_seq_0 = f_hap0_fasta.fetch(reference=ref_name, start=left_0 - left_extend, end=left_0)
    left_seq_1 = f_hap1_fasta.fetch(reference=ref_name, start=left_1 - left_extend, end=left_1)
    left_var = -1
    for idx in range(left_extend-1, 0, -1):
        if left_seq_0[idx] != left_seq_1[idx]:
            left_var = left_extend-idx
            break
    right_seq_0 = f_hap0_fasta.fetch(reference=ref_name, start=right_0, end=right_0 + right_extend)
    right_seq_1 = f_hap1_fasta.fetch(reference=ref_name, start=right_1, end=right_1 + right_extend)
    right_var = -1
    for idx in range(right_extend):
        if right_seq_0[idx] != right_seq_1[idx]:
            right_var = idx
            break
    return left_var, right_var



def check_coordinate(
        dict_ref_alts   :dict,
        f_hap0_fasta    :pysam.FastaFile,
        f_hap1_fasta    :pysam.FastaFile,
        dict_ref_consensus_map0: dict,
        dict_ref_consensus_map1: dict,
        dict_set_conflict_vars:  dict,
        extend_limit    :int=100,
        padding         :int=5
        ) -> dict:
    """
    Make sure the mapping point result in the same sequence as shown in the vcf file
    dict_effective_variant {}
        - key: var_start (at reference coordinate)
        - values: [flag_side, len_extend] # len_extend can be either right or left
    """
    dict_effective_variant = {}
    count_discrepency = 0
    for ref_name, dict_var_seq in dict_ref_alts.items():
        set_conflict = dict_set_conflict_vars[ref_name]
        for var_start, pair_var_seq in dict_var_seq.items():
            if var_start in set_conflict:
                continue
        
            seq_hap0 = pair_var_seq[0]
            seq_hap1 = pair_var_seq[1]
            pos_map0 = dict_ref_consensus_map0[ref_name][var_start]
            pos_map1 = dict_ref_consensus_map1[ref_name][var_start]
            long_hap0 = f_hap0_fasta.fetch(reference=ref_name, start=pos_map0-padding, end=pos_map0 + len(seq_hap0)+padding)
            long_hap1 = f_hap1_fasta.fetch(reference=ref_name, start=pos_map1-padding, end=pos_map1 + len(seq_hap1)+padding)

            if long_hap0 != long_hap1:
                if (long_hap0 in long_hap1) or (long_hap1 in long_hap0):
                    flag_side = left_right_check(long_hap0, long_hap1) # check which side the repetitive be
                    if flag_side == 0: # right side
                        # get additional extend_limit (default 100) bp from the reference
                        extend_hap0 = f_hap0_fasta.fetch(reference=ref_name, start=pos_map0+len(seq_hap0)+padding, end=pos_map0+len(seq_hap0)+extend_limit)
                        extend_hap1 = f_hap1_fasta.fetch(reference=ref_name, start=pos_map1+len(seq_hap1)+padding, end=pos_map1+len(seq_hap1)+extend_limit)
                        effect_hap0, effect_hap1, len_extend = extend_ref_seq_padding(long_hap0, long_hap1, extend_hap0, extend_hap1, True, padding)
                        if len_extend:
                            left_var, right_var = nearest_left_right_var(pos_map0, pos_map0+len(seq_hap0), f_hap0_fasta, \
                                                                         pos_map1, pos_map1+len(seq_hap1), f_hap1_fasta, ref_name, 40, 40+len_extend)
                            dict_effective_variant[var_start] = (0, len_extend, left_var, right_var)
                        else:
                            print("--- 0 EFFECTIVE VARIANT too long at", var_start, seq_hap0, seq_hap1)
                    elif flag_side == 1: # left side
                        extend_hap0 = f_hap0_fasta.fetch(reference=ref_name, start=pos_map0-extend_limit-padding, end=pos_map0-padding)
                        extend_hap1 = f_hap1_fasta.fetch(reference=ref_name, start=pos_map1-extend_limit-padding, end=pos_map1-padding)
                        effect_hap0, effect_hap1, len_extend = extend_ref_seq_padding(long_hap0, long_hap1, extend_hap0, extend_hap1, False, padding)
                        if len_extend:
                            left_var, right_var = nearest_left_right_var(pos_map0, pos_map0+len(seq_hap0), f_hap0_fasta, \
                                                                         pos_map1, pos_map1+len(seq_hap1), f_hap1_fasta, ref_name, 40+len_extend, 40)
                            dict_effective_variant[var_start] = (1, len_extend, left_var, right_var)
                        else:
                            print("--- 1 EFFECTIVE VARIANT too long at", var_start, seq_hap0, seq_hap1)
                    else: # both sides are extensible
                        extend_hap0 = f_hap0_fasta.fetch(reference=ref_name, start=pos_map0+len(seq_hap0)+padding, end=pos_map0+len(seq_hap0)+extend_limit)
                        extend_hap1 = f_hap1_fasta.fetch(reference=ref_name, start=pos_map1+len(seq_hap1)+padding, end=pos_map1+len(seq_hap1)+extend_limit)
                        r_effect_hap0, r_effect_hap1, r_len_extend = extend_ref_seq_padding(long_hap0, long_hap1, extend_hap0, extend_hap1, True, padding)
                        
                        extend_hap0 = f_hap0_fasta.fetch(reference=ref_name, start=pos_map0-extend_limit-padding, end=pos_map0-padding)
                        extend_hap1 = f_hap1_fasta.fetch(reference=ref_name, start=pos_map1-extend_limit-padding, end=pos_map1-padding)
                        l_effect_hap0, l_effect_hap1, l_len_extend = extend_ref_seq_padding(long_hap0, long_hap1, extend_hap0, extend_hap1, False, padding)
                        flag_extend = -1
                        if l_len_extend == 0: # right anyway
                            if r_len_extend == 0:
                                print("--- 2 EFFECTIVE VARIANT ENCOUNTER at", var_start, seq_hap0, seq_hap1, "L", l_len_extend)
                            else:
                                flag_extend=0
                        elif r_len_extend == 0: # left anyway
                            flag_extend=1
                        elif r_len_extend < l_len_extend: # right is better
                            flag_extend=0
                        else: # left is better
                            flag_extend=1

                        if flag_extend == 0:
                            left_var, right_var = nearest_left_right_var(pos_map0, pos_map0+len(seq_hap0), f_hap0_fasta, \
                                                                         pos_map1, pos_map1+len(seq_hap1), f_hap1_fasta, ref_name, 40, 40+r_len_extend)
                            dict_effective_variant[var_start] = (0, r_len_extend, left_var, right_var)
                        elif flag_extend == 1:
                            left_var, right_var = nearest_left_right_var(pos_map0, pos_map0+len(seq_hap0), f_hap0_fasta, \
                                                                         pos_map1, pos_map1+len(seq_hap1), f_hap1_fasta, ref_name, 40+l_len_extend, 40)
                            dict_effective_variant[var_start] = (1, l_len_extend, left_var, right_var)

            fetch_hap_0 = long_hap0[padding:-padding] 
            fetch_hap_1 = long_hap1[padding:-padding]
            if seq_hap0.upper() != fetch_hap_0.upper() and seq_hap0 != '*':
                print("Discrepency at", ref_name, str(var_start), str(pos_map0), "haplotype 0! Expect", seq_hap0, ", get", fetch_hap_0, "...")
                count_discrepency += 1
            if seq_hap1.upper() != fetch_hap_1.upper() and seq_hap1 != '*':
                print("Discrepency at", ref_name, str(var_start), str(pos_map1), "haplotype 1! Expect", seq_hap1, ", get", fetch_hap_1, "...")
                count_discrepency += 1
    print("Total Discrepency:", count_discrepency)
    return dict_effective_variant


def variant_map(
        fn_chain                :str,
        dict_ref_alts           :dict,
        dict_set_conflict_vars  :dict
        ) -> tuple:
    """
    Using chain file to build the variant map
    mapping from reference to target genome coordinate
    """
    dict_ref_consensus_map = {}
    for ref_name in dict_ref_alts.keys():
        dict_ref_consensus_map[ref_name] = {}
    
    # Read and parse the chain file
    dict_chain_info = {}
    key_tuple = None
    fc = open(fn_chain, 'r')
    for line in fc:
        fields = line.strip().split()
        if len(fields) > 0 and fields[0] == "chain":
            key_tuple = tuple(fields)
            dict_chain_info[key_tuple] = []
        else:
            dict_chain_info[key_tuple].append(fields)
    fc.close()
    
    for key_tuple, list_info in dict_chain_info.items():
        assert(key_tuple[4] == key_tuple[9])
        t_start  = int(key_tuple[5])
        assert(t_start == 0) # in this version, we only support one whole genome
        t_stop   = int(key_tuple[6])
        h_start  = int(key_tuple[10])
        ref_name = key_tuple[2]
        
        list_var_start  = sorted(dict_ref_alts[ref_name].keys())
        set_conflict = dict_set_conflict_vars[ref_name]
        idx_chain = 0
        pos_chain = t_start + int(list_info[idx_chain][0])
        offset = 0
        for var_start in list_var_start:
            if var_start in set_conflict:
                continue
            elif var_start < t_start:
                continue
            elif var_start > t_stop:
                break

            if var_start < pos_chain:
                dict_ref_consensus_map[ref_name][var_start] = var_start + offset
            else:
                while pos_chain <= var_start:
                    pos_chain += int(list_info[idx_chain][1])
                    offset -= int(list_info[idx_chain][1])
                    offset += int(list_info[idx_chain][2])
                    idx_chain += 1
                    pos_chain += int(list_info[idx_chain][0])
                dict_ref_consensus_map[ref_name][var_start] = var_start + offset
    return dict_ref_consensus_map
            

def count_haps(
        dict_ref_alts   :dict,
        f_sam0          :pysam.AlignmentFile,
        f_sam1          :pysam.AlignmentFile,
        dict_ref_consensus_map0 :dict,
        dict_ref_consensus_map1 :dict,
        dict_set_conflict_vars  :dict,
        debug           :bool=False
        ) -> dict:
    """
    Count the number of reads in each golden haplotype sam covering the variants
    """
    dict_ref_var_count = {}
    for ref_name, dict_vars in dict_ref_alts.items():
        dict_ref_var_count[ref_name] = {}
        set_conflict = dict_set_conflict_vars[ref_name]
        for var_start, hap_seqs in dict_vars.items():
            if var_start in set_conflict:
                continue
            if hap_seqs[2] == hap_seqs[3]: # if the var is homozygous
                continue
            hap0_start = dict_ref_consensus_map0[ref_name][var_start]
            hap0_stop  = hap0_start + len(hap_seqs[0])
            hap1_start = dict_ref_consensus_map1[ref_name][var_start]
            hap1_stop  = hap1_start + len(hap_seqs[1])
            
            # read numbers overlapping the variants
            count0 = f_sam0.count(contig=ref_name, start=hap0_start, stop=hap0_stop)
            count1 = f_sam1.count(contig=ref_name, start=hap1_start, stop=hap1_stop)
            if debug:
                print(ref_name, var_start, ':\n\thapA (' + str(count0) + "): ", end="")
                for read in f_sam0.fetch(contig=ref_name, start=hap0_start, stop=hap0_stop):
                    print(read.query_name, end=", ")
                print("\n\thapB (" + str(count1) + "): ", end="")
                for read in f_sam1.fetch(contig=ref_name, start=hap1_start, stop=hap1_stop):
                    print(read.query_name, end=", ")
                print("\n", end="")
                
            dict_ref_var_count[ref_name][var_start] = (count0,count1)
    return dict_ref_var_count


def get_bound(
        hap0_start,
        hap0_stop,
        hap1_start,
        hap1_stop,
        len_hap0,
        len_hap1,
        flag_side,
        len_extend
    ) -> tuple:
    """
    return: (eff_lbound0, eff_rbound0, eff_lbound1, eff_rbound1)
    """
    min_len = min(len_hap0, len_hap1)
    if flag_side == 0: # extend to the right
        if len_hap0 < len_hap1:
            return (hap0_stop+len_extend, \
                    hap0_start, \
                    hap1_start+min_len+len_extend, \
                    hap1_stop) #min(hap1_stop, hap1_start+min_len+len_extend))
        else:
            return (hap0_start+min_len+len_extend, \
                    hap0_stop, \
                    hap1_stop+len_extend, \
                    hap1_start)
                    #min(hap0_stop, hap0_start+min_len+len_extend), \
    else: # extend to the left
        if len_hap0 < len_hap1:
            return (hap0_stop, \
                    hap0_start-len_extend, \
                    hap1_start, \
                    hap1_stop-min_len-len_extend)
                    #max(hap1_start, hap1_stop-min_len-len_extend), \
        else:
            return (hap0_start, \
                    hap0_stop-min_len-len_extend, \
                    hap1_stop, \
                    hap1_start-len_extend)
                    #max(hap0_start, hap0_stop-min_len-len_extend), \



def count_haps_n_report_name(
        dict_ref_alts   :dict,
        f_sam0          :pysam.AlignmentFile,
        f_sam1          :pysam.AlignmentFile,
        dict_ref_consensus_map0 :dict,
        dict_ref_consensus_map1 :dict,
        dict_set_conflict_vars  :dict,
        dict_effective_var      :dict,
        padding         :int=5,
        debug           :bool=False
        ) -> dict:
    """
    Count the number of reads in each golden haplotype sam covering the variants
    """
    dict_ref_var_count = {}
    dict_ref_var_name  = {}
    for ref_name, dict_vars in dict_ref_alts.items():
        dict_ref_var_count[ref_name] = {}
        dict_ref_var_name [ref_name] = {}
        set_conflict = dict_set_conflict_vars[ref_name]

        len_dict_vars = len(dict_vars)
        for idx, (var_start, hap_seqs) in enumerate(dict_vars.items()):
            if var_start in set_conflict:
                continue
            if hap_seqs[2] == hap_seqs[3]: # if the var is homozygous
                continue
            hap0_start = dict_ref_consensus_map0[ref_name][var_start]
            hap0_stop  = hap0_start + len(hap_seqs[0])
            hap1_start = dict_ref_consensus_map1[ref_name][var_start]
            hap1_stop  = hap1_start + len(hap_seqs[1])
            
            #if var_start < 6611800:
            #    continue
            # read numbers overlapping the variants
            if dict_effective_var.get(var_start): # if the site has larger effective var size
                flag_side, len_extend, left_var, right_var = dict_effective_var[var_start]
                min_len = min(len(hap_seqs[0]), len(hap_seqs[1]))
                if flag_side == 0: # right extend
                    eff_start0 = hap0_start
                    eff_stop0  = hap0_start + len(hap_seqs[0]) + len_extend
                    eff_start1 = hap1_start
                    eff_stop1  = hap1_start + len(hap_seqs[1]) + len_extend
                else:
                    eff_start0 = hap0_stop - len(hap_seqs[0]) - len_extend
                    eff_stop0  = hap0_stop
                    eff_start1 = hap1_stop - len(hap_seqs[1]) - len_extend
                    eff_stop1  = hap1_stop
                eff_lbound0, eff_rbound0, eff_lbound1, eff_rbound1 = get_bound(hap0_start, hap0_stop, hap1_start, hap1_stop, \
                                                                               len(hap_seqs[0]), len(hap_seqs[1]), flag_side, len_extend)
                # compensate for the nearby variants
                if left_var != -1:
                    eff_lbound0 = min(eff_lbound0, hap0_start-left_var)
                    eff_lbound1 = min(eff_lbound1, hap1_start-left_var)
                if right_var != -1:
                    eff_rbound0 = max(eff_rbound0, hap0_stop+right_var)
                    eff_rbound1 = max(eff_rbound1, hap1_stop+right_var)
                
                read_segment0 = f_sam0.fetch(contig=ref_name, start=eff_start0, stop=eff_stop0)
                set_expand0 = set()
                set_inside0 = set()
                for read in read_segment0:
                    if read.reference_end >= eff_lbound0 and read.reference_start <= eff_rbound0:
                        set_expand0.add((read.query_name, read.is_read2))
                    else:
                        set_inside0.add((read.query_name, read.is_read2))

                read_segment1 = f_sam1.fetch(contig=ref_name, start=eff_start1, stop=eff_stop1)
                set_expand1 = set()
                set_inside1 = set()
                #print(hap1_start, hap1_start+len(hap_seqs[0]), len_extend, left_var, right_var, flag_side)
                #print(eff_lbound1, eff_rbound1)
                for read in read_segment1:
                    if read.reference_end >= eff_lbound1 and read.reference_start <= eff_rbound1:
                        set_expand1.add((read.query_name, read.is_read2))
                    else:
                        set_inside1.add((read.query_name, read.is_read2))
                
                """
                if var_start == 6611841:
                    print(flag_side, len_extend)
                    print(hap_seqs[0], hap_seqs[1])
                    print(hap0_start, hap0_stop)
                    print(hap1_start, hap1_stop)

                    print(set_expand1)
                    print(set_inside1)
                    print(hap0_start, eff_start0, hap0_stop, eff_stop0)
                    print(hap1_start, eff_start1, hap1_stop, eff_stop1)
                    print(len(set_expand0) + len(set_expand1))
                        
                read_segment0_start = f_sam0.fetch(contig=ref_name, start=eff_start0)
                read_segment0_stop  = f_sam0.fetch(contig=ref_name, start=eff_stop0)
                read_segment1_start = f_sam1.fetch(contig=ref_name, start=eff_start1)
                read_segment1_stop  = f_sam1.fetch(contig=ref_name, start=eff_stop1)
                name_set0_start = set([read.query_name for read in read_segment0_start])
                name_set0_stop  = set([read.query_name for read in read_segment0_stop])
                name_set1_start = set([read.query_name for read in read_segment1_start])
                name_set1_stop  = set([read.query_name for read in read_segment1_stop])
                
                name_set0 = name_set0_start.intersection(name_set0_stop)
                name_set1 = name_set1_start.intersection(name_set1_stop)
                count0 = len(name_set0)
                count1 = len(name_set1)
                #symmetric_difference
                print(var_start)
                print(name_set0_start.symmetric_difference(name_set0_stop))
                print(name_set1_start.symmetric_difference(name_set1_stop))"""
                
                count0 = len(set_expand0) 
                count1 = len(set_expand1) 
                dict_ref_var_count[ref_name][var_start] = (count0, count1)
                dict_ref_var_name [ref_name][var_start] = (set_expand0, set_expand1, count0, count1, set_inside0, set_inside1)
            else:
                read_segment0 = f_sam0.fetch(contig=ref_name, start=hap0_start, stop=hap0_stop)
                name_set0    = set([(read.query_name, read.is_read2) for read in read_segment0])
                count0 = len(name_set0)
                read_segment1 = f_sam1.fetch(contig=ref_name, start=hap1_start, stop=hap1_stop)
                name_set1    = set([(read.query_name, read.is_read2) for read in read_segment1])
                count1 = len(name_set1)
                
                dict_ref_var_count[ref_name][var_start] = (count0, count1)
                dict_ref_var_name [ref_name][var_start] = (name_set0, name_set1, count0, count1)
            
            if debug:
                print(ref_name, var_start, ':\n\thapA (' + str(count0) + "): ", end="")
                for read in f_sam0.fetch(contig=ref_name, start=hap0_start, stop=hap0_stop):
                    print(read.query_name, end=", ")
                print("\n\thapB (" + str(count1) + "): ", end="")
                for read in f_sam1.fetch(contig=ref_name, start=hap1_start, stop=hap1_stop):
                    print(read.query_name, end=", ")
                print("\n", end="")
    return dict_ref_var_count, dict_ref_var_name


def output_report(
        f_vcf               :pysam.VariantFile,
        dict_ref_var_count  :dict,
        fn_output           :str
        ) -> None:
    """
    ourput report
    """
    f_all = open(fn_output, 'w')
    f_gap = open(fn_output + '.gap', 'w')
    f_SNP = open(fn_output + '.SNP', 'w')
    f_all.write("CHR\tHET_SITE\tGOLDEN_DISTRIBUTION\tREF_COUNT\tALT_COUNT\tGAP\n")
    f_gap.write("CHR\tHET_SITE\tGOLDEN_DISTRIBUTION\tREF_COUNT\tALT_COUNT\n")
    f_SNP.write("CHR\tHET_SITE\tGOLDEN_DISTRIBUTION\tREF_COUNT\tALT_COUNT\n")
    for var in f_vcf:
        hap_0, hap_1 = var.samples[0]['GT']
        if hap_0 != 0 and hap_1 != 0:
            continue
        ref_name = var.contig
        if dict_ref_var_count[ref_name].get(var.start): # Exist legal variant
            count0, count1 = dict_ref_var_count[ref_name][var.start]
            len_var = 0
            if hap_0 == 0:
                read_distribution = count0/max((count0+count1),0.001)
                distring = format(read_distribution, '.8f') + '\t' + str(count0) + '\t' + str(count1)
                len_var = len(var.alts[hap_1-1])
            else:
                read_distribution = count1/max((count0+count1),0.001)
                distring = format(read_distribution, '.8f') + '\t' + str(count1) + '\t' + str(count0)
                len_var = len(var.alts[hap_0-1])
            f_all.write(ref_name + '\t' + str(var.start+1) + '\t' + distring + '\t')
            if len(var.ref) != len_var:
                f_gap.write(ref_name + '\t' + str(var.start+1) + '\t' + distring + '\n')
                f_all.write('.\n')
            else:
                f_SNP.write(ref_name + '\t' + str(var.start+1) + '\t' + distring + '\n')
                f_all.write('\n')
    
    f_all.close()
    f_gap.close()
    f_SNP.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-c0', '--hap0_chain', help='hap0 chain file')
    parser.add_argument('-c1', '--hap1_chain', help='hap1 chain file')
    parser.add_argument('-f0', '--hap0_fasta', help='hap0 consensus fasta file')
    parser.add_argument('-f1', '--hap1_fasta', help='hap1 consensus fasta file')
    parser.add_argument('-s0', '--hap0_sam', help='hap0 sam file')
    parser.add_argument('-s1', '--hap1_sam', help='hap1 sam file')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    
    fn_vcf = args.vcf
    fn_chain0 = args.hap0_chain
    fn_chain1 = args.hap1_chain
    fn_hap0_fasta = args.hap0_fasta
    fn_hap1_fasta = args.hap1_fasta
    fn_sam0 = args.hap0_sam
    fn_sam1 = args.hap1_sam
    fn_output = args.out
    var_chain = 25
    
    f_vcf = pysam.VariantFile(fn_vcf)
    f_hap0_fasta = pysam.FastaFile(fn_hap0_fasta)
    f_hap1_fasta = pysam.FastaFile(fn_hap1_fasta)
    print("Start locating variants and the conflicting variants...")
    dict_set_conflict_vars, dict_ref_alts = variant_seq(
            f_vcf=f_vcf,
            f_fasta=f_hap0_fasta
            )
    # extend conflict set
    for ref_name in dict_set_conflict_vars.keys():
        for pos in list(dict_set_conflict_vars[ref_name]):
            for extend in range(pos-var_chain, pos+var_chain):
                dict_set_conflict_vars[ref_name].add(extend)
    print("Start building the mapping consensus coordinate...")
    dict_ref_consensus_map0 = variant_map(
            fn_chain=fn_chain0,
            dict_ref_alts=dict_ref_alts,
            dict_set_conflict_vars=dict_set_conflict_vars
            )
    dict_ref_consensus_map1 = variant_map(
            fn_chain=fn_chain1,
            dict_ref_alts=dict_ref_alts,
            dict_set_conflict_vars=dict_set_conflict_vars
            )
    # obsolete if you are confident
    print("Checking if the coordinate is correct...")
    dict_effective_var = check_coordinate(
            dict_ref_alts=dict_ref_alts,
            f_hap0_fasta=f_hap0_fasta,
            f_hap1_fasta=f_hap1_fasta,
            dict_ref_consensus_map0=dict_ref_consensus_map0,
            dict_ref_consensus_map1=dict_ref_consensus_map1,
            dict_set_conflict_vars=dict_set_conflict_vars,
            padding=10
            )
    print("Checking the simulation sam file covering of the variants")
    f_sam0 = pysam.AlignmentFile(fn_sam0)
    f_sam1 = pysam.AlignmentFile(fn_sam1)
    dict_ref_var_count, dict_ref_var_name = count_haps_n_report_name(
            dict_ref_alts=dict_ref_alts,
            f_sam0=f_sam0,
            f_sam1=f_sam1,
            dict_ref_consensus_map0=dict_ref_consensus_map0,
            dict_ref_consensus_map1=dict_ref_consensus_map1,
            dict_set_conflict_vars=dict_set_conflict_vars,
            dict_effective_var=dict_effective_var,
            padding=10,
            debug=False
            )
    f_vcf = pysam.VariantFile(fn_vcf)
    print("Start output report...")
    output_report(
            f_vcf=f_vcf,
            dict_ref_var_count=dict_ref_var_count,
            fn_output=fn_output)
    print("Dump golden read names pickle file...")
    with open(fn_output + '.pickle', 'wb') as f:
        pickle.dump(dict_ref_var_name, f)

