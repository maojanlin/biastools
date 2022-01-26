import argparse
import re
from indelpost import Variant, VariantAlignment
import pickle
import os.path
from os import path
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
            dict_ref_alts[ref_name][var.start] = [var_seq0, var_seq1]
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


def check_coordinate(
        dict_ref_alts   :dict,
        f_hap0_fasta    :pysam.FastaFile,
        f_hap1_fasta    :pysam.FastaFile,
        dict_ref_consensus_map0: dict,
        dict_ref_consensus_map1: dict,
        dict_set_conflict_vars:  dict
        ) -> None:
    """
    Make sure the mapping point result in the same sequence as shown in the vcf file
    """
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
            fetch_hap_0 = f_hap0_fasta.fetch(reference=ref_name, start=pos_map0, end=pos_map0 + len(seq_hap0))
            fetch_hap_1 = f_hap1_fasta.fetch(reference=ref_name, start=pos_map1, end=pos_map1 + len(seq_hap1))
            if seq_hap0 != fetch_hap_0:
                print("Discrepency at", str(var_start), str(pos_map0), "haplotype 0! Expect", seq_hap0, ", get", fetch_hap_0, "...")
                count_discrepency += 1
            if seq_hap1 != fetch_hap_1:
                print("Discrepency at", str(var_start), str(pos_map1), "haplotype 1! Expect", seq_hap1, ", get", fetch_hap_1, "...")
                count_discrepency += 1
    print("Total Discrepency:", count_discrepency)


def variant_map(
        fn_chain                :str,
        dict_ref_alts           :dict,
        dict_set_conflict_vars  :dict
        ) -> tuple:
    """
    Using chain file to build the variant map
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
        dict_set_conflict_vars  :dict
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
            hap0_start = dict_ref_consensus_map0[ref_name][var_start]
            hap0_stop  = hap0_start + len(hap_seqs[0])
            hap1_start = dict_ref_consensus_map1[ref_name][var_start]
            hap1_stop  = hap1_start + len(hap_seqs[1])
            
            # read numbers overlapping the variants
            count0 = f_sam0.count(contig=ref_name, start=hap0_start, stop=hap0_stop)
            count1 = f_sam1.count(contig=ref_name, start=hap1_start, stop=hap1_stop)
            dict_ref_var_count[ref_name][var_start] = (count0,count1)
    return dict_ref_var_count


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
                read_distribution = count0/(count0+count1)
                distring = format(read_distribution, '.8f') + '\t' + str(count0) + '\t' + str(count1)
                len_var = len(var.alts[hap_1-1])
            else:
                read_distribution = count1/(count0+count1)
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
    padding = 5
    
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
            for extend in range(pos-padding, pos+padding):
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
    check_coordinate(
            dict_ref_alts=dict_ref_alts,
            f_hap0_fasta=f_hap0_fasta,
            f_hap1_fasta=f_hap1_fasta,
            dict_ref_consensus_map0=dict_ref_consensus_map0,
            dict_ref_consensus_map1=dict_ref_consensus_map1,
            dict_set_conflict_vars=dict_set_conflict_vars
            )
    print("Checking the simulation sam file covering of the variants")
    f_sam0 = pysam.AlignmentFile(fn_sam0)
    f_sam1 = pysam.AlignmentFile(fn_sam1)
    dict_ref_var_count = count_haps(
            dict_ref_alts=dict_ref_alts,
            f_sam0=f_sam0,
            f_sam1=f_sam1,
            dict_ref_consensus_map0=dict_ref_consensus_map0,
            dict_ref_consensus_map1=dict_ref_consensus_map1,
            dict_set_conflict_vars=dict_set_conflict_vars
            )
    f_vcf = pysam.VariantFile(fn_vcf)
    print("Start output report...")
    output_report(
            f_vcf=f_vcf,
            dict_ref_var_count=dict_ref_var_count,
            fn_output=fn_output)
