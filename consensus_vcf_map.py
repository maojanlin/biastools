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
        return 0, len(var.ref)
    else:
        alt = var.alts[genotype - 1]
        return len(var.ref) - len(alt), len(alt)


def variant_seq(
        f_vcf       :pysam.VariantFile,
        f_fasta     :pysam.FastaFile
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
    dict_ref_consensus_map = {}
    dict_set_conflict_vars = {}
    for ref_name in f_fasta.references:
        dict_ref_consensus_map[ref_name] = {}
        dict_set_conflict_vars[ref_name] = set()

    old_ref_name = ""
    for var in f_vcf:
        ref_name = var.contig
        if old_ref_name != ref_name: # changing the contig
            # Reset the parameters
            overlap0,    overlap1    =  0,  0
            coord_hap0,  coord_hap1  =  0,  0
            prev_start0, prev_start1 = -1, -1
            old_ref_name = ref_name
        
        dict_ref_consensus_map[ref_name][var.start] = [None, None]
        hap_0, hap_1 = var.samples[0]['GT']
        diff_hap0, len_var0 = len_var_seq(var, hap_0)
        diff_hap1, len_var1 = len_var_seq(var, hap_1)
        if var.start >= prev_start0 + overlap0 and var.start >= prev_start1 + overlap1: # checking if there are overlaps
            # hap0
            dict_ref_consensus_map[ref_name][var.start][0] = (var.start + coord_hap0, var.start + coord_hap0 + len_var0)
            coord_hap0 -= diff_hap0
            prev_start0 = var.start
            overlap0 = len_var0 - 1 if (diff_hap0 == 0) else diff_hap0
            # hap1
            dict_ref_consensus_map[ref_name][var.start][1] = (var.start + coord_hap1, var.start + coord_hap1 + len_var1)
            coord_hap1 -= diff_hap1
            prev_start1 = var.start
            overlap1 = len_var1 - 1 if (diff_hap1 == 0) else diff_hap1
        else: # overlapping variants are consider conflicts
            dict_set_conflict_vars[ref_name].add(prev_start1)
            dict_set_conflict_vars[ref_name].add(var.start)
        
    return dict_set_conflict_vars, dict_ref_consensus_map


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
        f_vcf           :pysam.VariantFile,
        f_hap0_fasta    :pysam.FastaFile,
        f_hap1_fasta    :pysam.FastaFile,
        dict_ref_consensus_map: dict,
        dict_set_conflict_vars: dict
        ) -> None:
    """
    Comment
    """
    for var in f_vcf:
        ref_name = var.contig
        if var.start in dict_set_conflict_vars[ref_name]:
            continue
        
        hap_0, hap_1 = var.samples[0]['GT']
        seq_hap0 = hap_seq(var, hap_0)
        seq_hap1 = hap_seq(var, hap_1)
        fetch_hap_0 = f_hap0_fasta.fetch(reference=ref_name, start=dict_ref_consensus_map[ref_name][var.start][0][0], end=dict_ref_consensus_map[ref_name][var.start][0][1])
        fetch_hap_1 = f_hap1_fasta.fetch(reference=ref_name, start=dict_ref_consensus_map[ref_name][var.start][1][0], end=dict_ref_consensus_map[ref_name][var.start][1][1])
        if seq_hap0 != fetch_hap_0:
            print("Discrepency at", str(var.start), "haplotype 0! Expect", seq_hap0, ", get", fetch_hap_0, "...")
        #print(var.start, var.start - dict_ref_consensus_map[ref_name][var.start][0][0]), var.start - dict_ref_consensus_map[ref_name][var.start][0][1])
        if seq_hap1 != fetch_hap_1:
            print("Discrepency at", str(var.start), "haplotype 1! Expect", seq_hap1, ", get", fetch_hap_1, "...")
        print(var.start, var.start - dict_ref_consensus_map[ref_name][var.start][0][0], var.start - dict_ref_consensus_map[ref_name][var.start][1][0])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    #parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f0', '--hap0_fasta', help='hap0 consensus fasta file')
    parser.add_argument('-f1', '--hap1_fasta', help='hap1 consensus fasta file')
    #parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    
    fn_vcf = args.vcf
    #fn_sam = args.sam
    fn_hap0_fasta = args.hap0_fasta
    fn_hap1_fasta = args.hap1_fasta
    #fn_output = args.out
    
    f_vcf   = pysam.VariantFile(fn_vcf)
    #f_sam   = pysam.AlignmentFile(fn_sam)
    f_hap0_fasta = pysam.FastaFile(fn_hap0_fasta)
    f_hap1_fasta = pysam.FastaFile(fn_hap1_fasta)
    print("Start building the mapping consensus coordinate...")
    dict_set_conflict_vars, dict_ref_consensus_map = variant_seq(
            f_vcf=f_vcf,
            f_fasta=f_hap0_fasta
            )
    # extend conflict set
    padding = 5
    for ref_name in dict_set_conflict_vars.keys():
        for pos in list(dict_set_conflict_vars[ref_name]):
            for extend in range(pos-padding, pos+padding):
                dict_set_conflict_vars[ref_name].add(extend)
    print("Checking if the coordinate is correct...")
    f_vcf   = pysam.VariantFile(fn_vcf)
    check_coordinate(
            f_vcf=f_vcf,
            f_hap0_fasta=f_hap0_fasta,
            f_hap1_fasta=f_hap1_fasta,
            dict_ref_consensus_map=dict_ref_consensus_map,
            dict_set_conflict_vars=dict_set_conflict_vars
            )
