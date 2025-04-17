import argparse
import multiprocessing
import pickle
import pysam
import numpy as np
from scipy.stats import chisquare
from typing import List, Tuple, Dict, Union
from collections import defaultdict
import logging
import edlib


# Constants
PADDING = 5
VAR_CHAIN = 25
EXTEND_LIMIT = 70

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VariantAnalyzer:
    def __init__(self, vcf_file: str, sam_file: str, fasta_file: str, 
                 golden_pickle: str = None, run_id: str = None, 
                 real_data: bool = False, output_file: str = "output.rpt", chromosome: str=None,
                 max_len_diff: int = 20):
        self.chromosome = chromosome
        self.f_vcf = pysam.VariantFile(vcf_file)
        self.f_sam = pysam.AlignmentFile(sam_file)
        self.f_fasta = pysam.FastaFile(fasta_file)
        self.golden_pickle = golden_pickle
        self.run_id = run_id
        self.real_data = real_data
        self.output_file = output_file
        self.dict_ref_haps = defaultdict(dict)
        self.dict_ref_cohorts = defaultdict(dict)
        self.dict_set_conflict_vars = defaultdict(set)
        self.dict_ref_bias = defaultdict(lambda: defaultdict(lambda: {'n_read':[0,0,0], 'n_var':[0,0,0,0], 'map_q':[0,0,0], 'distribute':[[],[],[],[]]}))
        self.dict_ref_var_pos = defaultdict(lambda: [[],[]])
        self.max_len_diff = max_len_diff
        print(f"start initializing {self.chromosome}")
        
        if not self.real_data and self.golden_pickle:
            with open(self.golden_pickle, "rb") as f:
                self.dict_ref_var_name = pickle.load(f)

    def analyze(self):
        self.build_variant_maps()
        self.compare_reads_to_variants()
        return self.generate_report()

    def build_variant_maps(self):
        logger.info("Start building the variant maps...")
        self.variant_seq()
        
        # Extend conflict set
        for ref_name, positions in self.dict_set_conflict_vars.items():
            for pos in list(positions):
                self.dict_set_conflict_vars[ref_name].update(range(pos - VAR_CHAIN, pos + VAR_CHAIN))

    def variant_seq(self):
        #for ref_name in self.f_fasta.references:

        list_var = list(self.f_vcf.fetch(self.chromosome))
        self.dict_ref_var_pos[self.chromosome] = [[var.start for var in list_var],[var.stop for var in list_var]]
        len_list_var = len(list_var)
        idx_var = 0
        while idx_var < len_list_var:
            idx_var = self.process_variant(list_var, idx_var)

    def process_variant(self, list_var: List[pysam.VariantRecord], idx_vcf: int):
        var = list_var[idx_vcf] # iterate of the variants
        ref_name = var.contig
        
        # We first treat the var as single variant, extend if it is repetitive
        var_start = var.start - PADDING
        var_stop  = var.stop  + PADDING
        ref_seq = self.f_fasta.fetch(reference=ref_name, start= var_start, end = var_stop)
        hap_0, hap_1 = var.samples[0]['GT']
        # standard making hap0 and hap1 sequence
        seq_hap0,diff_hap0,_ = self.switch_var_seq(var, ref_seq, var_start, hap_0)
        seq_hap1,diff_hap1,_ = self.switch_var_seq(var, ref_seq, var_start, hap_1)
        
        # EXTENDING the PADDING if hap0 and hap1 is indistinguishable
        if hap_0 != hap_1: # make sure it is a 0/1 haplotypes
            if (seq_hap0 in seq_hap1) or (seq_hap1 in seq_hap0):
                flag_side = self.left_right_check(seq_hap0, seq_hap1) # check which side the repetitive be
                if flag_side == 0: # right side
                    # get additional extend_limit (default 70) bp from the reference
                    ref_extend = self.f_fasta.fetch(reference=ref_name, start=var_stop, end=var_stop+EXTEND_LIMIT+PADDING)
                    seq_hap0, seq_hap1, len_extend = self.extend_ref_seq_padding(seq_hap0, seq_hap1, ref_extend, ref_extend, True, PADDING)
                    var_stop += len_extend
                    # Keep the sign of vcf_length for effective_length
                elif flag_side == 1: # left side
                    ref_extend = self.f_fasta.fetch(reference=ref_name, start=var_start-EXTEND_LIMIT-PADDING, end=var_start)
                    seq_hap0, seq_hap1, len_extend = self.extend_ref_seq_padding(seq_hap0, seq_hap1, ref_extend, ref_extend, False, PADDING)
                    var_start -= len_extend
                else: # both sides are extensible
                    r_ref_extend = self.f_fasta.fetch(reference=ref_name, start=var_stop, end=var_stop+EXTEND_LIMIT+PADDING)
                    r_seq_hap0, r_seq_hap1, r_len_extend = self.extend_ref_seq_padding(seq_hap0, seq_hap1, r_ref_extend, r_ref_extend, True, PADDING)
                    l_ref_extend = self.f_fasta.fetch(reference=ref_name, start=var_start-EXTEND_LIMIT-PADDING, end=var_start)
                    l_seq_hap0, l_seq_hap1, l_len_extend = self.extend_ref_seq_padding(seq_hap0, seq_hap1, l_ref_extend, l_ref_extend, False, PADDING)
                    if l_len_extend == 0: # right anyway
                        seq_hap0, seq_hap1, var_stop = r_seq_hap0, r_seq_hap1, var_stop + r_len_extend
                        len_extend = r_len_extend
                    elif r_len_extend == 0: # left anyway
                        seq_hap0, seq_hap1, var_start = l_seq_hap0, l_seq_hap1, var_start - l_len_extend
                        len_extend = l_len_extend
                    elif r_len_extend < l_len_extend: # right is better
                        seq_hap0, seq_hap1, var_stop = r_seq_hap0, r_seq_hap1, var_stop + r_len_extend
                        len_extend = r_len_extend
                    else: # left is better
                        seq_hap0, seq_hap1, var_start = l_seq_hap0, l_seq_hap1, var_start - l_len_extend
                        len_extend = l_len_extend
                
                if len_extend == False:
                    cohort_vars = list(self.f_vcf.fetch(ref_name, var_start+PADDING, var_stop-PADDING+EXTEND_LIMIT))
                    if (flag_side != 1) and (len(cohort_vars) > 1):
                        var_stop += EXTEND_LIMIT
                    else:
                        logger.warning(f"WARNING! Variant at contig: {ref_name}, pos: {var.start} exceed repetitive length limit: {EXTEND_LIMIT}, give up analyzing it." )
                        self.dict_set_conflict_vars[ref_name].add(var.start)
        

        # check if there are nearby variants between var_start to var_stop
        cohort_vars = list(self.f_vcf.fetch(ref_name, var_start+PADDING-VAR_CHAIN, var_stop-PADDING+VAR_CHAIN))
        if len(cohort_vars) == 1: # single variant within a region
            if self.dict_ref_haps[ref_name].get((var.start)):
                logger.warning(f"WARNING! Duplicate single variant at contig: {ref_name}, pos: {var.start}")
            self.dict_ref_haps[ref_name][(var.start)] = (var_start+PADDING, var_stop-PADDING, seq_hap0, seq_hap1) # left and right anchor point, hap0 and hap1 sequence
            #self.dict_ref_var_pos[ref_name][0].append(var.start)
            #self.dict_ref_var_pos[ref_name][1].append(var.stop)
            return idx_vcf + 1 # While Loop Management
        else: # the case where multiple variants are in the chaining area
            # Expanding to the chaining variants' chaining area
            cohort_start = min(var.start-VAR_CHAIN, min([v.start-VAR_CHAIN for v in cohort_vars]))
            cohort_maxstop = var.stop+VAR_CHAIN
            for v in cohort_vars:
                len_extend = self.check_right_repeat(v)
                cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) + len_extend + VAR_CHAIN for a in v.alleles])) # if no extend, len_extend would be 0
            # Iterate until there are no variants in the chaining area
            while cohort_vars != list(self.f_vcf.fetch(ref_name, cohort_start, cohort_maxstop)):
                cohort_vars = list(self.f_vcf.fetch(ref_name, cohort_start, cohort_maxstop))
                cohort_start = min(cohort_start, min([v.start-VAR_CHAIN for v in cohort_vars]))
                for v in cohort_vars:
                    len_extend = self.check_right_repeat(v)
                    cohort_maxstop = max(cohort_maxstop, max([v.start + len(a) + len_extend + VAR_CHAIN for a in v.alleles])) # if no extend, len_extend would be 0 
            # End of expansion

            # Iterative parameters
            # get the whole cohort sequence, and cut each smaller section for dict_ref_haps
            ref_seq = self.f_fasta.fetch(reference=ref_name, start=cohort_start, end=cohort_maxstop)
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
                    seq_hap0, diff_hap0, len_var= self.switch_var_seq(c_var, seq_hap0, adj_hap0, hap_0)
                    prev_start0 = c_var.start
                    overlap0 = len_var - 1 if (diff_hap0 == 0) else diff_hap0
                    list_start_hap[0].append(c_var.start - adj_hap0)
                    list_len_hap[0].append(len_var)
                else: # overlapping variants are consider conflicts
                    list_start_hap[0].append(-1)    # house keeping
                    list_len_hap[0].append(-1)      # house keeping
                    conflict_flag = True            # conflicts in the cohort
                    self.dict_set_conflict_vars[ref_name].add(prev_start0)
                    self.dict_set_conflict_vars[ref_name].add(c_var.start)
                if c_var.start > prev_start1 + overlap1:
                    adj_hap1 += diff_hap1
                    seq_hap1, diff_hap1, len_var = self.switch_var_seq(c_var, seq_hap1, adj_hap1, hap_1)
                    prev_start1 = c_var.start
                    overlap1 = len_var - 1 if (diff_hap1 == 0) else diff_hap1
                    list_start_hap[1].append(c_var.start - adj_hap1)
                    list_len_hap[1].append(len_var)
                else:
                    list_start_hap[1].append(-1)
                    list_len_hap[1].append(-1)
                    conflict_flag = True
                    self.dict_set_conflict_vars[ref_name].add(prev_start1)
                    self.dict_set_conflict_vars[ref_name].add(c_var.start)
            # complete the seq_hap0 and seq_hap1, do the assignment of dict_ref_seq
            for idx, c_var in enumerate(cohort_vars):
                start0 = list_start_hap[0][idx]
                start1 = list_start_hap[1][idx]
                seq_0 = seq_hap0[start0 - PADDING:start0 + list_len_hap[0][idx] + PADDING]
                seq_1 = seq_hap1[start1 - PADDING:start1 + list_len_hap[1][idx] + PADDING]
                # EXTEND if hap0 and hap1 is indistinguishable
                if (seq_0 != seq_1) and ((seq_0 in seq_1) or (seq_1 in seq_0)): # make sure it is a 0/1 haplotypes and got repetitive issue
                    c_var_start = c_var.start
                    c_var_stop  = c_var.stop
                    idx_hap0_extend = start0 + list_len_hap[0][idx] + PADDING
                    idx_hap1_extend = start1 + list_len_hap[1][idx] + PADDING
                    flag_side = self.left_right_check(seq_0, seq_1) # check which side the repetitive be
                    if flag_side == 1: # left side only
                        seq_0, seq_1, len_extend = self.extend_ref_seq_padding(seq_0, seq_1, seq_hap0[idx_hap0_extend:], seq_hap1[idx_hap1_extend:], False, PADDING)
                        c_var_start -= len_extend
                    else: # flag_side == 0 or 2: # right side
                        seq_0, seq_1, len_extend = self.extend_ref_seq_padding(seq_0, seq_1, seq_hap0[idx_hap0_extend:], seq_hap1[idx_hap1_extend:], True, PADDING)
                        c_var_stop += len_extend
                    if self.dict_ref_haps[ref_name].get((c_var.start)):
                        logger.warning(f"WARNING! Duplicate cohort variant at contig: {ref_name}, pos: {c_var.start}")
                        idx_vcf -= 1
                    self.dict_ref_haps[ref_name][(c_var.start)] = (c_var_start, c_var_stop, seq_0, seq_1)
                    if len_extend == False:
                        logger.warning(f"WARNING! Variant at contig: {ref_name}, pos: {c_var.start} exceed repetitive length limit: {EXTEND_LIMIT}, give up analyzing it. (in cohort)")
                        self.dict_set_conflict_vars[ref_name].add(c_var.start)
                else:
                    if self.dict_ref_haps[ref_name].get((c_var.start)):
                        logger.warning(f"WARNING! Duplicate cohort variant at contig: {ref_name}, pos: {c_var.start}")
                        idx_vcf -= 1
                    self.dict_ref_haps[ref_name][(c_var.start)] = (c_var.start, c_var.stop, seq_0, seq_1)
                    #self.dict_ref_var_pos[ref_name][0].append(c_var.start)
                    #self.dict_ref_var_pos[ref_name][1].append(c_var.stop)
            if not conflict_flag: # only generate the cohort if there are no conflict alleles
                seq_0 = seq_hap0[VAR_CHAIN-PADDING:start0 + list_len_hap[0][idx] + PADDING]
                seq_1 = seq_hap1[VAR_CHAIN-PADDING:start1 + list_len_hap[1][idx] + PADDING]
                max_cohort_stop = cohort_vars[-1].stop

                if (seq_0 != seq_1) and ((seq_0 in seq_1) or (seq_1 in seq_0)): # make sure it is a 0/1 haplotypes and got repetitive issue
                    flag_side = self.left_right_check(seq_0, seq_1) # check which side the repetitive be
                    if flag_side != 1: # right side or both side only
                        seq_0, seq_1, len_extend = self.extend_ref_seq_padding(seq_0, seq_1, seq_hap0[start0 + list_len_hap[0][idx] + PADDING:], \
                                                                                seq_hap1[start1 + list_len_hap[1][idx] + PADDING:], True, PADDING)
                        max_cohort_stop += len_extend
                for idy, c_var in enumerate(cohort_vars): # the start and end position of each var in the cohort
                    lpad_0 = list_start_hap[0][idy] - (VAR_CHAIN-PADDING)
                    lpad_1 = list_start_hap[1][idy] - (VAR_CHAIN-PADDING)
                    rpad_0 = len(seq_0) - lpad_0 - list_len_hap[0][idy]
                    rpad_1 = len(seq_1) - lpad_1 - list_len_hap[1][idy]
                    self.dict_ref_cohorts[ref_name][(c_var.start)] = (var.start, max_cohort_stop, seq_0, seq_1, lpad_0, lpad_1, rpad_0, rpad_1) # c_var should be the last in cohort
            return idx_vcf + len(cohort_vars)


    def handle_repetitive_sequences(self, ref_name: str, var: pysam.VariantRecord, seq_hap0: str, seq_hap1: str, var_start: int, var_stop: int) -> Tuple[str, str, int, int]:
        flag_side = self.left_right_check(seq_hap0, seq_hap1)
        if flag_side in [0, 2]:  # right side or both sides
            ref_extend = self.f_fasta.fetch(reference=ref_name, start=var_stop, end=var_stop+EXTEND_LIMIT+PADDING)
            seq_hap0, seq_hap1, len_extend = self.extend_ref_seq_padding(seq_hap0, seq_hap1, ref_extend, ref_extend, True, PADDING)
            var_stop += len_extend
        if flag_side in [1, 2]:  # left side or both sides
            ref_extend = self.f_fasta.fetch(reference=ref_name, start=var_start-EXTEND_LIMIT-PADDING, end=var_start)
            seq_hap0, seq_hap1, len_extend = self.extend_ref_seq_padding(seq_hap0, seq_hap1, ref_extend, ref_extend, False, PADDING)
            var_start -= len_extend
        return seq_hap0, seq_hap1, var_start, var_stop

    def process_cohort_variants(self, ref_name: str, var: pysam.VariantRecord, cohort_variants: List[pysam.VariantRecord]):
        cohort_start = cohort_variants[0].start - VAR_CHAIN
        cohort_stop = var.stop + VAR_CHAIN
        for c_var in cohort_variants:
            len_extend = self.check_right_repeat(c_var)
            if len_extend:
                cohort_stop = max(cohort_stop, c_var.stop + len_extend)

        
        
        # Process and store cohort information
        # This part needs to be implemented based on your specific requirements for cohort handling

    def check_right_repeat(self, c_var: pysam.VariantRecord) -> Union[bool, int]:
        hap_0, hap_1 = c_var.samples[0]['GT']
        # standard making hap0 and hap1 sequence
        if hap_0 != hap_1: # make sure it is a 0/1 haplotypes
            c_var_start = c_var.start - PADDING
            c_var_stop  = c_var.stop  + PADDING
            ref_seq = self.f_fasta.fetch(reference=c_var.contig, start= c_var_start, end = c_var_stop)
            c_seq_hap0, *_ = self.switch_var_seq(c_var, ref_seq, c_var_start, hap_0)
            c_seq_hap1, *_ = self.switch_var_seq(c_var, ref_seq, c_var_start, hap_1)
            if len(c_seq_hap0) != len(c_seq_hap1):
                if (c_seq_hap0 in c_seq_hap1) or (c_seq_hap1 in c_seq_hap0):
                    # get additional extend_limit (default 70) bp from the reference
                    ref_extend = self.f_fasta.fetch(reference=c_var.contig, start= c_var_stop, end = c_var_stop+EXTEND_LIMIT)
                    c_seq_hap0_extend, c_seq_hap1_extend, len_extend = self.extend_ref_seq(c_seq_hap0, c_seq_hap1, ref_extend, ref_extend, True)
                    return len_extend # right extension successful
        return False 

    def compare_reads_to_variants(self):
        logger.info("Start comparing reads to the variant map...")
        for segment in self.f_sam.fetch(self.chromosome):
            if segment.is_unmapped:
                continue
            self.process_read(segment)

    def process_read(self, segment: pysam.AlignedSegment):
        ref_name     = segment.reference_name
        seq_name     = segment.query_name
        flag_read_n  = segment.is_read2
        pos_start    = segment.reference_start # start position in genome coordiante, need +1 for vcf coordinate
        pos_end      = segment.reference_end
        cigar_tuples = segment.cigartuples
        mapq         = segment.mapping_quality
        read_seq     = segment.query_alignment_sequence # aligned sequence without SoftClip part
        if cigar_tuples == None:
            logger.warning(f"WARNING! {seq_name} is an aligned read without CIGAR!")
            return
        
        if self.real_data == False:
            rg_tag = segment.get_tag("RG")
            if '_' in rg_tag:
                chr_tag, hap_tag = rg_tag.split('_')
            else:
                chr_tag = None
                hap_tag = rg_tag

        # find the related variants without f_vcf.fetch
        direct_var_start, related_var_start = self.find_related_vars(ref_name, pos_start, pos_end)
        #fetching the sequence in the read_seq regarding to the variant
        for base_var_start in related_var_start:
            if base_var_start in self.dict_set_conflict_vars[ref_name]: # neglecting the conflict variant sites
                continue

            match_flag_0 = match_flag_1 = 0
            
            
            # Try cohort alignment first with exact matching
            if self.dict_ref_cohorts[ref_name].get(base_var_start):
                cohort_start, cohort_stop, cohort_seq0, cohort_seq1, lpad_0, lpad_1, rpad_0, rpad_1 = self.dict_ref_cohorts[ref_name][base_var_start]
                
                # Try exact matches for hap0 at both anchor points
                match_flag_0 = self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                               cohort_start, read_seq, cohort_seq0,
                                               cigar_tuples, PADDING, lpad_0+1, rpad_0+1, True)
                if match_flag_0 != 1:
                    match_flag_0 = max(match_flag_0, self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                                   cohort_stop, read_seq, cohort_seq0,
                                                   cigar_tuples, PADDING, lpad_0+1, rpad_0+1, False))

                # Try exact matches for hap1 at both anchor points
                match_flag_1 = self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                               cohort_start, read_seq, cohort_seq1,
                                               cigar_tuples, PADDING, lpad_1+1, rpad_1+1, True)
                if match_flag_1 != 1:
                    match_flag_1 = max(match_flag_1, self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                                   cohort_stop, read_seq, cohort_seq1,
                                                   cigar_tuples, PADDING, lpad_1+1, rpad_1+1, False))

            # Try local alignment if cohort didn't give definitive result
            #if match_flag_0 != 1 and match_flag_1 != 1:
            if match_flag_0 == match_flag_1:
                var_start, var_stop, seq_hap0, seq_hap1, *_ = self.dict_ref_haps[ref_name][base_var_start]
                
                # Try exact matches for both haplotypes at both anchor points
                match_flag_0 = self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                               var_start, read_seq, seq_hap0,
                                               cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
                if match_flag_0 != 1:
                    match_flag_0 = max(match_flag_0, self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                                   var_stop, read_seq, seq_hap0,
                                                   cigar_tuples, PADDING, PADDING+1, PADDING+1, False))

                match_flag_1 = self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                               var_start, read_seq, seq_hap1,
                                               cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
                if match_flag_1 != 1:
                    match_flag_1 = max(match_flag_1, self.match_to_hap(seq_name, ref_name, pos_start, pos_end,
                                                   var_stop, read_seq, seq_hap1,
                                                   cigar_tuples, PADDING, PADDING+1, PADDING+1, False))

                # If both exact matches failed, try edit distance at both anchor points
                if (match_flag_0 != -1 and match_flag_1 != -1) and (match_flag_0 != 1 and match_flag_1 != 1) \
                    and abs(len(seq_hap0) - len(seq_hap1)) <= self.max_len_diff:
                    edit_dist_0 = min(
                        self.get_edit_distance(read_seq, seq_hap0, pos_start, var_start,
                                             cigar_tuples, PADDING, PADDING+1, PADDING+1, True),
                        self.get_edit_distance(read_seq, seq_hap0, pos_start, var_stop,
                                             cigar_tuples, PADDING, PADDING+1, PADDING+1, False)
                    )
                    edit_dist_1 = min(
                        self.get_edit_distance(read_seq, seq_hap1, pos_start, var_start,
                                             cigar_tuples, PADDING, PADDING+1, PADDING+1, True),
                        self.get_edit_distance(read_seq, seq_hap1, pos_start, var_stop,
                                             cigar_tuples, PADDING, PADDING+1, PADDING+1, False)
                    )
                    
                    # Assign based on minimum edit distance
                    if edit_dist_0/len(seq_hap0) < edit_dist_1/len(seq_hap1) and edit_dist_0 < 5: #len(seq_hap0) * 0.2:
                        match_flag_0 = 1
                    elif edit_dist_1/len(seq_hap1) < edit_dist_0/len(seq_hap0) and edit_dist_1 < 5: #len(seq_hap1) * 0.2:
                        match_flag_1 = 1



            # Update statistics based on matching results
            if base_var_start not in direct_var_start:
                if not ((match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_1 == 1 and  match_flag_0 !=1)):
                    continue
            if match_flag_0 == -1 and match_flag_1 == -1:
                continue
            if match_flag_0 == 1 and match_flag_1 == 1:
                self.dict_ref_bias[ref_name][base_var_start]['n_var'][2] += 1
            elif match_flag_0 == 1:
                self.dict_ref_bias[ref_name][base_var_start]['n_var'][0] += 1
                # record the starting position of each read cover the variant
                self.dict_ref_bias[ref_name][base_var_start]['distribute'][0].append(pos_start)
                self.dict_ref_bias[ref_name][base_var_start]['distribute'][2].append(pos_end)
            elif match_flag_1 == 1:
                self.dict_ref_bias[ref_name][base_var_start]['n_var'][1] += 1
                # record the starting position of each read cover the variant
                self.dict_ref_bias[ref_name][base_var_start]['distribute'][1].append(pos_start)
                self.dict_ref_bias[ref_name][base_var_start]['distribute'][3].append(pos_end)
            else:
                self.dict_ref_bias[ref_name][base_var_start]['n_var'][3] += 1
            
            # standard updating of read number and mapping quality
            if self.real_data == True: # no golden information
                self.dict_ref_bias[ref_name][base_var_start]['n_read'][0] += 1
                self.dict_ref_bias[ref_name][base_var_start]['map_q'][0]  += mapq
            else:
                if self.run_id != None and self.run_id != chr_tag: # not the same chromosome
                    self.dict_ref_bias[ref_name][base_var_start]['n_read'][2] += 1
                    self.dict_ref_bias[ref_name][base_var_start]['map_q'][2] += 1
                elif self.dict_ref_var_name[ref_name].get(base_var_start) == None:
                    continue
                elif 'hapA' == hap_tag: # hapA
                    if len(self.dict_ref_var_name[ref_name][base_var_start]) == 6 and \
                       (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][4]: # read is covering the variant but not stretch
                        self.dict_ref_bias[ref_name][base_var_start]['map_q'][2] += mapq         
                        if (match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_0 != 1 and match_flag_1 == 1):
                            #print("!!!!!!!\t\t", var_start, seq_name, flag_read_n, 'hapA', match_flag_0 == 1)
                            pass
                        continue
                    if (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][0]: # check if the read name is in the golden set
                        self.dict_ref_bias[ref_name][base_var_start]['n_read'][0] += 1
                        self.dict_ref_bias[ref_name][base_var_start]['map_q'][0]  += mapq
                    else:
                        self.dict_ref_bias[ref_name][base_var_start]['n_read'][2] += 1
                        self.dict_ref_bias[ref_name][base_var_start]['map_q'][2] += 1
                elif 'hapB' == hap_tag: # hapB
                    if len(self.dict_ref_var_name[ref_name][base_var_start]) == 6 and \
                       (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][5]: # read is covering the variant but not stretch
                        self.dict_ref_bias[ref_name][base_var_start]['map_q'][2] += mapq         
                        if (match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_0 != 1 and match_flag_1 == 1):
                            #print("!!!!!!!\t\t", var_start, seq_name, flag_read_n, 'hapB', match_flag_1 == 1)
                            pass
                        continue
                    if (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][1]: # check if the read name is in the golden set
                        self.dict_ref_bias[ref_name][base_var_start]['n_read'][1] += 1
                        self.dict_ref_bias[ref_name][base_var_start]['map_q'][1]  += mapq
                    else:
                        self.dict_ref_bias[ref_name][base_var_start]['n_read'][2] += 1
                        self.dict_ref_bias[ref_name][base_var_start]['map_q'][2] += 1
                else:
                    logger.warning("WARNING, there is a read without haplotype information!!")

    def find_related_vars(self, ref_name: str, pos_start: int, pos_end: int) -> Tuple[List[int], List[int]]:
        direct_var_start = self.binary_search_variants(ref_name, pos_start, pos_end)
        related_var_start = direct_var_start
        if len(direct_var_start) > 0:
            new_start = pos_start
            new_end   = pos_end
            if self.dict_ref_cohorts[ref_name].get(direct_var_start[0]):
                new_start = self.dict_ref_cohorts[ref_name][direct_var_start[0]][0]  # cohort start
            if self.dict_ref_cohorts[ref_name].get(direct_var_start[-1]):
                new_end   = self.dict_ref_cohorts[ref_name][direct_var_start[-1]][1] # cohort end
            if new_start < pos_start or new_end > pos_end:
                related_var_start = self.binary_search_variants(ref_name, new_start, new_end)
        return direct_var_start, related_var_start
    
    def binary_search_variants(self, ref_name: str, pos_start: int, pos_end: int) -> List[int]:
        list_start, list_end = self.dict_ref_var_pos[ref_name]
        # binary search the range that pos_start > list_end[idx_2] and pos_end < list_start[idx_1], return list_start[idx_1:idx_2]
        idx_start = self.bisect_left(list_end, pos_start)
        idx_end   = self.bisect_right(list_start, pos_end)
        return list_start[idx_start:idx_end]
    
    def bisect_left(self, list_end: List[int], pos_start: int) -> int:
        left, right = 0, len(list_end)
        while left < right:
            mid = (left + right) // 2
            if list_end[mid] < pos_start:
                left = mid + 1
            else:
                right = mid
        return left

    def bisect_right(self, list_start: List[int], pos_end: int) -> int:
        left, right = 0, len(list_start)
        while left < right:
            mid = (left + right) // 2
            if list_start[mid] > pos_end:
                right = mid
            else:
                left = mid + 1
        return left


    def update_bias_data(self, ref_name: str, var_start: int, match_flag_0: int, match_flag_1: int, 
                         pos_start: int, pos_end: int, mapq: int, segment: pysam.AlignedSegment):
        bias_data = self.dict_ref_bias[ref_name][var_start]
        
        if match_flag_0 == 1 and match_flag_1 == 1:
            bias_data['n_var'][2] += 1
        elif match_flag_0 == 1:
            bias_data['n_var'][0] += 1
            bias_data['distribute'][0].append(pos_start)
            bias_data['distribute'][2].append(pos_end)
        elif match_flag_1 == 1:
            bias_data['n_var'][1] += 1
            bias_data['distribute'][1].append(pos_start)
            bias_data['distribute'][3].append(pos_end)
        else:
            bias_data['n_var'][3] += 1

        if self.real_data:
            bias_data['n_read'][0] += 1
            bias_data['map_q'][0] += mapq
        else:
            self.update_simulated_data(ref_name, var_start, segment, bias_data, mapq)

    def update_simulated_data(self, ref_name: str, var_start: int, segment: pysam.AlignedSegment, 
                              bias_data: Dict, mapq: int):
        rg_tag = segment.get_tag("RG")
        chr_tag, hap_tag = rg_tag.split('_') if '_' in rg_tag else (None, rg_tag)
        
        if self.run_id is not None and self.run_id != chr_tag:
            bias_data['n_read'][2] += 1
            bias_data['map_q'][2] += mapq
        elif self.dict_ref_var_name[ref_name].get(var_start) is None:
            return
        elif hap_tag == 'hapA':
            self.update_hap_data(ref_name, var_start, segment, bias_data, mapq, 0)
        elif hap_tag == 'hapB':
            self.update_hap_data(ref_name, var_start, segment, bias_data, mapq, 1)
        else:
            logger.warning("There is a read without haplotype information!!")

    def update_hap_data(self, ref_name: str, var_start: int, segment: pysam.AlignedSegment, 
                        bias_data: Dict, mapq: int, hap_index: int):
        var_name_data = self.dict_ref_var_name[ref_name][var_start]
        if len(var_name_data) == 6 and (segment.query_name, segment.is_read2) in var_name_data[4 + hap_index]:
            bias_data['map_q'][2] += mapq
        elif (segment.query_name, segment.is_read2) in var_name_data[hap_index]:
            bias_data['n_read'][hap_index] += 1
            bias_data['map_q'][hap_index] += mapq
        else:
            bias_data['n_read'][2] += 1
            bias_data['map_q'][2] += mapq

    def generate_report(self):
        logger.info("Start generating report...")
        return self.output_report()

    def output_report(self):
        all_output = []
        gap_output = []
        snp_output = []

        self.f_vcf.reset()
        for var in self.f_vcf.fetch(self.chromosome):
            all_line, gap_line, snp_line = self.format_variant_data(var)
            if all_line == None:
                continue
            all_output.append(all_line)
            if gap_line:
                gap_output.append(gap_line)
            if snp_line:
                snp_output.append(snp_line)

        return {
            'all': ''.join(all_output),
            'gap': ''.join(gap_output),
            'SNP': ''.join(snp_output)
        }

    def format_variant_data(self, var):
        ref_name = var.contig
        hap = var.samples[0]['GT']
        if (hap[0] != 0 and hap[1] != 0) or (hap[0] == 0 and hap[1] == 0):
            return None, None, None
        if hap[0] == 0:
            idx_ref, idx_alt = 0, 1
        else:
            idx_ref, idx_alt = 1, 0
        
        if var.start in self.dict_set_conflict_vars[ref_name]:
            return None, None, None

        bias_data = self.dict_ref_bias[ref_name][var.start]
        n_read = bias_data['n_read']
        n_var = bias_data['n_var']
        map_q = bias_data['map_q']
        p_value = self.chi_square_test(var.start, bias_data['distribute'][idx_alt])
        p_value = min(p_value, self.chi_square_test(var.start, bias_data['distribute'][idx_ref]))

        output_string = f"{ref_name}\t{var.start+1}\t{sum(n_read)}\t{self.get_division(sum(map_q[:2]), sum(n_read[:2]))}\t{p_value:.4f}\t"
        output_string += f"{self.get_division(n_var[idx_ref]+n_var[2]*0.5, sum(n_var[:3]))}\t{n_var[idx_ref]}\t{n_var[idx_alt]}\t{n_var[2]}\t{n_var[3]}"

        if not self.real_data:
            output_string += f"\t{self.get_division(n_read[idx_ref], sum(n_read[:2]))}\t{n_read[idx_ref]}\t{n_read[idx_alt]}\t{n_read[2]}"
            read_info = self.dict_ref_var_name[ref_name][var.start]
            output_string += f"\t{self.get_division(read_info[idx_ref+2], sum(read_info[2:4]))}\t{read_info[idx_ref+2]}\t{read_info[idx_alt+2]}"

        all_line = output_string + "\t"
        gap_line = None
        snp_line = None

        if len(var.ref) == len(var.alts[hap[idx_alt] - 1]):
            all_line += "\n"
            snp_line = output_string + "\n"
        else:
            all_line += ".\n"
            gap_line = output_string + "\n"
        return all_line, gap_line, snp_line

    def write_variant_data(self, var: pysam.VariantRecord, f_all, f_gap, f_SNP):
        ref_name = var.contig
        hap = var.samples[0]['GT']
        if (hap[0] != 0 and hap[1] != 0) or (hap[0] == 0 and hap[1] == 0):
            return
        if hap[0] == 0:
            idx_ref, idx_alt = 0, 1
        else:
            idx_ref, idx_alt = 1, 0
        
        if var.start in self.dict_set_conflict_vars[ref_name]:
            return

        bias_data = self.dict_ref_bias[ref_name][var.start]
        n_read = bias_data['n_read']
        n_var = bias_data['n_var']
        map_q = bias_data['map_q']
        p_value = self.chi_square_test(var.start, bias_data['distribute'][idx_alt])
        p_value = min(p_value, self.chi_square_test(var.start, bias_data['distribute'][idx_ref]))

        output_string = f"{ref_name}\t{var.start+1}\t{sum(n_read)}\t{self.get_division(sum(map_q[:2]), sum(n_read[:2]))}\t{p_value:.4f}\t"
        output_string += f"{self.get_division(n_var[idx_ref]+n_var[2]*0.5, sum(n_var[:3]))}\t{n_var[idx_ref]}\t{n_var[idx_alt]}\t{n_var[2]}\t{n_var[3]}"

        if not self.real_data:
            output_string += f"\t{self.get_division(n_read[idx_ref], sum(n_read[:2]))}\t{n_read[idx_ref]}\t{n_read[idx_alt]}\t{n_read[2]}"
            read_info = self.dict_ref_var_name[ref_name][var.start]
            output_string += f"\t{self.get_division(read_info[idx_ref+2], sum(read_info[2:4]))}\t{read_info[idx_ref+2]}\t{read_info[idx_alt+2]}"

        if len(var.ref) == len(var.alts[hap[idx_alt] - 1]):
            f_all.write(output_string + "\t\n")
            f_SNP.write(output_string + "\n")
        else:
            f_all.write(output_string + "\t.\n")
            f_gap.write(output_string + "\n")

    def chi_square_test(self, var_start: int, list_pos_start: List[int]) -> float:
        if len(list_pos_start) < 2:
            return 0
        bucket_num = 5
        bucket_len = int(100 / bucket_num)
        list_count = np.zeros(bucket_num)
        input_idx = np.minimum((var_start - np.array(list_pos_start)) // bucket_len, bucket_num - 1)
        try:
            np.add.at(list_count, input_idx, 1)
        except IndexError:
            print(var_start, list_pos_start)
        _, p_value = chisquare(list_count)
        return 0 if np.isnan(p_value) else p_value

    def get_division(self, num_1: float, num_2: float) -> str:
        if num_2 == 0:
            return 'nan'
        else:
            return format(num_1 / num_2, '.4f')

    def switch_var_seq(self, var: pysam.VariantRecord, ref: str, start: int, genotype: int) -> Tuple[str, int, int]:
        if genotype == 0:
            return ref, 0, len(var.ref)
        else:
            alt = var.alts[genotype - 1]
            return ref[:var.start-start] + alt + ref[var.stop-start:], len(var.ref) - len(alt), len(alt)

    def left_right_check(self, seq_hap0: str, seq_hap1: str) -> int:
        assert seq_hap0 != seq_hap1
        assert (seq_hap0 in seq_hap1) or (seq_hap1 in seq_hap0)
        len_0, len_1 = len(seq_hap0), len(seq_hap1)
        if len_0 > len_1:
            if seq_hap0[:len_1] == seq_hap1:
                return 0  # right side repetitive
            elif seq_hap0[-len_1:] == seq_hap1:
                return 1  # left side repetitive
        else:
            if seq_hap1[:len_0] == seq_hap0:
                return 0  # right side repetitive
            elif seq_hap1[-len_0:] == seq_hap0:
                return 1  # left side repetitive
        return 2  # in the middle

    def extend_ref_seq_padding(self, seq_hap0: str, seq_hap1: str, ref_extend_0: str, ref_extend_1: str, 
                               flag_right: bool = True, padding: int = PADDING) -> Tuple[str, str, int]:
        if flag_right:
            seq_hap0_extend, seq_hap1_extend, len_extend = self.extend_ref_seq(seq_hap0, seq_hap1, ref_extend_0[:-padding], ref_extend_1[:-padding], flag_right)
            if len_extend:
                return seq_hap0_extend + ref_extend_0[len_extend:len_extend+padding], seq_hap1_extend + ref_extend_1[len_extend:len_extend+padding], len_extend+padding
            else:
                return seq_hap0, seq_hap1, False
        else:
            seq_hap0_extend, seq_hap1_extend, len_extend = self.extend_ref_seq(seq_hap0, seq_hap1, ref_extend_0[padding:], ref_extend_1[padding:], flag_right)
            if len_extend:
                return ref_extend_0[-len_extend-padding:-len_extend] + seq_hap0_extend, ref_extend_1[-len_extend-padding:-len_extend] + seq_hap1_extend, len_extend+padding
            else:
                return seq_hap0, seq_hap1, False

    def extend_ref_seq(self, seq_hap0: str, seq_hap1: str, ref_extend_0: str, ref_extend_1: str, flag_right: bool = True) -> Tuple[str, str, int]:
        seq_hap0_extend, seq_hap1_extend = seq_hap0, seq_hap1
        assert (seq_hap0_extend in seq_hap1_extend) or (seq_hap1_extend in seq_hap0_extend)
        len_iterate = min(len(ref_extend_0), len(ref_extend_1))
        if flag_right:
            for idx in range(len_iterate):
                seq_hap0_extend += ref_extend_0[idx]
                seq_hap1_extend += ref_extend_1[idx]
                if not ((seq_hap0_extend in seq_hap1_extend) or (seq_hap1_extend in seq_hap0_extend)):
                    return seq_hap0_extend, seq_hap1_extend, idx+1
        else:
            for idx in range(len_iterate):
                seq_hap0_extend = ref_extend_0[-idx-1] + seq_hap0_extend
                seq_hap1_extend = ref_extend_1[-idx-1] + seq_hap1_extend
                if not ((seq_hap0_extend in seq_hap1_extend) or (seq_hap1_extend in seq_hap0_extend)):
                    return seq_hap0_extend, seq_hap1_extend, idx+1
        return seq_hap0_extend, seq_hap1_extend, False

    def match_to_hap(self, seq_name: str, contig: str, read_start: int, read_end: int, var_start: int,
                     seq_read: str, seq_hap: str, cigar_tuples: List[Tuple[int, int]],
                     padding: int, l_min_req: int, r_min_req: int, start_flag: bool = True) -> int:
        if read_start > var_start or read_end < var_start:
            return -1
        r_start = self.locate_by_cigar(read_start, var_start, cigar_tuples)
        
        if start_flag:
            l_bound = r_start - padding
            r_bound = l_bound + len(seq_hap)
        else:
            r_bound = r_start + padding
            l_bound = r_bound - len(seq_hap)

        min_match = 0
        if l_bound < 0:
            seq_hap = seq_hap[-l_bound:]
            l_bound = 0
            min_match = r_min_req
        if r_bound > len(seq_read):
            seq_hap = seq_hap[:len(seq_read)-r_bound]
            r_bound = len(seq_read)
            if min_match != 0:
                #logger.warning(f"The read is totally contained in the variant!! {seq_name} at {contig}, {var_start}")
                pass
            min_match = l_min_req
        if r_bound - l_bound < min_match:
            return -1
        if seq_read[l_bound:r_bound].upper() == seq_hap.upper():
            return 1
        else:
            return 0

    def locate_by_cigar(self, read_start: int, target_pos: int, cigar_tuples: List[Tuple[int, int]]) -> int:
        ref_cursor = read_start
        read_cursor = 0
        for code, runs in cigar_tuples:
            if code in (0, 7, 8):  # M or = or X
                ref_cursor += runs
                if ref_cursor > target_pos:
                    return read_cursor + (runs - ref_cursor + target_pos)
                read_cursor += runs
            elif code == 1:  # I
                if ref_cursor > target_pos:
                    return read_cursor
                read_cursor += runs
            elif code == 2:  # D
                ref_cursor += runs
                if ref_cursor > target_pos:
                    return read_cursor
            elif code in (4, 5):  # S or H, pysam already parsed
                pass
            else:
                logger.error(f"Unexpected CIGAR code {code}")
        return read_cursor

    def get_edit_distance(self, seq_read: str, seq_hap: str, pos_start: int, var_start: int,
                         cigar_tuples: List[Tuple[int, int]], padding: int,
                         l_min_req: int, r_min_req: int, start_flag: bool = True) -> int:
        r_start = self.locate_by_cigar(pos_start, var_start, cigar_tuples)
        
        if start_flag:
            l_bound = r_start - padding
            r_bound = l_bound + len(seq_hap)
        else:
            r_bound = r_start + padding
            l_bound = r_bound - len(seq_hap)

        # Boundary checks
        if l_bound < 0:
            seq_hap_adj = seq_hap[-l_bound:]
            l_bound = 0
            if len(seq_hap_adj) < r_min_req:
                return float('inf')
        if r_bound > len(seq_read):
            seq_hap_adj = seq_hap[:len(seq_read)-r_bound]
            r_bound = len(seq_read)
            if len(seq_hap_adj) < l_min_req:
                return float('inf')
        
        read_segment = seq_read[l_bound:r_bound].upper()
        hap_segment = seq_hap.upper()
        
        result = edlib.align(read_segment, hap_segment, task='distance')
        return result['editDistance']


def analyze_read(list_info):
    path_vcf, path_bam, path_ref, golden_pickle, run_id, real_data, out, chromosome, max_len_diff = list_info
    analyzer = VariantAnalyzer(path_vcf, path_bam, path_ref, 
                                   golden_pickle, run_id, 
                                   real_data, out, chromosome, max_len_diff)
    return analyzer.analyze()

def write_report(path_out, result, real_data: bool):
    f_all = open(path_out + ".all", "w")
    f_gap = open(path_out + ".gap", "w")
    f_snp = open(path_out + ".snp", "w")

    headers = get_report_headers(real_data)
    f_all.write(headers['all'])
    f_gap.write(headers['gap'])
    f_snp.write(headers['SNP'])

    for element in result:
        f_all.write(element['all'])
        f_gap.write(element['gap'])
        f_snp.write(element['SNP'])
    f_all.close()
    f_gap.close()
    f_snp.close()

def get_report_headers(real_data: bool) -> Dict[str, str]:
    if real_data:
        header = "CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tBOTH\tOTHER\tGAP\n"
        return {'all': header, 'gap': header, 'SNP': header}
    else:
        header = "CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tBOTH\tOTHER\tMAP_BALANCE\tMAP_REF\tMAP_ALT\tMIS_MAP\tSIM_BALANCE\tSIM_REF\tSIM_ALT\tGAP\n"
        return {'all': header, 'gap': header, 'SNP': header}




def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze reference bias in genomic sequencing data.")
    parser.add_argument('-v', '--vcf', help='VCF file', required=True)
    parser.add_argument('-s', '--sam', help='SAM file', required=True)
    parser.add_argument('-f', '--fasta', help='Reference FASTA file', required=True)
    parser.add_argument('-r', '--real_data', help='Turn off hap_information warning for real data', action='store_true')
    parser.add_argument('-p', '--golden_pickle', help='Pickle file containing the golden information for report reference')
    parser.add_argument('-i', '--run_id', help='Tag for run_id, can be used to indicate chromosome number')
    parser.add_argument('-t', '--thread', help='Number of threads', type=int, default=8)
    parser.add_argument('-o', '--out', help='Output file', required=True)
    parser.add_argument('-m', '--max_len_diff', type=int, default=20,
                       help='Maximum length difference between haplotypes for edit distance calculation [20]')
    return parser.parse_args()


def main():
    args = parse_arguments()
    bam_file = args.sam
    vcf_file = args.vcf
    ref_file = args.fasta
    pickle_file = args.golden_pickle
    run_id = args.run_id
    real_data = args.real_data
    out = args.out
    thread = args.thread
    max_len_diff = args.max_len_diff
    list_chromosome = []
    for chromosome in pysam.VariantFile(vcf_file).header.contigs.keys():
        if 'random' in chromosome or 'chrUn' in chromosome or 'alt' in chromosome or 'chrEBV' in chromosome or 'chrM' in chromosome:
            continue
        list_chromosome.append(chromosome)
    
    list_info = [(vcf_file, bam_file, ref_file, pickle_file, run_id, real_data, out, chromosome, max_len_diff) for chromosome in list_chromosome]
    with multiprocessing.Pool(thread) as pool:
        result = pool.map(analyze_read, list_info)
    write_report(out, result, real_data)



if __name__ == "__main__":
    main()
