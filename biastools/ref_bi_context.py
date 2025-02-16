import argparse
import pickle
import pysam
import numpy as np
from scipy.stats import chisquare
from typing import List, Tuple, Dict, Union, Set
from collections import defaultdict
import logging
import concurrent.futures
import multiprocessing

from contextlib import contextmanager
import threading
import time
import traceback
import sys

import pstats
from line_profiler import profile



def dump_thread_stacks():
    """
    Dumps the stack traces of all threads.
    """
    for thread_id, stack in sys._current_frames().items():
        print(f"\nThread ID: {thread_id}")
        for filename, lineno, name, line in traceback.extract_stack(stack):
            print(f"File: {filename}, Line: {lineno}, Function: {name}")
            if line:
                print(f"  {line.strip()}")

@contextmanager
def track_lock(lock, name=''):
    start = time.time()
    logger.info(f"{name}: Attempting to acquire...")
    lock.acquire()
    try:
        logger.info(f"{name}: Acquired after {time.time() - start:.4f} seconds")
        #dump_thread_stacks()
        yield
    finally:
        lock.release()
        logger.info(f"{name}: Released")

lock = threading.Lock()

# Constants
PADDING = 5
VAR_CHAIN = 25
EXTEND_LIMIT = 70

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_chromosome(chromosome: str, sam_file: str, dict_haps: Dict, dict_cohorts: Dict, 
                       dict_set_conflict_vars: Set, dict_var_pos: Dict, real_data: bool, run_id: str = None) -> Dict:
    logger.info(f"Processing chromosome: {chromosome}")
    f_sam = pysam.AlignmentFile(sam_file + f".{chromosome}.bam", "rb")
    dict_bias = defaultdict(lambda: {'n_read': [0, 0, 0], 'n_var': [0, 0, 0, 0], 'map_q': [0, 0, 0], 'distribute': [[], [], [], []]})
    
    for segment in f_sam.fetch(chromosome):
        if segment.is_unmapped:
            continue
        process_read(segment, chromosome, dict_haps, dict_cohorts, dict_set_conflict_vars, dict_var_pos, dict_bias, real_data, run_id)
    
    return dict_bias

def process_read(self, segment: pysam.AlignedSegment):
    """
        ref_name = segment.reference_name
        pos_start = segment.reference_start
        pos_end = segment.reference_end
        cigar_tuples = segment.cigartuples
        mapq = segment.mapping_quality
        read_seq = segment.query_alignment_sequence

        related_vars = list(self.f_vcf.fetch(ref_name, pos_start, pos_end))
        
        for var in related_vars:
            if var.start in self.dict_set_conflict_vars[ref_name]:
                continue

            match_flag_0 = self.match_to_hap(segment.query_name, pos_start, pos_end, var.start, 
                                             read_seq, self.dict_ref_haps[ref_name][var.start][2], 
                                             cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
            match_flag_1 = self.match_to_hap(segment.query_name, pos_start, pos_end, var.start, 
                                             read_seq, self.dict_ref_haps[ref_name][var.start][3], 
                                             cigar_tuples, PADDING, PADDING+1, PADDING+1, True)

            self.update_bias_data(ref_name, var.start, match_flag_0, match_flag_1, pos_start, pos_end, mapq, segment)"""
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
        if base_var_start in self.dict_set_conflict_vars: # neglecting the conflict variant sites
            continue

        # 1 for match, 0 for unmatch, -1 for not cover
        match_flag_0 = 0
        match_flag_1 = 0
            
        # 1. Cohort alignment
        if self.dict_cohorts.get(base_var_start): # Anchor Left
            cohort_start, cohort_stop, cohort_seq0, cohort_seq1, lpad_0, lpad_1, rpad_0, rpad_1 = self.dict_cohorts[base_var_start] 
            match_flag_0 = self.match_to_hap(seq_name, pos_start, pos_end, cohort_start, read_seq, cohort_seq0, cigar_tuples, PADDING, lpad_0+1, rpad_0+1, True)
            if match_flag_0 != 1: # Left doesn't work, anchor Right
                match_flag_0 = max(match_flag_0, self.match_to_hap(seq_name, pos_start, pos_end, cohort_stop, read_seq, cohort_seq0, cigar_tuples, PADDING, lpad_0+1, rpad_0+1, False))
            match_flag_1 = self.match_to_hap(seq_name, pos_start, pos_end, cohort_start, read_seq, cohort_seq1, cigar_tuples, PADDING, lpad_1+1, rpad_1+1, True)
            if match_flag_1 != 1: # Left doesn't work, anchor Right
                match_flag_1 = max(match_flag_1, self.match_to_hap(seq_name, pos_start, pos_end, cohort_stop, read_seq, cohort_seq1, cigar_tuples, PADDING, lpad_1+1, rpad_1+1, False))
            
        # 2. Local alignment
        if match_flag_0 == match_flag_1: # both or others
            var_start, var_stop, seq_hap0, seq_hap1 = self.dict_haps[base_var_start]
            match_flag_0 = self.match_to_hap(seq_name, pos_start, pos_end, var_start, read_seq, seq_hap0, cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
            if match_flag_0 != 1:
                match_flag_0 = max(match_flag_0, self.match_to_hap(seq_name, pos_start, pos_end, var_stop, read_seq, seq_hap0, cigar_tuples, PADDING, PADDING+1, PADDING+1, False))
            match_flag_1 = self.match_to_hap(seq_name, pos_start, pos_end, var_start, read_seq, seq_hap1, cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
            if match_flag_1 != 1:
                match_flag_1 = max(match_flag_1, self.match_to_hap(seq_name, pos_start, pos_end, var_stop, read_seq, seq_hap1, cigar_tuples, PADDING, PADDING+1, PADDING+1, False))
            
        # 3. Assign Values
        if base_var_start not in direct_var_start:
            if not ((match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_1 == 1 and  match_flag_0 !=1)):
                continue
        if match_flag_0 == -1 and match_flag_1 == -1:
            continue
        if match_flag_0 == 1 and match_flag_1 == 1:
            self.dict_bias[base_var_start]['n_var'][2] += 1
        elif match_flag_0 == 1:
            self.dict_bias[base_var_start]['n_var'][0] += 1
            # record the starting position of each read cover the variant
            self.dict_bias[base_var_start]['distribute'][0].append(pos_start)
            self.dict_bias[base_var_start]['distribute'][2].append(pos_end)
        elif match_flag_1 == 1:
            self.dict_bias[base_var_start]['n_var'][1] += 1
            # record the starting position of each read cover the variant
            self.dict_bias[base_var_start]['distribute'][1].append(pos_start)
            self.dict_bias[base_var_start]['distribute'][3].append(pos_end)
        else:
            self.dict_bias[base_var_start]['n_var'][3] += 1
            
        # standard updating of read number and mapping quality
        if self.real_data == True: # no golden information
            self.dict_bias[base_var_start]['n_read'][0] += 1
            self.dict_bias[base_var_start]['map_q'][0]  += mapq
        else:
            if self.run_id != None and self.run_id != chr_tag: # not the same chromosome
                self.dict_bias[base_var_start]['n_read'][2] += 1
                self.dict_bias[base_var_start]['map_q'][2] += 1
            elif self.dict_ref_var_name[ref_name].get(base_var_start) == None:
                continue
            elif 'hapA' == hap_tag: # hapA
                if len(self.dict_ref_var_name[ref_name][base_var_start]) == 6 and \
                   (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][4]: # read is covering the variant but not stretch
                    self.dict_bias[base_var_start]['map_q'][2] += mapq         
                    if (match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_0 != 1 and match_flag_1 == 1):
                        #print("!!!!!!!\t\t", var_start, seq_name, flag_read_n, 'hapA', match_flag_0 == 1)
                        pass
                    continue
                if (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][0]: # check if the read name is in the golden set
                    self.dict_bias[base_var_start]['n_read'][0] += 1
                    self.dict_bias[base_var_start]['map_q'][0]  += mapq
                else:
                    self.dict_bias[base_var_start]['n_read'][2] += 1
                    self.dict_bias[base_var_start]['map_q'][2] += 1
            elif 'hapB' == hap_tag: # hapB
                if len(self.dict_ref_var_name[ref_name][base_var_start]) == 6 and \
                   (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][base_var_start][5]: # read is covering the variant but not stretch
                    self.dict_bias[base_var_start]['map_q'][2] += mapq         
                    if (match_flag_0 == 1 and match_flag_1 != 1) or (match_flag_0 != 1 and match_flag_1 == 1):
                        pass
                    continue
                if (seq_name, flag_read_n) in self.dict_ref_var_name[ref_name][var.start][1]: # check if the read name is in the golden set
                    self.dict_bias[base_var_start]['n_read'][1] += 1
                    self.dict_bias[base_var_start]['map_q'][1]  += mapq
                else:
                    self.dict_bias[base_var_start]['n_read'][2] += 1
                    self.dict_bias[base_var_start]['map_q'][2] += 1
            else:
                logger.warning("WARNING, there is a read without haplotype information!!")

def find_related_vars(ref_name: str, pos_start: int, pos_end: int, dict_var_pos: Dict) -> Tuple[List[int], List[int]]:
    direct_var_start = binary_search_variants(dict_var_pos[ref_name], pos_start, pos_end)
    related_var_start = direct_var_start
    if len(direct_var_start) > 0:
        new_start = pos_start
        new_end = pos_end
        if dict_var_pos[ref_name][0][direct_var_start[0]] < pos_start:
            new_start = dict_var_pos[ref_name][0][direct_var_start[0]]
        if dict_var_pos[ref_name][1][direct_var_start[-1]] > pos_end:
            new_end = dict_var_pos[ref_name][1][direct_var_start[-1]]
        if new_start < pos_start or new_end > pos_end:
            related_var_start = binary_search_variants(dict_var_pos[ref_name], new_start, new_end)
    return direct_var_start, related_var_start

def binary_search_variants(var_pos: List[List[int]], pos_start: int, pos_end: int) -> List[int]:
    list_start, list_end = var_pos
    idx_start = bisect_left(list_end, pos_start)
    idx_end = bisect_right(list_start, pos_end)
    return list_start[idx_start:idx_end]

def bisect_left(list_end: List[int], pos_start: int) -> int:
    left, right = 0, len(list_end)
    while left < right:
        mid = (left + right) // 2
        if list_end[mid] < pos_start:
            left = mid + 1
        else:
            right = mid
    return left

def bisect_right(list_start: List[int], pos_end: int) -> int:
    left, right = 0, len(list_start)
    while left < right:
        mid = (left + right) // 2
        if list_start[mid] > pos_end:
            right = mid
        else:
            left = mid + 1
    return left

def process_variant(base_var_start: int, segment: pysam.AlignedSegment, dict_haps: Dict, dict_cohorts: Dict, 
                    pos_start: int, pos_end: int, read_seq: str, cigar_tuples: List[Tuple[int, int]]) -> Tuple[int, int]:
    match_flag_0 = 0
    match_flag_1 = 0

    if base_var_start in dict_cohorts:
        cohort_start, cohort_stop, cohort_seq0, cohort_seq1, lpad_0, lpad_1, rpad_0, rpad_1 = dict_cohorts[base_var_start]
        match_flag_0 = match_to_hap(segment.query_name, pos_start, pos_end, cohort_start, read_seq, cohort_seq0, cigar_tuples, PADDING, lpad_0+1, rpad_0+1, True)
        if match_flag_0 != 1:
            match_flag_0 = max(match_flag_0, match_to_hap(segment.query_name, pos_start, pos_end, cohort_stop, read_seq, cohort_seq0, cigar_tuples, PADDING, lpad_0+1, rpad_0+1, False))
        match_flag_1 = match_to_hap(segment.query_name, pos_start, pos_end, cohort_start, read_seq, cohort_seq1, cigar_tuples, PADDING, lpad_1+1, rpad_1+1, True)
        if match_flag_1 != 1:
            match_flag_1 = max(match_flag_1, match_to_hap(segment.query_name, pos_start, pos_end, cohort_stop, read_seq, cohort_seq1, cigar_tuples, PADDING, lpad_1+1, rpad_1+1, False))

    if match_flag_0 == match_flag_1:
        var_start, var_stop, seq_hap0, seq_hap1 = dict_haps[base_var_start]
        match_flag_0 = match_to_hap(segment.query_name, pos_start, pos_end, var_start, read_seq, seq_hap0, cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
        if match_flag_0 != 1:
            match_flag_0 = max(match_flag_0, match_to_hap(segment.query_name, pos_start, pos_end, var_stop, read_seq, seq_hap0, cigar_tuples, PADDING, PADDING+1, PADDING+1, False))
        match_flag_1 = match_to_hap(segment.query_name, pos_start, pos_end, var_start, read_seq, seq_hap1, cigar_tuples, PADDING, PADDING+1, PADDING+1, True)
        if match_flag_1 != 1:
            match_flag_1 = max(match_flag_1, match_to_hap(segment.query_name, pos_start, pos_end, var_stop, read_seq, seq_hap1, cigar_tuples, PADDING, PADDING+1, PADDING+1, False))

    return match_flag_0, match_flag_1

def match_to_hap(seq_name: str, read_start: int, read_end: int, var_start: int,
                 seq_read: str, seq_hap: str, cigar_tuples: List[Tuple[int, int]],
                 padding: int, l_min_req: int, r_min_req: int, start_flag: bool = True) -> int:
    if read_start > var_start or read_end < var_start:
        return -1
    r_start = locate_by_cigar(read_start, var_start, cigar_tuples)
    
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
            logger.warning(f"The read is totally contained in the variant!! {seq_name} at {var_start}")
        min_match = l_min_req
    if r_bound - l_bound < min_match:
        return -1
    if seq_read[l_bound:r_bound].upper() == seq_hap.upper():
        return 1
    else:
        return 0

def locate_by_cigar(read_start: int, target_pos: int, cigar_tuples: List[Tuple[int, int]]) -> int:
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

def update_bias_data(dict_bias: Dict, base_var_start: int, match_flag_0: int, match_flag_1: int, 
                     pos_start: int, pos_end: int, mapq: int, real_data: bool, run_id: str, 
                     chr_tag: str, hap_tag: str, seq_name: str, flag_read_n: bool):
    bias_data = dict_bias[base_var_start]
    
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

    if real_data:
        bias_data['n_read'][0] += 1
        bias_data['map_q'][0] += mapq
    else:
        update_simulated_data(bias_data, run_id, chr_tag, hap_tag, seq_name, flag_read_n, mapq)

def update_simulated_data(bias_data: Dict, run_id: str, chr_tag: str, hap_tag: str, 
                          seq_name: str, flag_read_n: bool, mapq: int):
    if run_id is not None and run_id != chr_tag:
        bias_data['n_read'][2] += 1
        bias_data['map_q'][2] += mapq
    elif hap_tag == 'hapA':
        update_hap_data(bias_data, seq_name, flag_read_n, mapq, 0)
    elif hap_tag == 'hapB':
        update_hap_data(bias_data, seq_name, flag_read_n, mapq, 1)
    else:
        logger.warning("There is a read without haplotype information!!")

def update_hap_data(bias_data: Dict, seq_name: str, flag_read_n: bool, mapq: int, hap_index: int):
    bias_data['n_read'][hap_index] += 1
    bias_data['map_q'][hap_index] += mapq






class VariantAnalyzer:
    def __init__(self, vcf_file: str, sam_file: str, fasta_file: str, 
                 golden_pickle: str = None, run_id: str = None, 
                 real_data: bool = False, output_file: str = "output.rpt"):
        self.f_vcf = pysam.VariantFile(vcf_file)
        self.f_sam = pysam.AlignmentFile(sam_file)
        self.f_name = sam_file
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
        
        if not self.real_data and self.golden_pickle:
            with open(self.golden_pickle, "rb") as f:
                self.dict_ref_var_name = pickle.load(f)

    """
    def analyze(self):
        self.build_variant_maps()
        self.lock = threading.Lock()

        # TODO try multiprocessing instead of threading
        list_list_reads = []
        for chrom in self.dict_ref_haps.keys():
            tmp_f_sam = pysam.AlignmentFile(self.f_name + f".{chrom}" + ".bam", "r")
            list_list_reads.append(tmp_f_sam.fetch(chrom))

        with concurrent.futures.ThreadPoolExecutor() as executor:
            #futures = {executor.submit(self.compare_reads_to_variants, chrom) for chrom in self.f_vcf.header.contigs}
            futures = {executor.submit(self.process_chromosome, chrom, list_list_reads[idx], self.dict_ref_haps[chrom].copy(), self.dict_ref_cohorts[chrom].copy(), \
                                       self.dict_set_conflict_vars[chrom].copy(), self.dict_ref_var_pos[chrom].copy(), self.real_data, self.run_id) \
                       for idx, chrom in enumerate(self.dict_ref_haps.keys())}
            for future in concurrent.futures.as_completed(futures):
                #self.compare_reads_to_variants()
                print(future, "executed!")
                print(future.result())
        self.generate_report()

    def process_chromosome(self, chromosome, list_reads, dict_ref_haps, dict_ref_cohorts, dict_set_conflict_vars, dict_ref_var_pos, real_data, run_id):
        #with track_lock(self.lock, name=f"Chromosome {chromosome} Lock"):
        processor = ChromosomeProcessor(chromosome, list_reads, dict_ref_haps, dict_ref_cohorts, dict_set_conflict_vars, dict_ref_var_pos, real_data, run_id)
        processor.compare_reads_to_variants(chromosome)
        self.dict_ref_bias[chromosome] = processor.dict_bias.copy()
        #logger.info(f"Finished processing for chromosome: {chromosome}")"""

    def analyze(self):
        self.build_variant_maps()
        self.compare_reads_to_variants()
        self.generate_report()    

    def compare_reads_to_variants(self):
        logger.info("Start comparing reads to the variant map...")
        with multiprocessing.get_context("spawn").Pool() as pool:
            results = []
            for chromosome in self.dict_ref_haps.keys():
                results.append(pool.apply_async(VariantAnalyzer.process_chromosome, (
                    chromosome, self.f_name, self.dict_ref_haps[chromosome], 
                    self.dict_ref_cohorts[chromosome], self.dict_set_conflict_vars[chromosome], 
                    self.dict_ref_var_pos[chromosome], self.real_data, self.run_id
                )))
            #pool.close()
            #pool.join()
            for result in results:
                chromosome_bias = result.get()
                self.dict_ref_bias.update(chromosome_bias)

    def build_variant_maps(self):
        logger.info("Start building the variant maps...")
        self.variant_seq()
        
        # Extend conflict set
        for ref_name, positions in self.dict_set_conflict_vars.items():
            for pos in list(positions):
                self.dict_set_conflict_vars[ref_name].update(range(pos - VAR_CHAIN, pos + VAR_CHAIN))

    def variant_seq(self):
        for ref_name in self.f_fasta.references:
            list_var = list(self.f_vcf.fetch(ref_name))
            self.dict_ref_var_pos[ref_name] = [[var.start for var in list_var],[var.stop for var in list_var]]
            len_list_var = len(list_var)
            idx_var = 0
            while idx_var < len_list_var:
                idx_var = self.process_variant(ref_name, list_var, idx_var)

    def process_variant(self, ref_name: str, list_var: List[pysam.VariantRecord], idx_vcf: int):
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


    def generate_report(self):
        logger.info("Start output report...")
        self.output_report()

    def output_report(self):
        with open(self.output_file, 'w') as f_all, \
             open(self.output_file + '.gap', 'w') as f_gap, \
             open(self.output_file + '.SNP', 'w') as f_SNP:
            
            headers = self.get_report_headers()
            f_all.write(headers['all'])
            f_gap.write(headers['gap'])
            f_SNP.write(headers['SNP'])

            self.f_vcf.reset()
            for var in self.f_vcf:
                self.write_variant_data(var, f_all, f_gap, f_SNP)

    def get_report_headers(self) -> Dict[str, str]:
        if self.real_data:
            header = "CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tBOTH\tOTHER\tGAP\n"
            return {'all': header, 'gap': header, 'SNP': header}
        else:
            header = "CHR\tHET_SITE\tNUM_READS\tAVG_MAPQ\tEVEN_P_VALUE\tBALANCE\tREF\tALT\tBOTH\tOTHER\tMAP_BALANCE\tMAP_REF\tMAP_ALT\tMIS_MAP\tSIM_BALANCE\tSIM_REF\tSIM_ALT\tGAP\n"
            return {'all': header, 'gap': header, 'SNP': header}

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
        np.add.at(list_count, input_idx, 1)
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


def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze reference bias in genomic sequencing data.")
    parser.add_argument('-v', '--vcf', help='VCF file', required=True)
    parser.add_argument('-s', '--sam', help='SAM file', required=True)
    parser.add_argument('-f', '--fasta', help='Reference FASTA file', required=True)
    parser.add_argument('-r', '--real_data', help='Turn off hap_information warning for real data', action='store_true')
    parser.add_argument('-p', '--golden_pickle', help='Pickle file containing the golden information for report reference')
    parser.add_argument('-t', '--run_id', help='Tag for run_id, can be used to indicate chromosome number')
    parser.add_argument('-o', '--out', help='Output file', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    try:
        args = parse_arguments()
        analyzer = VariantAnalyzer(args.vcf, args.sam, args.fasta, 
                                   args.golden_pickle, args.run_id, 
                                   args.real_data, args.out)
        analyzer.analyze()

        #p = pstats.Stats('profile_output.prof')
        #p.sort_stats('cumulative').print_stats(20)
        #p.print_callers('acquire')

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise