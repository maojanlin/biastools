import pysam
import gzip

import numpy as np
import os
import argparse
import pickle


def output_wig(
        output_name :str,
        data_name   :str,
        list_data   :list
        ) -> None:
    """
    output single wig file
    """
    f_o = gzip.open(output_name, 'wt')
    for array_info in list_data:
        ref_name, wig_start, array_wig = array_info
        wig_end = wig_start + len(array_wig)
        
        f_o.write("browser position " + ref_name + ":" + str(wig_start) + "-" + str(wig_end) + '\n')
        f_o.write("browser hide all\n")
        f_o.write("track type=wiggle_0 name=\"" + data_name + "\" description=\"variableStep format\"  visibility=hide autoScale=on" + \
                "color=50,150,255 graphType=points priority=10\n")
        f_o.write("variableStep chrom=" + ref_name + '\n')
        for idx, depth in enumerate(array_wig):
            f_o.write(str(wig_start+idx) + ' ' + str(round(depth, 2)) + '\n')
    f_o.close()



def report_wig(
        fn_output           :str,
        dict_3D_measures    :dict,
        ) -> None:
    """
    output the wig format for read_depth, var_density, and dip_density
    this whole process take times
    """
    # list_info is composed of array_RD, array_VD, array_ND, array_score
    # (ref_name, region_begin, array_info)
    list_info = [[],[],[],[]]
    #ref_name, region_begin, array_read_depth, array_var_density, array_dip_density, array_score, array_score_sum = wig_info
    for ref_name, dict_array in dict_3D_measures.items():
        for start_pos, array_info in dict_array.items():
            array_RD, array_VD, array_ND, array_score = array_info
            list_info[0].append((ref_name, start_pos, array_RD))
            list_info[1].append((ref_name, start_pos, array_VD))
            list_info[2].append((ref_name, start_pos, array_ND))
            list_info[3].append((ref_name, start_pos, array_score))

    output_wig(
        output_name = (fn_output + '.read_depth.wig.gz'),
        data_name   = 'avg_read_depth',
        list_data   = list_info[0]
        )
    output_wig(
        output_name = (fn_output + '.var_density.wig.gz'),
        data_name   = 'var_density',
        list_data   = list_info[1]
        )
    output_wig(
        output_name = (fn_output + '.dip_density.wig.gz'),
        data_name   = 'non_diploid_density',
        list_data   = list_info[2]
        )
    output_wig(
        output_name = (fn_output + '.score_sum.wig.gz'),
        data_name   = '3D_scoring_sum',
        list_data   = list_info[3]
        )


def scanning_bias(
        f_gvcf      :pysam.VariantRecord
        ) -> dict:
    """
    Scanning the fn_gvcf to find the region with 
        - high read depth,
        - high density of variants, or
        - non diploid evidence.
    return the raw numbers
    """
    # Extract the read_depth and variant informations
    ref_name = None # record the reference name
    last_pos = -2   # record the last mpileup position
    start_pos = None # record the starting position of each region
    dict_ref_info  = {}
    for var in f_gvcf:
        if ref_name != var.contig: # new chromosome
            ref_name = var.contig
            dict_ref_info[ref_name] = {}

            start_pos = var.start
            dict_ref_info[ref_name][start_pos] = {'depth':[], 'var':[]}
        elif var.start > last_pos + 1: # the same chromsome, new position
            start_pos = var.start
            dict_ref_info[ref_name][start_pos] = {'depth':[], 'var':[]}
        elif var.start == last_pos: # duplicate position, pop the last read depth info
            dict_ref_info[ref_name][start_pos]['depth'].pop() 
        last_pos = var.start
        
        ref_name = var.contig
        total_depth = var.samples[0]['DP']
        
        # store the read depth
        dict_ref_info[ref_name][start_pos]['depth'].append((var.start, total_depth))

        alt_depth = None
        if var.samples[0].get('AD'):
            alt_depth = list(var.samples[0]['AD'])
        else:
            alt_depth = [0, total_depth]
        
        # calculate diploid score
        list_alleles = list(var.alleles)
        if sum(alt_depth) != total_depth: # often happens at indels
            alt_depth.append(total_depth - sum(alt_depth))
            list_alleles.append('Others')
        list_alt_depth = sorted(alt_depth, reverse=True)
        #max_alt_depth = list_alt_depth[0]
        num_var = 0
        for idx, depth in enumerate(list_alt_depth):
            if depth > total_depth*15/100: # consider as variant, exclude the 0,0 case
                num_var = idx + 1
            else:
                break
        if num_var > 1:
            nonDip_flag = False
            if num_var > 2 or list_alt_depth[1]*2 < list_alt_depth[0]:
                nonDip_flag = True                                                                                           
            dict_ref_info[ref_name][start_pos]['var'].append([var.start, total_depth, list_alt_depth[:num_var], nonDip_flag, \
                                                              alt_depth, list_alleles])
                                                              # -> for debug purpose
    return dict_ref_info


def boundary_compensate(
        target_array    :np.array,
        window_size     :int
    ) -> np.array:
    """
    compensate for padding zeros
    """
    if len(target_array) < window_size:
        return target_array

    half_window = int(window_size/2)
    # compensate left side
    for idx in range(half_window):
        target_array[idx] *= (window_size / (half_window+idx))
    # compensate right side
    for idx in range(-1, -half_window-1, -1):
        target_array[idx] *= (window_size / (half_window-idx-1))
    return target_array


def calculate_measures(
        dict_ref_info   :dict,
        window_size     :int=400
    ) -> dict:
    """
    Take the raw data and calculate 
        - the moving average of read_depth
        - over window number of variants
        - over window number of non_diploid site
    """
    # Two parameters we have:
    # list_depth
    # list_var_sites

    # Analyze the density of the variants
    dict_3D_measures = {}
    for ref_name, dict_start_pos in dict_ref_info.items():
        dict_3D_measures[ref_name] = {}
        for start_pos, dict_var_info in dict_start_pos.items():
            list_depth     = dict_var_info['depth']
            list_var_sites = dict_var_info['var']

            half_window = round(window_size/2)
            # Counting average readepth over the window (moving average)
            region_begin = list_depth[0][0]
            region_end   = list_depth[-1][0] + 1
            assert(start_pos == region_begin)
            array_read_depth = np.zeros(region_end - region_begin + window_size)
            for site_info in list_depth:
                index = site_info[0] - region_begin
                depth = site_info[1]
                array_read_depth[index:index+window_size] += depth
            array_read_depth /= window_size
            array_read_depth = array_read_depth[half_window:-half_window]
            array_read_depth = boundary_compensate(array_read_depth, window_size)

            # Calculate variant density over the window
            array_var_density = np.zeros(region_end - region_begin + window_size)
            array_dip_density = np.zeros(region_end - region_begin + window_size)
            for site_info in list_var_sites:
                index = site_info[0] - region_begin
                nonDip_flag = site_info[3]
            
                array_var_density[index:index+window_size] += 1
                if nonDip_flag:
                    array_dip_density[index:index+window_size] += 1
            array_var_density = array_var_density[half_window:-half_window]
            array_dip_density = array_dip_density[half_window:-half_window]
            #array_var_density = boundary_compensate(array_var_density, window_size)
            #array_dip_density = boundary_compensate(array_dip_density, window_size)

            dict_3D_measures[ref_name][region_begin] = [array_read_depth, array_var_density, array_dip_density]
    return dict_3D_measures


def link_bias_region_and_report(
        array_score     :np.array,
        region_begin    :int,
        ref_name        :str,
        f_ob            ,
        f_os            ,
        threshold_1     :int=3,
        threshold_2     :int=5,
        link_dist       :int=1000
    ) -> tuple:
    """
    Find and link the bias region according to thresholds
    report files:
        - bed file: bias region
        - bed file: suspicious region
        - csv file: detailed report of bias and suspicious region
    """
    list_region = []
    pos_start = -1
    pos_stop  = -link_dist -1
    for idx, score in enumerate(array_score):
        if score > threshold_1:
            if idx > pos_stop + link_dist:
                list_region.append((pos_start, pos_stop+1))
                #print(idx, pos_start, pos_stop+1)
                pos_start = idx
                pos_stop  = idx
            else:
                pos_stop = idx
    if len(list_region) == 0 or list_region[-1] != (pos_start, pos_stop+1):
        list_region.append((pos_start, pos_stop+1))
    list_region = list_region[1:] # first region is decoy
    
    # report bias region and suspicious region
    list_bias = []
    list_suspicious = []
    for pos_start, pos_stop in list_region:
        max_score = max(array_score[pos_start:pos_stop])
        avg_score = np.mean(array_score[pos_start:pos_stop])
        if max_score > threshold_2:
            list_bias.append((pos_start + region_begin, pos_stop + region_begin, max_score, avg_score))
        else:
            list_suspicious.append((pos_start + region_begin, pos_stop + region_begin, max_score, avg_score))
    if f_ob:
        for segment in list_bias:
            f_ob.write(ref_name + '\t' + str(segment[0]) + '\t' + str(segment[1]) + '\tlen:' + str(segment[1]-segment[0]) + ',max:' + str(round(segment[2],2)) + ',avg:' + str(round(segment[3],2)) + '\n')
    if f_os:
        for segment in list_suspicious:
            f_os.write(ref_name + '\t' + str(segment[0]) + '\t' + str(segment[1]) + '\tlen:' + str(segment[1]-segment[0]) + ',max:' + str(round(segment[2],2)) + ',avg:' + str(round(segment[3],2)) + '\n')
    for segment in sorted(list_bias, key=lambda ele: (ele[1]-ele[0])*ele[2]*ele[3], reverse=True)[:5]:
        print(ref_name + ' ' + str(segment[0]) + ' ' + str(segment[1]) + ' len:' + str(segment[1]-segment[0]) + ',max:' + str(round(segment[2],2)) + ',avg:' + str(round(segment[3],2)))
        pass
    return list_bias, list_suspicious


def calculate_3D_score(
        dict_3D_measures :dict,
        fn_out_report    :str,
        list_statistics  :list
    ) -> tuple:
    """
    Take in the 3D measures and output the 3D score
    """
    avg_RD, std_RD, avg_VD, std_VD, avg_ND, std_ND = list_statistics

    f_ob = open(fn_out_report + '.bias.bed', 'w')
    f_os = open(fn_out_report + '.suspicious.bed', 'w')
    f_ob.write('#chrom\tchromStart\tchromEnd\tname\n')
    f_os.write('#chrom\tchromStart\tchromEnd\tname\n')

    link_dist = 1000
    for ref_name, dict_region_begin in dict_3D_measures.items():
        old_region_begin = -link_dist
        old_array = np.array([])
        for region_begin, array_info in sorted(dict_region_begin.items()):
            array_read_depth, array_var_density, array_dip_density = array_info
            #print(region_begin, region_begin+len(array_info[0]))

            #array_score_product = np.round(array_read_depth/avg_RD) * (array_var_density/avg_VD+0.1) * (array_dip_density/avg_ND+0.1)
            #array_score_product = np.where(array_score_product > 30, 30, array_score_product)
            """
            array_score_sum = np.round(array_read_depth/avg_RD) + (array_var_density/avg_VD) + (array_dip_density/avg_ND)
            array_score_sum = np.where(array_read_depth > avg_RD/2, array_score_sum, 0)
            """
            array_Z_score_RD = (array_read_depth-avg_RD)/std_RD - 1
            array_Z_score_RD = np.where(array_Z_score_RD > 0, array_Z_score_RD, 0)
            array_Z_score_VD = (array_var_density-avg_VD)/std_VD
            array_Z_score_VD = np.where(array_Z_score_VD > 0, array_Z_score_VD, 0)
            array_Z_score_ND = (array_dip_density-avg_ND)/std_ND
            array_Z_score_ND = np.where(array_Z_score_ND > 0, array_Z_score_ND, 0)
            array_score_sum  = array_Z_score_RD + array_Z_score_VD + array_Z_score_ND
            array_score_product  = array_Z_score_RD * (array_Z_score_VD + array_Z_score_ND)
            #array_score_sum = (array_read_depth-avg_RD)/std_RD
            #array_score_sum = (array_var_density-avg_VD)/std_VD
            #array_score_sum = (array_dip_density-avg_ND)/std_ND
            #array_score_sum = array_read_depth/avg_RD
            #array_score_sum = array_var_density/avg_VD

            #array_score_sum = np.where(array_score_sum > 0, array_score_sum, 0)
            #array_score_sum = np.where(array_score_sum > 30, 30, array_score_sum)
            #link_bias_region_and_report(array_score_sum, region_begin, ref_name, f_ob, f_os)
            #link_bias_region_and_report(array_score_product, region_begin, ref_name, f_ob, f_os,20,30,1000)
            #link_bias_region_and_report(array_score_sum, region_begin, ref_name, f_ob, f_os,3,5,link_dist)
            #print(old_region_begin, old_region_begin+len(old_array), region_begin)
            dict_3D_measures[ref_name][region_begin].append(array_score_sum)
            if old_region_begin + len(old_array) + link_dist > region_begin:
                assert(old_region_begin + len(old_array) < region_begin)
                # Connect
                diff = region_begin - old_region_begin - len(old_array)
                old_array = np.concatenate((old_array, np.zeros(diff), array_score_sum))
            else:
                if old_region_begin != -1000:
                    link_bias_region_and_report(old_array, old_region_begin, ref_name, f_ob, f_os,3,5,link_dist)
                #dict_3D_measures[ref_name][old_region_begin].append(old_array)
                old_region_begin = region_begin
                old_array        = array_score_sum
        link_bias_region_and_report(old_array, old_region_begin, ref_name, f_ob, f_os,3,5,link_dist)
    f_ob.close()
    f_os.close()
    
    # report the region with low Read depth
    f_or = open(fn_out_report + '.lowRd.bed', 'w')
    f_or.write('#chrom\tchromStart\tchromEnd\tname\n')
    rd_thresh = min(int(avg_RD/5),10)
    for ref_name, dict_region_begin in dict_3D_measures.items():
        global_start = []
        global_stop  = []
        for region_begin, array_info in sorted(dict_region_begin.items()):
            array_read_depth, *_ = array_info

            bool_low = array_read_depth < rd_thresh
            #print(bool_low)
            bool_low_shift = np.concatenate(([False], bool_low))[:-1]
            bool_start = bool_low > bool_low_shift
            bool_stop  = bool_low < bool_low_shift
            
            list_start = [idx+region_begin for idx, x in enumerate(bool_start) if x]
            list_stop  = [idx+region_begin for idx, x in enumerate(bool_stop ) if x]
            
            if len(list_start) == len(list_stop):
                list_start.append(region_begin + len(array_read_depth))
            if global_start == []:
                global_start = list_start
                global_stop  = list_stop
            else:
                if list_start[0] == region_begin:
                    global_start += list_start[1:]
                    global_stop  += list_stop
                else:
                    global_stop += [region_begin-1]
                    global_start += list_start
                    global_stop  += list_stop
        global_start = global_start[:-1]
        assert(len(global_start) == len(global_stop))
        for idx in range(len(global_start)):
            st = global_start[idx]
            ed = global_stop[idx]
            f_or.write(ref_name + '\t' + str(st) + '\t' + str(ed) + '\tlen:' + str(ed-st) + '\n')
    f_or.close()


def get_baseline(
        fn_baseline :str
    ) -> list:
    """
    Take and parse the last line of fn_baseline
    """
    f = open(fn_baseline, 'r')
    for line in f:
        pass
    f.close()
    _, avg_RD, std_RD, avg_VD, std_VD, avg_ND, std_ND = line.split()
    return [float(avg_RD), float(std_RD), float(avg_VD), float(std_VD), float(avg_ND), float(std_ND)]


def calculate_avg(
        dict_3D_measures :dict,
    ):
    total_read_depth  = np.array([])
    total_var_density = np.array([])
    total_dip_density = np.array([])
    for ref_name, dict_array in dict_3D_measures.items():
        for start_pos, array_info in dict_array.items():
            array_read_depth, array_var_density, array_dip_density = array_info
            positive_var      = array_var_density[array_var_density != 0]
            positive_dip      = array_dip_density[array_var_density != 0]
            
            total_read_depth  = np.concatenate((total_read_depth , array_read_depth))
            total_var_density = np.concatenate((total_var_density, positive_var))
            total_dip_density = np.concatenate((total_dip_density, positive_dip))
    return [np.mean(total_read_depth), np.std(total_read_depth), np.mean(total_var_density), \
            np.std(total_var_density), np.mean(total_dip_density), np.std(total_dip_density)]





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gvcf_file', help='the gvcf file of a specific region')
    parser.add_argument('-w', '--window_size', help='window size for average depth, density analysis', type=int, default=400)
    parser.add_argument('-rd', '--read_depth', help='the average sequence read depth')
    parser.add_argument('-b', '--baseline', help='the baseline report generate by sample_baseline.py')
    parser.add_argument('-s', '--sample', action='store_true', help='sample for the baseline')
    parser.add_argument('-o',  '--out_report', help='scanning bed file and reports')
    parser.add_argument('-wig', '--out_wig', help='flag for wig output', action='store_true')
    args = parser.parse_args()
    
    fn_gvcf       = args.gvcf_file
    rd_thresh     = args.read_depth
    window_size   = args.window_size
    fn_baseline   = args.baseline
    flag_sample   = args.sample
    fn_out_report = args.out_report
    flag_wig      = args.out_wig

    f_gvcf = pysam.VariantFile(fn_gvcf)
    # load or calculate the 3D measures depending on pickle file existance
    if os.path.exists(fn_gvcf + '.pickle'):
        print("Pickle file", fn_gvcf + '.pickle', 'exist, load it instead of recalculate...')
        f_i = open(fn_gvcf + '.pickle', 'rb')
        dict_3D_measures = pickle.load(f_i)
        f_i.close()
    else:
        print("Process the mpileup file", fn_gvcf + '...')
        dict_ref_info = scanning_bias(f_gvcf=f_gvcf)
        dict_3D_measures = calculate_measures(
            dict_ref_info=dict_ref_info,
            window_size=window_size
            )
        print("Store the measures information as", fn_gvcf + '.pickle...')
        f_o = open(fn_gvcf + '.pickle', 'wb')
        pickle.dump(dict_3D_measures, f_o)
        f_o.close()
    
    # Load or calculate the baseline of the measures
    if fn_baseline:
        # avg_RD, std_RD, avg_VD, std_VD, avg_ND, std_ND
        list_statistics = get_baseline(fn_baseline)
    elif flag_sample:
        list_statistics = calculate_avg(dict_3D_measures)
    else:
        list_statistics = [30, 10, 0.7, 1.6, 0.3, 1.2]
    if rd_thresh:
        list_statistics[0] = rd_thresh
    

    print("Calculate 3D scoring and output bed...")
    calculate_3D_score(dict_3D_measures, fn_out_report, list_statistics)

    if flag_wig: # output wig files if -ow option
        print("Output wig format...")
        report_wig(
            fn_output=fn_out_report,
            dict_3D_measures=dict_3D_measures
            )
