import pysam
import gzip

import numpy as np
import os
import argparse
import pickle


def output_wig(
        output_name :str,
        ref_name    :str,
        data_name   :str,
        wig_start   :int,
        wig_end     :int,
        array_wig   :np.array
        ) -> None:
    """
    output single wig file
    """
    f_o = gzip.open(output_name, 'wt')
    f_o.write("browser position " + ref_name + ":" + str(wig_start) + "-" + str(wig_end) + '\n')
    f_o.write("browser hide all\n")
    f_o.write("track type=wiggle_0 name=\"" + data_name + "\" description=\"variableStep format\"  visibility=hide autoScale=on" + \
            "color=50,150,255 graphType=points priority=10\n")
    f_o.write("variableStep chrom=" + ref_name + '\n')
    for idx, depth in enumerate(array_wig):
        f_o.write(str(wig_start+idx) + ' ' + str(round(depth, 2)) + '\n')
    f_o.close()



def report_wig(
        fn_output   :str,
        wig_info    :tuple,
        ) -> None:
    """
    output the wig format for read_depth, var_density, and dip_density
    this whole process take times
    """
    ref_name, region_begin, array_read_depth, array_var_density, array_dip_density, array_score, array_score_sum = wig_info
    output_wig(
        output_name = (fn_output + '.read_depth.wig.gz'),
        ref_name    = ref_name,
        data_name   = 'avg_read_depth',
        wig_start   = region_begin,
        wig_end     = region_begin+len(array_read_depth),
        array_wig   = array_read_depth
        )
    output_wig(
        output_name = (fn_output + '.var_density.wig.gz'),
        ref_name    = ref_name,
        data_name   = 'var_density',
        wig_start   = region_begin,
        wig_end     = region_begin+len(array_var_density),
        array_wig   = array_var_density
        )
    output_wig(
        output_name = (fn_output + '.dip_density.wig.gz'),
        ref_name    = ref_name,
        data_name   = 'non_diploid_density',
        wig_start   = region_begin,
        wig_end     = region_begin+len(array_dip_density),
        array_wig   = array_dip_density
        )
    output_wig(
        output_name = (fn_output + '.score.wig.gz'),
        ref_name    = ref_name,
        data_name   = '3D_scoring_product',
        wig_start   = region_begin,
        wig_end     = region_begin+len(array_score),
        array_wig   = array_score
        )
    output_wig(
        output_name = (fn_output + '.score_sum.wig.gz'),
        ref_name    = ref_name,
        data_name   = '3D_scoring_sum',
        wig_start   = region_begin,
        wig_end     = region_begin+len(array_score_sum),
        array_wig   = array_score_sum
        )


def count_density(
        list_var_sites  :list,
        window_size     :int=400
        ) -> list:
    """
    count how many variants and non-diploid variants are nearby (with in window size)
    """
    half_window = round(window_size/2)
    density_var = np.zeros(len(list_var_sites))
    density_dip = np.zeros(len(list_var_sites))
    for idx, var_info in enumerate(list_var_sites):
        idy = idx + 1
        pos = var_info[1]
        while idy < len(list_var_sites):
            pos_check = list_var_sites[idy][1]
            if pos_check <= pos + half_window:
                # how many variants are nearby
                density_var[idx] += 1
                density_var[idy] += 1
                # how many non-diploid variants are nearby
                density_dip[idx] += list_var_sites[idy][3]
                density_dip[idy] += var_info[3]
            else:
                break
            idy += 1
    return density_var, density_dip


def score_dep_dip_den(
        list_var_sites  :list,
        density_var     :list,
        density_dip     :list,
        debug           :bool=False
        ) -> list:
    """
    count the final scoring based on the depth / density of var / proportion of non-diploid var
    """
    list_bias_vars = []
    for idx, var_info in enumerate(list_var_sites):
        # Depth score
        score = [var_info[-1][0],0,0]
        # Density score
        if density_var[idx] > 10: 
            score[2] = 1
        elif density_var[idx] > 3:
            score[2] = 0.5
        # Diploid score
        if density_dip[idx]*2 > density_var[idx]:
            score[1] = 1
        elif density_dip[idx]*10 > density_var[idx]*3:
            score[1] =0.5
        # final score
        list_bias_vars.append(sum(score) >= 2)
        if debug:
            if sum(score) >= 2:
                print('*', round(density_dip[idx],1), density_var[idx], var_info)
            else:
                print(' ', round(density_dip[idx],1), density_var[idx], var_info)
    return list_bias_vars


def scanning_bias(
        f_gvcf      :pysam.VariantRecord,
        window_size :int=400
        ) -> dict:
    """
    Scanning the fn_gvcf to find the region with 
        - high read depth,
        - high density of variants, or
        - non diploid evidence.
    """
    # Extract the read_depth and variant informations
    ref_name = None         # record the reference name
    list_depth = []         # record the read_depth per site
    list_var_sites = []     # record the variant information per var
    dict_ref_info  = {}
    for var in f_gvcf:
        ref_name = var.contig
        if dict_ref_info.get(ref_name) == None:
            dict_ref_info[ref_name] = {'depth':[], 'var':[]}
        total_depth = var.samples[0]['DP']
        
        # store the read depth
        dict_ref_info[ref_name]['depth'].append((var.start, total_depth))

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
            if depth >= total_depth*15/100: # consider as variant
                num_var = idx + 1
            else:
                break
        if num_var > 1:
            nonDip_flag = False
            if num_var > 2 or list_alt_depth[1]*2 <= list_alt_depth[0]:
                nonDip_flag = True                                                              # -> for debug purpose
            dict_ref_info[ref_name]['var'].append([var.start, total_depth, list_alt_depth[:num_var], nonDip_flag, alt_depth, list_alleles])

    # Two parameters we have:
    # list_depth
    # list_var_sites

    # Analyze the density of the variants
    dict_3D_measures = {}
    for ref_name, dict_var_info in dict_ref_info.items():
        list_depth     = dict_var_info['depth']
        list_var_sites = dict_var_info['var']
        density_var, density_dip = count_density(list_var_sites, window_size)

        half_window = round(window_size/2)
        # Counting average readepth over the window
        region_begin = list_depth[0][0]
        region_end   = list_depth[-1][0]
        array_read_depth = np.zeros(region_end - region_begin + window_size)
        for site_info in list_depth:
            index = site_info[0] - region_begin
            depth = site_info[1]
            array_read_depth[index:index+window_size] += depth
        array_read_depth /= window_size
        array_read_depth = array_read_depth[half_window:-half_window]

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

        dict_3D_measures[ref_name] = {}
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
        if max(array_score[pos_start:pos_stop]) > threshold_2:
            list_bias.append((pos_start + region_begin, pos_stop + region_begin))
        else:
            list_suspicious.append((pos_start + region_begin, pos_stop + region_begin))
    if f_ob:
        for segment in list_bias:
            f_ob.write(ref_name + ' ' + str(segment[0]) + ' ' + str(segment[1]) + '\n')
    if f_os:
        for segment in list_suspicious:
            f_os.write(ref_name + ' ' + str(segment[0]) + ' ' + str(segment[1]) + '\n')
    return list_bias, list_suspicious



def calculate_3D_score(
        dict_3D_measures :dict,
        fn_out_report    :str
    ) -> tuple:
    """
    Take in the 3D measures and output the 3D score
    """
    f_ob = open(fn_out_report + '.bias.bed', 'w')
    f_os = open(fn_out_report + '.suspicious.bed', 'w')
    f_ob.write('#chrom chromStart chromEnd\n')
    f_os.write('#chrom chromStart chromEnd\n')
    for ref_name, dict_region_begin in dict_3D_measures.items():
        for region_begin, array_info in dict_region_begin.items():
            array_read_depth, array_var_density, array_dip_density = array_info
            
            # calculate the mean for normalization 
            #TODO should be replaced by biastools.baseline
            avg_read_depth   = np.mean(array_read_depth)
            positive_avg_var = sum(array_var_density*(array_var_density!=0))/sum(array_var_density!=0)
            positive_avg_dip = sum(array_dip_density*(array_dip_density!=0))/sum(array_dip_density!=0)
            
            #print(avg_read_depth, positive_avg_var, positive_avg_dip)
            array_score_product = np.round(array_read_depth/avg_read_depth) * (array_var_density/positive_avg_var+0.1) * (array_dip_density/positive_avg_dip+0.1)
            array_score_product = np.where(array_score_product > 30, 30, array_score_product)

            array_score_sum = np.round(array_read_depth/avg_read_depth) + (array_var_density/positive_avg_var) + (array_dip_density/positive_avg_dip)
            array_score_sum = np.where(array_score_sum > 30, 30, array_score_sum)
            link_bias_region_and_report(array_score_sum, region_begin, ref_name, f_ob, f_os)
    f_ob.close()
    f_os.close()

    #return array_score_product, array_score_sum
    return ref_name, region_begin, array_read_depth, array_var_density, array_dip_density, array_score_product, array_score_sum






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gvcf_file', help='the gvcf file of a specific region')
    parser.add_argument('-w', '--window_size', help='window size for average depth, density analysis', type=int, default=400)
    parser.add_argument('-rd', '--read_depth', help='the average sequence read depth', type=int, default=30)
    parser.add_argument('-o',  '--out_report', help='scanning bed file and reports')
    parser.add_argument('-wig', '--out_wig', help='flag for wig output', action='store_true')
    args = parser.parse_args()
    
    fn_gvcf       = args.gvcf_file
    rd_thresh     = args.read_depth
    window_size   = args.window_size
    fn_out_report = args.out_report
    flag_wig      = args.out_wig

    f_gvcf = pysam.VariantFile(fn_gvcf)
    # load or calculate the 3D measures depending on pickle file existance
    if os.path.exists(fn_out_report + '.pickle'):
        print("Pickle file", fn_out_report + '.pickle', 'exist, load it instead of recalculate...')
        f_i = open(fn_out_report + '.pickle', 'rb')
        dict_3D_measures = pickle.load(f_i)
        f_i.close()
    else:
        print("Process the mpileup file", fn_gvcf + '...')
        dict_3D_measures = scanning_bias(
            f_gvcf=f_gvcf,
            window_size=window_size
            )
        print("Store the measures information as", fn_out_report + '.pickle...')
        f_o = open(fn_out_report + '.pickle', 'wb')
        pickle.dump(dict_3D_measures, f_o)
        f_o.close()
    
    print("Calculate 3D scoring and output bed...")
    tuple_3dScore = calculate_3D_score(dict_3D_measures, fn_out_report)

    print("Output wig format...")
    if flag_wig: # output wig files if -ow option
        report_wig(
            fn_output=fn_out_report,
            wig_info=tuple_3dScore
            )
