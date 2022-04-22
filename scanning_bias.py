import argparse
import pysam
import numpy as np


def count_density(
        list_var_sites  :list
        ) -> list:
    """
    count how many variants and non-diploid variants are nearby
    """
    density_var = np.zeros(len(list_var_sites))
    density_dip = np.zeros(len(list_var_sites))
    for idx, var_info in enumerate(list_var_sites):
        idy = idx + 1
        pos = var_info[1]
        while idy < len(list_var_sites):
            pos_check = list_var_sites[idy][1]
            if pos_check <= pos + 200:
                # how many variants are nearby
                density_var[idx] += 1
                density_var[idy] += 1
                # how many non-diploid variants are nearby
                density_dip[idx] += list_var_sites[idy][-1][1]
                density_dip[idy] += var_info[-1][1]
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
        rd_thresh   :int
        ) -> dict:
    """
    Scanning the fn_gvcf to find the region with 
        - high read depth,
        - high density of variants, or
        - non diploid evidence.
    """
    # Scan the mpileup file to find all the potential variant sites
    list_var_sites = []
    for var in f_gvcf:
        ref_name = var.contig
        total_depth = var.samples[0]['DP']
        alt_depth = None
        if var.samples[0].get('AD'):
            alt_depth = list(var.samples[0]['AD'])
        else:
            alt_depth = [0, total_depth]

        dep_dip_plod = [0,0,1] # score of read_depth, deploid, ploid
        # calculate depth score
        if total_depth > 3*rd_thresh: # high read depth
            dep_dip_plod[0] = 1
        elif total_depth > 2*rd_thresh: # slightly high read depth
            dep_dip_plod[0] = 0.5
        
        # calculate diploid score
        list_alleles = list(var.alleles)
        if sum(alt_depth) != total_depth: # often happens at indels
            alt_depth.append(total_depth - sum(alt_depth))
            list_alleles.append('Others')
        list_alt_depth = sorted(alt_depth, reverse=True)
        max_alt_depth = list_alt_depth[0]
        for idx, depth in enumerate(list_alt_depth[1:]):
            if depth*5 >= max_alt_depth: # consider as variant
                dep_dip_plod[2] = idx + 1
                if idx >= 1: # more than diploid
                    dep_dip_plod[1] = 1
                elif depth*2 <= max_alt_depth: # stronly not diploid
                    dep_dip_plod[1] = 1
                elif depth*3 <= max_alt_depth*2: # slightly not diploid
                    dep_dip_plod[1] = 0.5
                else: # diploid
                    dep_dip_plod[1] = 0.01
        if dep_dip_plod[1] > 0:
            list_var_sites.append([ref_name, var.start, total_depth, alt_depth, list_alleles, dep_dip_plod])
            
    # Analyze the density of the variants
    density_var, density_dip = count_density(list_var_sites)

    # calculate bias score`/
    list_bias_vars = score_dep_dip_den(list_var_sites, density_var, density_dip)
    
    list_segment = []
    start_pos = list_var_sites[0][1]
    end_pos   = list_var_sites[0][1]
    num_var = 1
    num_dip = 1 if (list_var_sites[0][-1][1] >= 0.5) else 0
    tmp_num_var, tmp_num_dip = 0, 0
    for idx, var_info in enumerate(list_var_sites[1:]):
        ref_name = var_info[0]
        site_pos = var_info[1]
        dep_dip_plod = var_info[-1]
        non_diploid  = 1 if (dep_dip_plod[1] >= 0.5) else 0
        if list_bias_vars[idx] == True:
            if site_pos -end_pos < 3000: # segment connection distance
                end_pos = site_pos
                num_var += tmp_num_var + 1
                num_dip += tmp_num_dip + non_diploid
            else:
                list_segment.append((ref_name, start_pos, end_pos, num_var, num_dip))
                start_pos = site_pos
                end_pos   = site_pos
                num_var   = 1
                num_dip   = non_diploid
            tmp_num_var = 0
            tmp_num_dip = 0
        else:
            tmp_num_var += 1
            tmp_num_dip += non_diploid
    if len(list_segment) == 0 or list_segment[-1][1] != start_pos:
        list_segment.append((ref_name, start_pos, end_pos, num_var, num_dip))

    # output the segments
    list_suspicious = []
    print("The bias region:")
    print("Ref\tRegion\t\t\tLength\tNum_Var\tNum_non-diplod\tVar_per_1000\tNon-diploid(%)")
    for segment in list_segment:
        if segment[2] - segment[1] + 1 > 10000:
            len_segment = segment[2] - segment[1] + 1
            print(segment[0], str(segment[1]) + '-' + str(segment[2]), len_segment, segment[3], segment[4], round(segment[3]/len_segment*1000), round(segment[4]/segment[3]*100,1), sep='\t')
        else:
            list_suspicious.append(segment)
    print("The suspicous region:")
    for segment in list_suspicious:
        len_segment = segment[2] - segment[1] + 1
        if segment[3] <= 3:
            print(segment[0], str(segment[1]) + '-' + str(segment[2]), len_segment, segment[3], segment[4], sep='\t')
        else:
            print(segment[0], str(segment[1]) + '-' + str(segment[2]), len_segment, segment[3], segment[4], round(segment[3]/len_segment*1000), round(segment[4]/segment[3]*100,1), sep='\t')
    print("Number of variantst:", len(list_var_sites))
    print("Number of non-diploid variants:", sum([(var_info[-1][1] >= 0.5) for var_info in list_var_sites]))



    return None





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gvcf_file', help='the gvcf file of a specific region')
    parser.add_argument('-rd', '--read_depth', help='the average sequence read depth', type=int, default=30)
    parser.add_argument('-o', '--out', help='scanning report')
    args = parser.parse_args()
    
    fn_gvcf   = args.gvcf_file
    rd_thresh = args.read_depth
    fn_output = args.out

    f_gvcf = pysam.VariantFile(fn_gvcf)
    scanning_bias(
        f_gvcf=f_gvcf,
        rd_thresh=rd_thresh
        )
