import argparse
import pysam
import numpy as np


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
        #genotype = var.samples[0]['GT']
        total_depth = var.samples[0]['DP']
        alt_depth = None
        if var.samples[0].get('AD'):
            alt_depth = list(var.samples[0]['AD'])
        else:
            alt_depth = [0, total_depth]

        dep_dip_den = [0,0,0] # score of read_depth, deploid, variant density
        # calculate depth score
        if total_depth > 3*rd_thresh: # high read depth
            dep_dip_den[0] = 1
        elif total_depth > 2*rd_thresh: # slightly high read depth
            dep_dip_den[0] = 0.5
        
        # calculate diploid score
        list_alleles = list(var.alleles)
        if sum(alt_depth) != total_depth:
            alt_depth.append(total_depth - sum(alt_depth))
            list_alleles.append('Others')
        list_alt_depth = sorted(alt_depth, reverse=True)
        max_alt_depth = list_alt_depth[0]
        for depth in list_alt_depth[1:]: 
            if depth*5 >= max_alt_depth: # consider as variant
                if depth*2 <= max_alt_depth: # stronly not diploid
                    dep_dip_den[1] = 1
                elif depth*3 <= max_alt_depth*2: # slightly not diploid
                    dep_dip_den[1] = 0.5
                else: # non diploid
                    dep_dip_den[1] = 0.1
        if dep_dip_den[1] > 0:
            list_var_sites.append([ref_name, var.start, total_depth, alt_depth, list_alleles, dep_dip_den])

    # Analyze the density of the variants
    density_map = np.zeros(len(list_var_sites))
    for idx, var_info in enumerate(list_var_sites):
        idy = idx + 1
        pos = var_info[1]
        while idy < len(list_var_sites):
            pos_check = list_var_sites[idy][1]
            if pos_check <= pos + 200:
                density_map[idx] += 1
                density_map[idy] += 1
            else:
                break
            idy += 1
        if density_map[idx] > 10:
            list_var_sites[-1][2] = 1
        elif density_map[idx] > 3:
            list_var_sites[-1][2] = 0.5

    # print the result
    for idx, var_info in enumerate(list_var_sites):
        # calculate score`/
        dep_dip_den = var_info[-1]
        if sum(dep_dip_den) >= 1.5:
            print('*', density_map[idx], var_info)
        else:
            print(' ', density_map[idx], var_info)

    print("number of sites being print:", len(list_var_sites))
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
