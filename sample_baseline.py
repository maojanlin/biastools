import pysam
import argparse
import random
import numpy as np
import os
from subprocess import call

from scanning_bias import scanning_bias 


def baseline(
        f_mpileup   :pysam.VariantRecord
        ) -> tuple:
    """
    Take in the sample mpileup, and output the average read_depth/variant Density/non Diploid portion
    """
    dict_3D_measures = scanning_bias(
        f_gvcf=f_mpileup
        )
    
    total_read_depth  = np.array([])
    total_var_density = np.array([])
    total_dip_density = np.array([])
    print('read_depth(mean/std)   var_density   non_diploid')
    for ref_name, dict_array in dict_3D_measures.items():
        for start_pos, array_info in dict_array.items():
            array_read_depth, array_var_density, array_dip_density = array_info
            
            avg_read_depth   = np.mean(array_read_depth)
            std_read_depth   = np.std(array_read_depth)
            #positive_avg_var = np.mean(array_var_density)
            #positive_std_var = np.std(array_var_density)
            #positive_avg_dip = np.mean(array_dip_density)
            #positive_std_dip = np.std(array_dip_density)
            positive_var     = array_var_density[array_var_density != 0]
            positive_dip     = array_dip_density[array_var_density != 0]
            positive_avg_var = np.mean(positive_var)
            positive_std_var = np.std(positive_var)
            positive_avg_dip = np.mean(positive_dip)
            positive_std_dip = np.std(positive_dip)

            print(ref_name, round(avg_read_depth,2), round(std_read_depth,2), '     ', round(positive_avg_var,2), round(positive_std_var,2), '     ', round(positive_avg_dip,2), round(positive_std_dip, 2))
            total_read_depth  = np.concatenate((total_read_depth , array_read_depth))
            total_var_density = np.concatenate((total_var_density, positive_var))
            total_dip_density = np.concatenate((total_dip_density, positive_dip))
    
    print(len(total_read_depth))
    print('total_read_depth (mean/std)')
    print(round(np.mean(total_read_depth),2), round(np.std(total_read_depth),2))
    print('total_var_density (mean/std)')
    print(round(np.mean(total_var_density),2), round(np.std(total_var_density),2))
    print('total_non_diploid (mean/std)')
    print(round(np.mean(total_dip_density),2), round(np.std(total_dip_density),2))
    return np.mean(total_read_depth), np.mean(total_var_density), np.mean(total_dip_density)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam_file', help='the bam file we want to sample')
    parser.add_argument('-f', '--reference_fasta', help='the reference fasta file for mpileup building')
    parser.add_argument('-o', '--sample_bed', help='the sampled 1/1000 bed file')
    parser.add_argument('-th', '--threshold_contig', help='the minimum contig length for sampling', type=int, default=10000000)
    parser.add_argument('--seed', help='seed for random sampling', type=int, default=0)
    parser.add_argument('-k', '--kill', help='kill all storage files', action='store_true')
    args = parser.parse_args()
    
    fn_bam    = args.bam_file
    fn_ref    = args.reference_fasta
    fn_sample = args.sample_bed
    min_len   = args.threshold_contig
    seed      = args.seed
    kill_flag = args.kill

    f_bam = pysam.AlignmentFile(fn_bam)

    random.seed(seed)
    fo = open(fn_sample + '.bed', 'w')
    for idx, name in enumerate(f_bam.header.references):
        contig_len = f_bam.header.lengths[idx]
        if contig_len > min_len:
            thousandth = int(contig_len / 1000)
            sample_start = random.randint(0, contig_len - thousandth)
            fo.write(name + ' ' + str(sample_start) + ' ' + str(sample_start+thousandth) + '\n')
    fo.close()

    # SAMTOOLS command for extract the sample region bam file
    if os.path.exists(fn_sample + '.bam') and not kill_flag:
        print(fn_sample + '.bam already exist.')
    else:
        command = ('samtools view -h ' + fn_bam + ' -L ' + fn_sample + '.bed -o ' + fn_sample + '.bam')
        print(command)
        call(command, shell=True)

    # BCFTOOLS command for mpileup the bam file
    if os.path.exists(fn_sample + '.mpileup') and not kill_flag:
        print(fn_sample + '.mpileup already exist.')
    else:
        command = ('bcftools mpileup --count-orphans --annotate FORMAT/AD,FORMAT/DP -f ' + fn_ref + ' --min-BQ 0 --min-MQ 0 ' \
                   + fn_sample + '.bam -o ' + fn_sample + '.mpileup')
        print(command)
        call(command, shell=True)

    f_mpileup = pysam.VariantFile(fn_sample + '.mpileup')
    baseline(f_mpileup)

    
