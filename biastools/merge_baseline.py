import pysam
import argparse
import random
import numpy as np
import os
from subprocess import call

from scanning_bias import scanning_bias, calculate_measures 
from sample_baseline import sample_select 
import pickle


def baseline(
        f_mpileup_1 :pysam.VariantRecord,
        f_mpileup_2 :pysam.VariantRecord,
        fn_sample   :str,
        window_size :int
        ) -> tuple:
    """
    Take in the sample mpileup, and output the average read_depth/variant Density/non Diploid portion
    """
    dict_ref_info_1 = scanning_bias(f_gvcf=f_mpileup_1)
    dict_3D_measures_1 = calculate_measures(
        dict_ref_info=dict_ref_info_1,
        window_size=window_size
        )
    dict_ref_info_2 = scanning_bias(f_gvcf=f_mpileup_2)
    dict_3D_measures_2 = calculate_measures(
        dict_ref_info=dict_ref_info_2,
        window_size=window_size
        )
    
    total_read_depth  = np.array([])
    total_var_density = np.array([])
    total_dip_density = np.array([])
    
    fo = open(fn_sample + '.baseline', 'w')
    fo.write('#chr pos segment_len RD_mean RD_std VD_mean VD_std ND_mean ND_std\n')
    for ref_name, dict_array in dict_3D_measures_1.items():
        for start_pos, array_info in dict_array.items():
            array_read_depth, array_var_density, array_dip_density = array_info
            
            avg_read_depth   = np.mean(array_read_depth)
            std_read_depth   = np.std(array_read_depth)

            fo.write(ref_name + ' ' + str(start_pos) + ' ' + str(len(array_read_depth)) + ' ')
            fo.write(str(round(avg_read_depth,2))   + ' ' + str(round(std_read_depth,2)) + ' ')
            #positive_var     = array_var_density[array_var_density != 0]
            positive_var     = array_var_density
            if len(positive_var) > 0:
                positive_avg_var = np.mean(positive_var)
                positive_std_var = np.std(positive_var)
                fo.write(str(round(positive_avg_var,2)) + ' ' + str(round(positive_std_var,2)) + ' ')
            
            #positive_dip     = array_dip_density[array_var_density != 0]
            positive_dip     = array_dip_density
            if len(positive_dip) > 0:
                positive_avg_dip = np.mean(positive_dip)
                positive_std_dip = np.std(positive_dip)
                fo.write(str(round(positive_avg_dip,2)) + ' ' + str( round(positive_std_dip, 2)) + '\n')
            else:
                fo.write('\n')
            
            total_read_depth  = np.concatenate((total_read_depth , array_read_depth))
            total_var_density = np.concatenate((total_var_density, positive_var))
            total_dip_density = np.concatenate((total_dip_density, positive_dip))
    for ref_name, dict_array in dict_3D_measures_2.items():
        for start_pos, array_info in dict_array.items():
            array_read_depth, array_var_density, array_dip_density = array_info
            
            avg_read_depth   = np.mean(array_read_depth)
            std_read_depth   = np.std(array_read_depth)

            fo.write(ref_name + ' ' + str(start_pos) + ' ' + str(len(array_read_depth)) + ' ')
            fo.write(str(round(avg_read_depth,2))   + ' ' + str(round(std_read_depth,2)) + ' ')
            #positive_var     = array_var_density[array_var_density != 0]
            positive_var     = array_var_density
            if len(positive_var) > 0:
                positive_avg_var = np.mean(positive_var)
                positive_std_var = np.std(positive_var)
                fo.write(str(round(positive_avg_var,2)) + ' ' + str(round(positive_std_var,2)) + ' ')
            
            #positive_dip     = array_dip_density[array_var_density != 0]
            positive_dip     = array_dip_density
            if len(positive_dip) > 0:
                positive_avg_dip = np.mean(positive_dip)
                positive_std_dip = np.std(positive_dip)
                fo.write(str(round(positive_avg_dip,2)) + ' ' + str( round(positive_std_dip, 2)) + '\n')
            else:
                fo.write('\n')
            
            total_read_depth  = np.concatenate((total_read_depth , array_read_depth))
            total_var_density = np.concatenate((total_var_density, positive_var))
            total_dip_density = np.concatenate((total_dip_density, positive_dip))

    fo.write('#total sample len: ' + str(len(total_read_depth)) + '\n')
    fo.write('#total_statistics:\n')
    fo.write('#chr pos segment_len RD_mean RD_std VD_mean VD_std ND_mean ND_std\n# ')
    fo.write(str(round(np.mean(total_read_depth),5))  + ' ' + str(round(np.std(total_read_depth),5)) + ' ')
    fo.write(str(round(np.mean(total_var_density),5)) + ' ' + str(round(np.std(total_var_density),5)) + ' ')
    fo.write(str(round(np.mean(total_dip_density),5)) + ' ' + str(round(np.std(total_dip_density),5)))
    fo.close()
    print("[Biastools] Generate " + fn_sample + '.baseline')
    return np.mean(total_read_depth), np.mean(total_var_density), np.mean(total_dip_density)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b1', '--bam_file_1', help='the bam file we want to sample')
    parser.add_argument('-b2', '--bam_file_2', help='the bam file we want to sample')
    parser.add_argument('-f', '--reference_fasta', help='the reference fasta file for mpileup building')
    parser.add_argument('-o', '--sample_bed', help='the sampled 1/1000 bed file')
    parser.add_argument('-w', '--window_size', help='window size for average depth, density analysis', type=int, default=400)
    parser.add_argument('-th', '--threshold_contig', help='the minimum contig length for sampling', type=int, default=10000000)
    parser.add_argument('--seed', help='seed for random sampling', type=int, default=0)
    args = parser.parse_args()
    
    fn_bam_1  = args.bam_file_1
    fn_bam_2  = args.bam_file_2
    fn_ref    = args.reference_fasta
    fn_sample = args.sample_bed
    min_len   = args.threshold_contig
    window_size = args.window_size
    seed      = args.seed

    # sample bed file according to the bam file information
    f_bam = pysam.AlignmentFile(fn_bam_1)
    sample_select(fn_sample, seed, min_len, f_bam)


    # SAMTOOLS command for extract the sample region bam file
    command = ('samtools view -h ' + fn_bam_1 + ' -L ' + fn_sample + '.bed -o ' + fn_bam_1 + '.sample.bam')
    print(command)
    call(command, shell=True)
    command = ('samtools view -h ' + fn_bam_2 + ' -L ' + fn_sample + '.bed -o ' + fn_bam_2 + '.sample.bam')
    print(command)
    call(command, shell=True)

    # BCFTOOLS command for mpileup the bam file
    command = ('bcftools mpileup --count-orphans --annotate FORMAT/AD,FORMAT/DP -f ' + fn_ref + ' --min-BQ 0 --min-MQ 0 ' \
            + fn_bam_1 + '.sample.bam -o ' + fn_bam_1 + '.sample.mpileup')
    print(command)
    call(command, shell=True)
    command = ('bcftools mpileup --count-orphans --annotate FORMAT/AD,FORMAT/DP -f ' + fn_ref + ' --min-BQ 0 --min-MQ 0 ' \
            + fn_bam_2 + '.sample.bam -o ' + fn_bam_2 + '.sample.mpileup')
    print(command)
    call(command, shell=True)
    
    f_mpileup_1 = pysam.VariantFile(fn_bam_1 + '.sample.mpileup')
    f_mpileup_2 = pysam.VariantFile(fn_bam_2 + '.sample.mpileup')
    
    print('[Biastools] Generate sample baseline')
    baseline(f_mpileup_1, f_mpileup_2, fn_sample, window_size)

    

