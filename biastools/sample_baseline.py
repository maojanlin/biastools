import pysam
import argparse
import random
import numpy as np
import os
from subprocess import call

from scanning_bias import scanning_bias, calculate_measures 


def baseline(
        f_mpileup   :pysam.VariantRecord,
        fn_sample   :str,
        window_size :int
        ) -> tuple:
    """
    Take in the sample mpileup, and output the average read_depth/variant Density/non Diploid portion
    """
    dict_ref_info = scanning_bias(f_gvcf=f_mpileup)
    dict_3D_measures = calculate_measures(
        dict_ref_info=dict_ref_info,
        window_size=window_size
        )
    
    total_read_depth  = np.array([])
    total_var_density = np.array([])
    total_dip_density = np.array([])
    
    fo = open(fn_sample + '.baseline', 'w')
    fo.write('#chr pos segment_len RD_mean RD_std VD_mean VD_std ND_mean ND_std\n')
    for ref_name, dict_array in dict_3D_measures.items():
        for start_pos, array_info in dict_array.items():
            array_read_depth, array_var_density, array_dip_density = array_info
            
            avg_read_depth   = np.mean(array_read_depth)
            std_read_depth   = np.std(array_read_depth)
            #positive_avg_var = np.mean(array_var_density)
            #positive_std_var = np.std(array_var_density)
            #positive_avg_dip = np.mean(array_dip_density)
            #positive_std_dip = np.std(array_dip_density)

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
            #total_var_density = np.concatenate((total_var_density, array_var_density))
            #total_dip_density = np.concatenate((total_dip_density, array_dip_density))
    
    fo.write('#total sample len: ' + str(len(total_read_depth)) + '\n')
    fo.write('#total_statistics:\n')
    fo.write('#chr pos segment_len RD_mean RD_std VD_mean VD_std ND_mean ND_std\n# ')
    fo.write(str(round(np.mean(total_read_depth),5))  + ' ' + str(round(np.std(total_read_depth),5)) + ' ')
    fo.write(str(round(np.mean(total_var_density),5)) + ' ' + str(round(np.std(total_var_density),5)) + ' ')
    fo.write(str(round(np.mean(total_dip_density),5)) + ' ' + str(round(np.std(total_dip_density),5)))
    fo.close()
    return np.mean(total_read_depth), np.mean(total_var_density), np.mean(total_dip_density)


def sample_select(
        fn_sample   :str,
        seed        :int,
        min_len     :int,
        f_bam       :pysam.AlignmentFile
    ):
    """
    Take out the contig length greater than min_len (threshold_contig)
    For each contig, takes 100 segments totally equal to 1/1000 of the contig length
    """
    random.seed(seed)
    fo = open(fn_sample + '.bed', 'w')
    write_flag = False
    for idx, name in enumerate(f_bam.header.references):
        contig_len = f_bam.header.lengths[idx]
        if contig_len > min_len:
            write_flag = True
            thousandth = int(contig_len / 100000)
            list_sample_start = random.sample(range(100000 - 1), 100)
            for sample_start in sorted(list_sample_start):
                fo.write(name + ' ' + str(sample_start*thousandth) + ' ' + str(sample_start*thousandth+thousandth) + '\n')
        elif write_flag == False:
            fo.write(name + ' 1 ' + str(contig_len) + '\n')
            write_flag = True
    fo.close()





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam_file', help='the bam file we want to sample')
    parser.add_argument('-f', '--reference_fasta', help='the reference fasta file for mpileup building')
    parser.add_argument('-o', '--sample_bed', help='the sampled 1/1000 bed file')
    parser.add_argument('-w', '--window_size', help='window size for average depth, density analysis', type=int, default=400)
    parser.add_argument('-th', '--threshold_contig', help='the minimum contig length for sampling', type=int, default=10000000)
    parser.add_argument('--seed', help='seed for random sampling', type=int, default=0)
    parser.add_argument('-k', '--kill', help='kill all storage files', action='store_true')
    args = parser.parse_args()
    
    fn_bam    = args.bam_file
    fn_ref    = args.reference_fasta
    fn_sample = args.sample_bed
    min_len   = args.threshold_contig
    window_size = args.window_size
    seed      = args.seed
    kill_flag = args.kill

    f_bam = pysam.AlignmentFile(fn_bam)

    # sample bed file according to the bam file information
    sample_select(fn_sample, seed, min_len, f_bam)


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
    baseline(f_mpileup, fn_sample, window_size)

    
