# Wrap up python file for the biastools 3rd module
import subprocess
import sys
import os
import argparse
from biastools.biastools import check_program_install, catch_assert


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--out', help="Path to output directory ['out_dir'].", default="out_dir")
    parser.add_argument('-g', '--genome', help="Path to the reference genome.")
    parser.add_argument('-i', '--bam', help="Path to the alignment bam file, should be SORTED.")
    parser.add_argument('-s', '--sample_id', help="Sample ID ['sample'].", default="sample")
    parser.add_argument('-r', '--run_id', help="Run ID ['run'].", default="run")
    # Process options
    parser.add_argument('--scan',        help='[1] Option to scan and report bias region.', action='store_true')
    parser.add_argument('--compare_bam', help='[2] Option to generate common baseline and compare.', action='store_true')
    parser.add_argument('--compare_rpt', help='[3] Option to directly compare two bias report.', action='store_true')

    parser.add_argument('-t', '--thread', help="Number of threads to use [max].", type=int)
    parser.add_argument('--force', help="running the program without checking prerequisite programs.", action='store_true')
    # [1]
    parser.add_argument('-w', '--wig', help="Generate the wig files for the three measures, VERY SLOW [False]", action='store_true')
    parser.add_argument('-R', '--range', help="The range in the bam file targeted for analysis.")
    # [2]
    parser.add_argument('-i2', '--bam2',     help="Path to the second alignment bam file want to compare, should be SORTED.")
    parser.add_argument('-m',  '--mpileup',  help="Path to the mpileup file of the first bam file.")
    parser.add_argument('-m2', '--mpileup2', help="Path to the mpileup file of the second bam file.")
    # [3]
    parser.add_argument('-b1', '--bed1', help="Path to the first bed file for comparison.")
    parser.add_argument('-b2', '--bed2', help="Path to the second bed file for comparison.")
    parser.add_argument('-l2', '--lowRd2', help="Path to the .lowRd.bed report of the second file.")
    args = parser.parse_args()
    
    path_output = args.out
    path_ref   = args.genome
    bam_file   = args.bam
    sample_id  = args.sample_id
    run_id     = args.run_id
    
    flag_scan        = args.scan
    flag_compare_bam = args.compare_bam
    flag_compare_rpt = args.compare_rpt
    try:
        assert flag_scan + flag_compare_bam + flag_compare_rpt >= 1 
    except AssertionError:
        catch_assert(parser, "At least one of the --scan/compare_bam/compare_rpt option should be specified.")

    flag_force = args.force
    thread = args.thread
    if thread == None:
        if sys.platform == "darwin":
            result = subprocess.run(["sysctl -n hw.ncpu"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        else:
            result = subprocess.run(["nproc"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        thread = int(result.stdout.strip())
    flag_wig  = args.wig
    Range     = args.range 
    bam_file2 = args.bam2
    mpileup_file = args.mpileup
    mpileup_file2 = args.mpileup2
    bed_file1 = args.bed1
    bed_file2 = args.bed2
    lowRd_file2 = args.lowRd2

    
    # Checking prerequisite programs are installed
    if flag_force != True:
        check_program_install(["bedtools", \
                               "samtools", \
                               "bcftools"])

    # Start running
    command = "mkdir -p " + path_output
    subprocess.call(command, shell=True)
    prefix = path_output + '/' + sample_id 
    path_module = os.path.dirname(__file__) + '/'
    if flag_scan:
        print("[Biastools] Scanning...")
        if os.path.exists(bam_file+'.bai'):
            pass
        else:
            command = ["samtools", "index", bam_file]
            subprocess.call(command)
        
        print("[BIASTOOLS] SAMPLE", bam_file, " as ", sample_id + ".baseline ...")
        command = ["python3", path_module+"sample_baseline.py", "-b", bam_file, "-f", path_ref, "-o", prefix+".sample"]
        print(' '.join(command))
        subprocess.call(command)
        
        if Range == None:
            print("[BIASTOOLS] Process the whole bam file...")
            target_bam = bam_file
        else:
            print("[BIASTOOLS] Extract reads from " + Range + "...")
            target_bam = prefix + '.range.bam'
            command = ["samtools", "view", " -h", bam_file, Range, "-o", target_bam, "-@", thread]
            print(' '.join(command))
            subprocess.call(command)

        print("[BIASTOOLS] Format the mpileup...")
        if os.path.exists(prefix+'.'+run_id+'.mpileup'):
            print(prefix+'.'+run_id+'.mpileup already exist!')
        else:
            command = ["bcftools", "mpileup", "--count-orphans", "--annotate", "FORMAT/AD,FORMAT/DP", \
                       "-f", path_ref, \
                       "--min-BQ", "0", \
                       "--min-MQ", "0", \
                       "--threads", str(thread), target_bam, "-o", prefix+'.'+run_id+'.mpileup']
            print(' '.join(command))
            subprocess.call(command)
        print("[BIASTOOLS] Scanning bias...")
        if flag_wig:
            command = ["python3", path_module+"scanning_bias.py", "-g", prefix+'.'+run_id+'.mpileup', "--sample", "-b", prefix+".sample.baseline", \
                       "-wig", "-o", prefix+'.'+run_id+'.scanning']
        else:
            command = ["python3", path_module+"scanning_bias.py", "-g", prefix+'.'+run_id+'.mpileup', "--sample", "-b", prefix+".sample.baseline", \
                       "-o", prefix+'.'+run_id+'.scanning']
        print(' '.join(command))
        subprocess.call(command)
    
    if flag_compare_bam:
        if os.path.exists(bam_file+'.bai'):
            pass
        else:
            command = ["samtools", "index", bam_file]
            subprocess.call(command)
        if os.path.exists(bam_file2+'.bai'):
            pass
        else:
            command = ["samtools", "index", bam_file2]
            subprocess.call(command)

        print("[Biastools] Generate common baseline...")
        baseline = prefix+"."+run_id+".combine"
        command = ["python3", path_module+"merge_baseline.py", "-b1", bam_file, "-b2", bam_file2, "-f", path_ref, "-o", baseline]
        #print(' '.join(command))
        subprocess.call(command)
        command = ' '.join(["python3", path_module+"scanning_bias.py", "-g", mpileup_file,  "-b", baseline+".baseline", "-o", baseline+".1.scanning", ">", prefix+"."+run_id+".log"])
        #print(command)
        subprocess.call(command, shell=True)
        command = ' '.join(["python3", path_module+"scanning_bias.py", "-g", mpileup_file2, "-b", baseline+".baseline", "-o", baseline+".2.scanning", ">", prefix+"."+run_id+".log"])
        #print(command)
        subprocess.call(command, shell=True)

        print("[Biastools] Compare two bam files with common baseline...")
        command = ' '.join(["bash", path_module+"biastools_compare.sh", path_output, sample_id, run_id, \
                            baseline+".1.scanning.bias.bed", \
                            baseline+".2.scanning.bias.bed", \
                            baseline+".2.scanning.lowRd.bed", \
                            path_module])
        print(command)
        subprocess.call(command, shell=True)
    if flag_compare_rpt:
        print("[Biastools] Compare two bed files...")
        command = ' '.join(["bash", path_module+"biastools_compare.sh", path_output, sample_id, run_id, bed_file1, bed_file2, lowRd_file2, path_module])
        print(command)
        subprocess.call(command, shell=True)




if __name__ == "__main__":
    main()
