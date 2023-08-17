# Wrap up python file for the biastools 1st and 2nd module
import subprocess
import sys
import os
import argparse
from shutil import which

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


def check_program_install(list_names):
    flag_violate = False
    for name in list_names:
        if is_tool(name) == False:
            print(name, "is a prerequisite program, please install it before running biastools")
            flag_violate = True
    if flag_violate:
        exit()


def bool2str(flag):
    if flag:
        return "1"
    else:
        return "0"




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--out', help="Path to output directory ['out_dir'].", default="out_dir")
    parser.add_argument('-g', '--genome', help="Path to the reference genome.")
    parser.add_argument('-v', '--vcf', help="Path to the personal vcf file.")
    parser.add_argument('-s', '--sample_id', help="Sample ID ['sample'].", default="sample")
    parser.add_argument('-r', '--run_id', help="Run ID ['run'].", default="run")
    # Process options
    parser.add_argument('--simulate', help='[1] Option to run biastools simulation.', action='store_true')
    parser.add_argument('--align',    help='[2] Option to run biastools align.', action='store_true')
    parser.add_argument('--analyze',  help='[3] Option to run biastools analyze.', action='store_true')
    parser.add_argument('--predict',  help='[4] Option to predict bias from analysis report.', action='store_true')

    parser.add_argument('-t', '--thread', help="Number of threads to use [max].", type=int)
    parser.add_argument('--force', help="running the program without checking prerequisite programs.", action='store_true')
    # [1]
    parser.add_argument('-x', '--coverage', help="Read coverage to simulate [30].", type=int, default=30)
    # [2]
    parser.add_argument('-a', '--aligner', help="Aligner to use (bowtie2|bwamem) [bowtie2]", default="bowtie2")
    parser.add_argument('-b', '--align_index', help="Path to the aligner index (target reference)")
    # [3]
    parser.add_argument('-n', '--naive', help= "Option to run the naive assignment method [False].", action='store_true')
    parser.add_argument('-R', '--real',  help= "Option for performing analysis on real data [False].", action='store_true')
    parser.add_argument('-d', '--boundary', help= "Boundary to plot the indel balance plot [20]", type=int, default=20)
    parser.add_argument('-lr', '--list_report', help= "List of bias report to plot the indel balance plot", nargs='+')
    parser.add_argument('-ld', '--list_run_id', help= "List of run ID for namings in the indel balance plot", nargs='+')
    # [4]
    parser.add_argument('-ps', '--sim_report',  help= "Path to the simulation report.")
    parser.add_argument('-pr', '--real_report', help= "Path to the real read report  [out_dir/sample.real.run.bias].")
    
    parser.add_argument('--scan',        help='[1] Option to scan and report bias region.', action='store_true')
    parser.add_argument('--compare_bam', help='[2] Option to generate common baseline and compare.', action='store_true')
    parser.add_argument('--compare_rpt', help='[3] Option to directly compare two bias report.', action='store_true')
    parser.add_argument('-i', '--bam', help="Path to the alignment bam file, should be SORTED.")
    # [1]
    parser.add_argument('-w',  '--wig', help="Generate the wig files for the three measures, VERY SLOW [False]", action='store_true')
    parser.add_argument('-rg', '--range', help="The range in the bam file targeted for analysis.")
    # [2]
    parser.add_argument('-i2', '--bam2',     help="Path to the second alignment bam file want to compare, should be SORTED.")
    parser.add_argument('-m',  '--mpileup',  help="Path to the mpileup file of the first bam file.")
    parser.add_argument('-m2', '--mpileup2', help="Path to the mpileup file of the second bam file.")
    # [3]
    parser.add_argument('-b1', '--bed1', help="Path to the first bed file for comparison.")
    parser.add_argument('-b2', '--bed2', help="Path to the second bed file for comparison.")
    parser.add_argument('-l2', '--lowRd2', help="Path to the .lowRd.bed report of the second file.")
    args = parser.parse_args()
    
    ##### Parameters for biastool_analysis
    path_output = args.out
    path_ref   = args.genome
    path_vcf   = args.vcf
    sample_id  = args.sample_id
    run_id     = args.run_id
    
    flag_simulate = args.simulate
    flag_align    = args.align
    flag_analyze  = args.analyze
    flag_predict  = args.predict
    

    flag_force = args.force
    thread = args.thread
    if thread == None:
        result = subprocess.run(["nproc"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        thread = int(result.stdout.strip())
    
    coverage = args.coverage
    aligner  = args.aligner
    assert aligner=="bowtie2" or aligner=="bwamem", "only bowtie2 and bwamem are supported."
    align_index = args.align_index

    flag_naive  = args.naive
    flag_real   = args.real
    boundary    = args.boundary
    list_report = args.list_report
    list_run_id = args.list_run_id
    if list_report:
        assert len(list_report) == len(list_run_id), " Number of list --list_report and --list_run_id entries are inconsistent."

    sim_report  = args.sim_report
    real_report = args.real_report
    if flag_predict: 
        assert real_report != None, "--real_report should be specified when using --predict"


    ##### Parameters for biastool_scan
    flag_scan        = args.scan
    flag_compare_bam = args.compare_bam
    flag_compare_rpt = args.compare_rpt
    assert flag_simulate + flag_align + flag_analyze + flag_predict >= 1 or \
           flag_scan + flag_compare_bam + flag_compare_rpt >= 1, """at least one of the --simulate/align/analyze/predict or
                                                                                        --scan/compare_bam/compare_rpt option should be specified."""
    assert ((flag_simulate + flag_align + flag_analyze + flag_predict >= 1) and \
           (flag_scan + flag_compare_bam + flag_compare_rpt >= 1)) == False, "the analyze options and scan options should not be specified at the same time!"
    
    flag_wig  = args.wig
    Range     = args.range 
    bam_file   = args.bam
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
                               "bcftools", \
                               "bwa", \
                               "bowtie2", \
                               "gzip", \
                               "tabix", \
                               "mason_simulator"])

    # Start running
    command = "mkdir -p " + path_output
    subprocess.call(command, shell=True)

    if flag_simulate:
        assert path_ref != None, "--genome should be specified when using --simulate"
        assert path_vcf != None, "--vcf should be specified when using --simulate"
        print("[Biastools] Simulate...")
        command = ' '.join(["bash biastools_simulation.sh", path_ref, path_vcf, path_output, sample_id, str(thread), str(coverage)])
        #print(command)
        subprocess.call(command, shell=True)
    if flag_align:
        assert path_ref != None, "--genome should be specified when using --align"
        assert path_vcf != None, "--vcf should be specified when using --align"
        if align_index == None:
            align_index = path_ref
        print("[Biastools] Align...")
        command = ' '.join(["bash biastools_align.sh", path_ref, path_vcf, path_output, sample_id, str(thread), aligner, align_index, run_id])
        #print(command)
        subprocess.call(command, shell=True)
    if flag_analyze:
        if list_report != None:
            print("[Biastools] Plot the indel balance plot for multiple bias reports...")
            if flag_real:
                subprocess.call(['python3', 'indel_balance_plot.py', "-lr", ' '.join(list_report), "-ln", ' '.join(list_run_id), \
                                            "-vcf", path_output+"/"+sample_id+".het.vcf.gz", "-bd", str(boundary), "-map", \
                                            "-out", path_output+"/"+sample_id+"."+run_id+".real", "-real"])
            else:
                subprocess.call(['python3', 'indel_balance_plot.py', "-lr"] + list_report + ["-ln"] + list_run_id + [ \
                                            "-vcf", path_output+"/"+sample_id+".het.vcf.gz", "-bd", str(boundary), "-map", \
                                            "-out", path_output+"/"+sample_id+"."+run_id+".sim"])
        else:
            assert path_ref != None, "--genome should be specified when using --analyze"
            assert path_vcf != None, "--vcf should be specified when using --analyze"
            print("[Biastools] Analyze and plot...")
            command = ' '.join(["bash biastools_analysis.sh", path_ref, path_vcf, path_output, sample_id, str(thread), run_id, bool2str(flag_real), \
                                bool2str(flag_naive), str(boundary)])
            #print(command)
            subprocess.call(command, shell=True)
    if flag_predict:
        print("[Biastools] Predict bias...")
        command = ' '.join(["bash biastools_predict.sh", path_output, sample_id, run_id, bool2str(flag_real), real_report, sim_report])
        #print(command)
        subprocess.call(command, shell=True)



    ######## Running biastools_scan
    prefix = path_output + '/' + sample_id 
    if flag_scan:
        print("[Biastools] Scanning...")
        if os.path.exists(bam_file+'.bai'):
            pass
        else:
            command = ["samtools", "index", bam_file]
            subprocess.call(command)
        
        print("[BIASTOOLS] SAMPLE", bam_file, " as ", sample_id + ".baseline ...")
        command = ["python3", "sample_baseline.py", "-b", bam_file, "-f", path_ref, "-o", prefix+".sample"]
        #print(' '.join(command))
        subprocess.call(command)
        
        if Range == None:
            print("[BIASTOOLS] Process the whole bam file...")
            target_bam = bam_file
        else:
            print("[BIASTOOLS] Extract reads from " + Range + "...")
            target_bam = prefix + '.range.bam'
            command = ["samtools", "view", " -h", bam_file, Range, "-o", target_bam, "-@", thread]
            #print(' '.join(command))
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
            #print(' '.join(command))
            subprocess.call(command)
        print("[BIASTOOLS] Scanning bias...")
        if flag_wig:
            command = ["python3", "scanning_bias.py", "-g", prefix+'.'+run_id+'.mpileup', "--sample", "-b", prefix+".sample.baseline", \
                       "-wig", "-o", prefix+'.'+run_id+'.scanning']
        else:
            command = ["python3", "scanning_bias.py", "-g", prefix+'.'+run_id+'.mpileup', "--sample", "-b", prefix+".sample.baseline", \
                       "-o", prefix+'.'+run_id+'.scanning']
        #print(' '.join(command))
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
        command = ["python3", "merge_baseline.py", "-b1", bam_file, "-b2", bam_file2, "-f", path_ref, "-o", baseline]
        #print(' '.join(command))
        subprocess.call(command)
        command = ' '.join(["python3", "scanning_bias.py", "-g", mpileup_file,  "-b", baseline+".baseline", "-o", baseline+".1.scanning", ">", prefix+"."+run_id+".log"])
        #print(command)
        subprocess.call(command, shell=True)
        command = ' '.join(["python3", "scanning_bias.py", "-g", mpileup_file2, "-b", baseline+".baseline", "-o", baseline+".2.scanning", ">", prefix+"."+run_id+".log"])
        #print(command)
        subprocess.call(command, shell=True)

        print("[Biastools] Compare two bam files with common baseline...")
        command = ' '.join(["bash biastools_compare.sh", path_output, sample_id, run_id, \
                            baseline+".1.scanning.bias.bed", \
                            baseline+".2.scanning.bias.bed", \
                            baseline+".2.scanning.lowRd.bed"])
        #print(command)
        subprocess.call(command, shell=True)
    if flag_compare_rpt:
        print("[Biastools] Compare two bed files...")
        command = ' '.join(["bash biastools_compare.sh", path_output, sample_id, run_id, bed_file1, bed_file2, lowRd_file2])
        #print(command)
        subprocess.call(command, shell=True)


    



    
