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
        print("Use --force option if you want to disable the prerequisite program check.")
        exit(1)


def bool2str(flag):
    if flag:
        return "1"
    else:
        return "0"


def catch_assert(parser, message):
    print('\n', message, '\n')
    parser.print_usage()
    exit(1)




def main():
    parser = argparse.ArgumentParser(description="Simulation/Alignment/Analyzing/Prediction module of the Biastools v0.1.1")
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
    parser.add_argument('-i', '--bam', help="Path to the alignment bam file, should be sorted [out_dir/sample.run_id.sorted.bam].")
    parser.add_argument('-n', '--naive', help= "Option to run the naive assignment method [False].", action='store_true')
    parser.add_argument('-R', '--real',  help= "Option for performing analysis on real data [False].", action='store_true')
    parser.add_argument('-d', '--boundary', help= "Boundary to plot the indel balance plot [20]", type=int, default=20)
    parser.add_argument('-lr', '--list_report', help= "List of bias report to plot the indel balance plot", nargs='+')
    parser.add_argument('-ld', '--list_run_id', help= "List of run ID for namings in the indel balance plot", nargs='+')
    # [4]
    parser.add_argument('-ps', '--sim_report',  help= "Path to the simulation report.")
    parser.add_argument('-pr', '--real_report', help= "Path to the real read report  [out_dir/sample.real.run.bias].")
    args = parser.parse_args()
    
    ##### Parameters for biastool_analysis
    path_output = args.out
    path_ref   = args.genome
    path_vcf   = args.vcf
    sample_id  = args.sample_id
    run_id     = args.run_id
    bam_file   = args.bam
    if bam_file == None:
        bam_file = path_output + '/' + sample_id + '.' + run_id + '.sorted.bam'
    
    flag_simulate = args.simulate
    flag_align    = args.align
    flag_analyze  = args.analyze
    flag_predict  = args.predict

    path_module = os.path.dirname(__file__) + '/'
    try:
        assert flag_simulate + flag_align + flag_analyze + flag_predict >= 1 
    except AssertionError:
        catch_assert(parser, "At least one of the --simulate/align/analyze/predict option should be specified.")

    flag_force = args.force
    thread = args.thread
    if thread == None:
        if sys.platform == "darwin":
            result = subprocess.run(["sysctl -n hw.ncpu"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        else:
            result = subprocess.run(["nproc"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        thread = int(result.stdout.strip())
    
    coverage = args.coverage
    aligner  = args.aligner
    align_index = args.align_index
    try:
        assert aligner=="bowtie2" or aligner=="bwamem" 
    except AssertionError:
        catch_assert(parser, "Only bowtie2 and bwamem are supported.")

    flag_naive  = args.naive
    flag_real   = args.real
    boundary    = args.boundary
    list_report = args.list_report
    list_run_id = args.list_run_id
    if list_report:
        try:
            assert len(list_report) == len(list_run_id)
        except AssertionError:
            catch_assert(parser, "Number of list --list_report and --list_run_id entries are inconsistent.")

    sim_report  = args.sim_report
    real_report = args.real_report
    if flag_predict: 
        try:
            assert real_report != None
        except AssertionError:
            catch_assert(parser, "<real_report> should be specified when using --predict")


    
    # Checking prerequisite programs are installed
    if flag_force != True:
        list_program = ["bedtools", \
                         "samtools", \
                         "bcftools", \
                         "gzip", \
                         "tabix"]
        if flag_align:
            list_program += ["bwa", "bowtie2"]
        if flag_simulate:
            list_program.append("mason_simulator")
        check_program_install( list_program ) 

    # Start running
    command = "mkdir -p " + path_output
    subprocess.call(command, shell=True)

    if flag_simulate:
        try:
            assert path_ref != None
            assert path_vcf != None 
        except AssertionError:
            catch_assert(parser, "<genome> and <vcf> should be specified when using --simulate")
        print("[Biastools] Simulate...")
        command = ' '.join(["bash", path_module+"biastools_simulation.sh", path_ref, path_vcf, path_output, sample_id, str(thread), str(coverage), path_module])
        #print(command)
        subprocess.call(command, shell=True)
    if flag_align:
        try:
            assert path_ref != None
            assert path_vcf != None
        except AssertionError:
            catch_assert(parser, "<genome> and <vcf> should be specified when using --align")
        if align_index == None:
            align_index = path_ref
        print("[Biastools] Align...")
        command = ' '.join(["bash", path_module+"biastools_align.sh", path_ref, path_vcf, path_output, sample_id, str(thread), aligner, align_index, run_id, path_module])
        #print(command)
        subprocess.call(command, shell=True)
    if flag_analyze:
        if list_report != None:
            print("[Biastools] Plot the indel balance plot for multiple bias reports...")
            if flag_real:
                subprocess.call(['python3', path_module+'indel_balance_plot.py', "-lr", ' '.join(list_report), "-ln", ' '.join(list_run_id), \
                                            "-vcf", path_output+"/"+sample_id+".het.vcf.gz", "-bd", str(boundary), "-map", \
                                            "-out", path_output+"/"+sample_id+"."+run_id+".real", "-real"])
            else:
                subprocess.call(['python3', path_module+'indel_balance_plot.py', "-lr"] + list_report + ["-ln"] + list_run_id + [ \
                                            "-vcf", path_output+"/"+sample_id+".het.vcf.gz", "-bd", str(boundary), "-map", \
                                            "-out", path_output+"/"+sample_id+"."+run_id+".sim"])
        else:
            try:
                assert path_ref != None
                assert path_vcf != None
            except AssertionError:
                catch_assert(parser, "<genome> and <vcf> should be specified when using --analyze")
            print("[Biastools] Analyze and plot...")
            command = ' '.join(["bash", path_module+"biastools_analysis.sh", path_ref, path_vcf, path_output, sample_id, str(thread), run_id, bool2str(flag_real), \
                                bool2str(flag_naive), str(boundary), path_module, bam_file])
            #print(command)
            subprocess.call(command, shell=True)
    if flag_predict:
        print("[Biastools] Predict bias...")
        command = ' '.join(["bash", path_module+"biastools_predict.sh", path_output, sample_id, run_id, bool2str(flag_real), real_report, sim_report, path_module])
        #print(command)
        subprocess.call(command, shell=True)




    
if __name__ == "__main__":
    main()



    
