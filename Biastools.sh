#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help() {
   # Display Help
   echo "Simulation/Alignment/Analyzing/Prediction module of the Biastools."
   echo
   echo "Syntax: Biastools.sh [-h] [options]  -1/2/3/4 -o <work_dir> -g <ref.fa> -v <vcf> -s 'sample_name' -r 'run0'"
   echo " Inputs:"
   echo "    -o path     Path to output directory ['out_dir']."
   echo "    -g path     Path to the reference genome."
   echo "    -v path     Path to the personal vcf file."
   echo "    -s string   Sample ID ['sample']."
   echo "    -r string   Run ID ['run']."
   echo " Process:"
   echo "    -1          Performing Simulation."
   echo "    -2          Performing Alignment."
   echo "    -3          Performing Analysis and Plotting."
   echo "    -4          Performing Prediction."
   echo " Options:"
   echo "    -h          Print Help information."
   echo "    -t INT      Number of threads to use [max]"
   echo "  @1[simulation]:"
   echo "    -x int      Read Coverage to simulate [30]."
   echo "  @2[aligment]"
   echo "    -a string   Aligner to use (bowtie2|bwamem) [bowtie2]"
   echo "    -b path     Path to the aligner index (target reference)"
   echo "  @3[analyzing and plot]:"
   echo "    -n          Option to run the naive assignment method [False]."
   echo "    -R          Option for performing analysis on real data [False]."
   echo "    -d INT      Boundary to plot the indel balance plot [20]"
   echo "    -l list     List of bias report to plot the indel balance plot"
   echo "    -L list     List of run ID for namings in the indel balance plot"
   echo "  @4[predict bias]:"
   echo "    -R          Option for performing analysis on real data [False]."
   echo "    -p path     Path to the simulation report."
   echo "    -P path     Path to the real read report  [out_dir/sample.real.run.bias]."
   echo
   exit 1
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Default Parameters
coverage=30
sample_id='sample'
path_output='out_dir'
ALN='bowtie2'
ALN_IDX='none'
report_real='none'
report_simulation='none'
run_id='run'
THR=$(nproc)
flag_real=0
flag_naive=0
list_report='none'
list_run_id='none'
boundary=20
while getopts "ho:g:v:s:x:r:a:b:t:p:P:l:L:d:Rn1234" option; do
   case $option in
      h) Help;;
      o) path_output=$OPTARG;;
      g) path_ref=$OPTARG;;
      v) path_vcf=$OPTARG;;
      s) sample_id=$OPTARG;;
      x) coverage=$OPTARG;;
      r) run_id=$OPTARG;;
      a) ALN=${OPTARG};;
      b) ALN_IDX=${OPTARG};;
      t) THR=${OPTARG};;
      p) report_simulation=${OPTARG};;
      P) report_real=${OPTARG};;
      l) list_report=${OPTARG};;
      L) list_run_id=${OPTARG};;
      d) boundary=${OPTARG};;
      R) flag_real=1;;
      n) flag_naive=1;;
      1) flag_simulation=1;;
      2) flag_alignment=1;;
      3) flag_analysis=1;;
      4) flag_predict=1;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done


if [[ ! ${ALN} =~ ^(bowtie2|bwamem)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem"
    exit 1
fi

mkdir -p ${path_output}
if [[ ${flag_simulation} == 1 ]]; then
    echo "[Biastools] Simulate..."
    bash biastools_simulation.sh ${path_ref} ${path_vcf} ${path_output} ${sample_id} ${THR} ${coverage}
fi
if [[ ${flag_alignment} == 1 ]]; then
    echo "[Biastools] Align..."
    bash biastools_align.sh ${path_ref} ${path_vcf} ${path_output} ${sample_id} ${THR} ${ALN} ${ALN_IDX} ${run_id}
fi
if [[ ${flag_analysis} == 1 ]]; then
    if [[ ${list_report} != 'none'  && ${list_run_id} != 'none' ]]; then
        echo "[Biastools] Plot the indel balance plot for multiple bias reports..."
        if [[ ${flag_real} == 1 ]]; then
            python3 indel_balance_plot.py  -lr ${list_report} \
                                           -ln ${list_run_id} \
                                           -vcf ${path_output}/${sample_id}.het.vcf.gz \
                                           -bd ${boundary} \
                                           -map \
                                           -out ${path_output}/${sample_id}.${run_id}.real \
                                           -real
        else
            python3 indel_balance_plot.py  -lr ${list_report} \
                                           -ln ${list_run_id} \
                                           -vcf ${path_output}/${sample_id}.het.vcf.gz \
                                           -bd ${boundary} \
                                           -map \
                                           -out ${path_output}/${sample_id}.${run_id}.sim
        fi
    else
        echo "[Biastools] Analyze and plot..."
        bash biastools_analysis.sh ${path_ref} ${path_vcf} ${path_output} ${sample_id} ${THR} ${run_id} ${flag_real} ${flag_naive} ${boundary}
    fi
fi
if [[ ${flag_predict} == 1 ]]; then
    echo "[Biastools] Predict bias..."
    bash biastools_predict.sh ${path_output} ${sample_id} ${run_id} ${flag_real} ${threshold} ${report_real} ${report_simulation}
fi

echo "Thank you for using Biastools!"
