#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help() {
   # Display Help
   echo "Scanning module of the Biastools."
   echo
   echo "Syntax: Biastools_scan.sh [-h] [options]  -1 -o <work_dir> -g <ref.fa> -s 'sample_name' -r 'run_id'"
   echo " Inputs:"
   echo "    -o path     Path to output directory ['out_dir']."
   echo "    -g path     Path to the reference genome."
   echo "    -i path     Path to the alignment bam file, should be SORTED."
   echo "    -s string   Sample ID ['sample']."
   echo "    -r string   Run ID ['run']."
   echo " Process:"
   echo "    -1          Performing Scanning."
   echo "    -2          Performing Common Basline Comparison"
   echo "    -3          Performing Comparison."
   echo " Options:"
   echo "    -h          Print Help information."
   echo "    -t INT      Number of threads to use [max]."
   echo "  @1[scanning]:"
   echo "    -w          Generate the wig files for the three measures, VERY SLOW [False]."
   echo "  @2[baseline and compare]"
   echo "    -i path     Path to the alignment bam file."
   echo "    -I path     Path to the alignment bam file."
   echo "    -m path     Path to the alignment mpileup file."
   echo "    -M path     Path to the alignment mpileup file."
   echo "  @3[compare]"
   echo "    -b path     Path to 'target' bed file."
   echo "    -B path     Path to 'improved' bed file."
   echo "    -l path     Path to the .lowRd.bed report of the 'improved' file."
   echo
   exit 1
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Default Parameters
sample_id='sample'
path_output='out_dir'
run_id='run'
THR=$(nproc)
Range=""
flag_wig=0
while getopts "ho:g:i:I:s:r:t:Rwb:B:l:m:M:123" option; do
   case $option in
      h) Help;;
      o) path_output=$OPTARG;;
      g) path_ref=$OPTARG;;
      i) bam_file=$OPTARG;;
      I) bam_file_2=$OPTARG;;
      s) sample_id=$OPTARG;;
      r) run_id=$OPTARG;;
      t) THR=${OPTARG};;
      R) Range=$OPTARG;;
      w) flag_wig=1;;
      b) target_bed=${OPTARG};;
      B) improve_bed=${OPTARG};;
      l) improve_lowRd=${OPTARG};;
      m) mpileup_1=${OPTARG};;
      M) mpileup_2=${OPTARG};;
      1) flag_scanning=1;;
      2) flag_compare_baseline=1;;
      3) flag_compare=1;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done


mkdir -p ${path_output}
prefix=${path_output}/${sample_id}
if [[ ${flag_scanning} == 1 ]]; then
    echo "[Biastools] Scanning..."
    if [[ ! -f ${bam_file}.bai ]]; then
        samtools index ${bam_file}
    fi
    
    echo "[BIASTOOLS] SAMPLE THE WHOLE ${bam_file} as ${sample_id}.baseline ..." 
    #TODO sample_baseline can be used together with the $Range, else total sample can be used
    python3 sample_baseline.py -b ${bam_file} -f ${path_ref} -o ${prefix}.sample
    
    if [[ ${Range} == "" ]]; then
        echo "[BIASTOOLS] Process the whole bam file..."
        target_bam=${bam_file}
    else
        echo "[BIASTOOLS] Extract reads from ${Range}..."
        target_bam=${prefix}.range.bam
        samtools view -h ${bam_file} ${Range} -o ${target_bam} -@ ${THR}
    fi
    
    echo "[BIASTOOLS] Format the mpileup..."
    if [[ ! -f ${prefix}.${run_id}.mpileup ]]; then
        bcftools mpileup \
          --count-orphans \
          --annotate FORMAT/AD,FORMAT/DP\
          -f "${path_ref}" \
          --min-BQ 0 \
          --min-MQ 0 \
          --threads ${THR} \
          ${target_bam} \
          -o ${prefix}.${run_id}.mpileup
    else
        echo "${prefix}.${run_id}.mpileup already exist!"
    fi
    
    echo "[BIASTOOLS] Scanning bias..."
    if [[ ${flag_wig} == 1 ]]; then
        python3 scanning_bias.py -g ${prefix}.${run_id}.mpileup --sample -b ${prefix}.sample.baseline -wig -o ${prefix}.${run_id}.scanning
    else
        python3 scanning_bias.py -g ${prefix}.${run_id}.mpileup --sample -b ${prefix}.sample.baseline -o ${prefix}.${run_id}.scanning
    fi
fi

if [[ ${flag_compare_baseline} == 1 ]]; then
    if [[ ! -f ${bam_file}.bai ]]; then
        samtools index ${bam_file}
    fi
    if [[ ! -f ${bam_file_2}.bai ]]; then
        samtools index ${bam_file_2}
    fi

    echo "[Biastools] Generate common baseline..."
    python3 merge_baseline.py -b1 ${bam_file} -b2 ${bam_file_2} -f ${path_ref} -o ${prefix}.combine
    python3 scanning_bias.py -g ${mpileup_1} -b ${prefix}.combine.baseline -o ${prefix}.${run_id}.combine.1.scanning > ${prefix}.${run_id}.log
    python3 scanning_bias.py -g ${mpileup_2} -b ${prefix}.combine.baseline -o ${prefix}.${run_id}.combine.2.scanning > ${prefix}.${run_id}.log

    echo "[Biastools] Compare two bam files with common baseline..."
    bash biastools_compare.sh ${path_output} ${sample_id} ${run_id} \
                              ${prefix}.${run_id}.combine.1.scanning.bias.bed \
                              ${prefix}.${run_id}.combine.2.scanning.bias.bed \
                              ${prefix}.${run_id}.combine.2.scanning.lowRd.bed
fi

if [[ ${flag_compare} == 1 ]]; then
    echo "[Biastools] Compare two bed files..."
    bash biastools_compare.sh ${path_output} ${sample_id} ${run_id} ${target_bed} ${improve_bed} ${improve_lowRd}
fi

echo "Thank you for using Biastools!"
