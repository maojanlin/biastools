path_ref=$1
path_vcf=$2
path_out=$3
sample_id=$4
THR=$5
run_id=$6
flag_real=$7
flag_naive=$8
boundary=$9
path_module=${10}
prefix=${path_out}/${sample_id}
bam_file=${11}


echo "[Biastools] Intersect the bam file and vcf file"
if [ ! -f ${prefix}.het.vcf.gz ]; then
    bcftools norm -f ${path_ref} ${path_vcf} -m +any -Oz -o ${prefix}.normalized.vcf.gz
    bcftools index ${prefix}.normalized.vcf.gz
    python3 ${path_module}filter_het_VCF.py -v ${prefix}.normalized.vcf.gz  -o ${prefix}.het.vcf.gz
    tabix -p vcf ${prefix}.het.vcf.gz
fi
if [ ! -f ${prefix}.${run_id}.sorted.het.bam ]; then
    bedtools intersect -a ${bam_file} -b ${prefix}.het.vcf.gz | samtools sort -@ ${THR} > ${prefix}.${run_id}.sorted.het.bam
    samtools index ${prefix}.${run_id}.sorted.het.bam
fi


echo "[Biastools] Reference bias analysis"
if [[ ${flag_naive} == 1 ]]; then
    assign_method=${path_module}"ref_bi_naive.py"
else
    assign_method=${path_module}"ref_bi_context.py"
fi

mkdir -p ${path_out}/${run_id}"_report"
r_prefix=${path_out}/${run_id}"_report"/${sample_id}
if [[ ${flag_real} == 1 ]]; then
    python3 ${assign_method} -s ${prefix}.${run_id}.sorted.het.bam \
                             -v ${prefix}.het.vcf.gz \
                             -f ${path_ref} \
                             -p ${prefix}.golden.rpt.pickle \
                             -o ${prefix}.${run_id}.real.bias \
                             --real
    # indel balance plot
    python3 ${path_module}indel_balance_plot.py  -lr ${prefix}.${run_id}.real.bias \
                                             -ln ${run_id} \
                                             -vcf ${prefix}.het.vcf.gz \
                                             -bd ${boundary} \
                                             -map \
                                             -out ${r_prefix}.${run_id}.real \
                                             -real
else
    python3 ${assign_method} -s ${prefix}.${run_id}.sorted.het.bam \
                             -v ${prefix}.het.vcf.gz \
                             -f ${path_ref} \
                             -p ${prefix}.golden.rpt.pickle \
                             -o ${prefix}.${run_id}.sim.bias
    
    # report the bias categories and report
    python3 ${path_module}golden_graph_report.py -mb ${prefix}.${run_id}.sim.bias.SNP -out ${r_prefix}.${run_id}.SNP
    python3 ${path_module}golden_graph_report.py -mb ${prefix}.${run_id}.sim.bias.gap -out ${r_prefix}.${run_id}.gap
    # plot the measures with NMB and NAB                      
    python3 ${path_module}golden_graph.py        -mb ${prefix}.${run_id}.sim.bias.SNP -out ${r_prefix}.${run_id}.SNP
    python3 ${path_module}golden_graph.py        -mb ${prefix}.${run_id}.sim.bias.gap -out ${r_prefix}.${run_id}.gap
    # indel balance plot
    python3 ${path_module}indel_balance_plot.py  -lr ${prefix}.${run_id}.sim.bias \
                                                 -ln ${run_id} \
                                                 -vcf ${prefix}.het.vcf.gz \
                                                 -bd ${boundary} \
                                                 -map \
                                                 -out ${r_prefix}.${run_id}.sim
fi
