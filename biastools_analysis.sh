path_ref=$1
path_vcf=$2
path_out=$3
sample_id=$4
THR=$5
run_id=$6
flag_real=$7
flag_naive=$8
boundary=$9
prefix=${path_out}/${sample_id}

echo "[Biastools] Reference bias analysis"

if [[ ${flag_naive} == 1 ]]; then
    assign_method="ref_bi_naive.py"
else
    assign_method="ref_bi_context.py"
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
    python3 indel_balance_plot.py  -lr ${prefix}.${run_id}.real.bias \
                                   -ln ${run_id} \
                                   -vcf ${prefix}.het.vcf.gz \
                                   -bd ${boundary} \
                                   -map \
                                   -out ${r_prefix}.${run_id}.real
else
    python3 ${assign_method} -s ${prefix}.${run_id}.sorted.het.bam \
                             -v ${prefix}.het.vcf.gz \
                             -f ${path_ref} \
                             -p ${prefix}.golden.rpt.pickle \
                             -o ${prefix}.${run_id}.sim.bias
    
    # report the bias categories and report
    python3 golden_graph_report.py -mb ${prefix}.${run_id}.sim.bias.SNP -out ${r_prefix}.${run_id}.SNP
    python3 golden_graph_report.py -mb ${prefix}.${run_id}.sim.bias.gap -out ${r_prefix}.${run_id}.gap
    # plot the measures with NMB and NAB                      
    python3 golden_graph.py        -mb ${prefix}.${run_id}.sim.bias.SNP -out ${r_prefix}.${run_id}.SNP
    python3 golden_graph.py        -mb ${prefix}.${run_id}.sim.bias.gap -out ${r_prefix}.${run_id}.gap
    # indel balance plot
    python3 indel_balance_plot.py  -lr ${prefix}.${run_id}.sim.bias \
                                   -ln ${run_id} \
                                   -vcf ${prefix}.het.vcf.gz \
                                   -bd ${boundary} \
                                   -map \
                                   -out ${r_prefix}.${run_id}.sim
fi
