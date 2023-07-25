path_ref=$1
path_vcf=$2
path_out=$3
sample_id=$4
THR=$5
run_id=$6
flag_real=$7
flag_naive=$8
prefix=${path_out}/${sample_id}

echo "[Biastools] Reference bias analysis"

if [[ ${flag_real} == 1 ]]; then
    if [[ ${flag_naive} == 1 ]]; then
        python3 ref_bi_naive.py -s ${prefix}.${run_id}.sorted.het.bam \
                                -v ${prefix}.het.vcf.gz \
                                -f ${path_ref} \
                                -p ${prefix}.golden.rpt.pickle \
                                -o ${prefix}.real.${run_id}.bias \
                                --real
    else
        python3 ref_bi_context.py -s ${prefix}.${run_id}.sorted.het.bam \
                                  -v ${prefix}.het.vcf.gz \
                                  -f ${path_ref} \
                                  -p ${prefix}.golden.rpt.pickle \
                                  -o ${prefix}.real.${run_id}.bias \
                                  --real
    fi
else
    if [[ ${flag_naive} == 1 ]]; then
        python3 ref_bi_naive.py -s ${prefix}.${run_id}.sorted.het.bam \
                                -v ${prefix}.het.vcf.gz \
                                -f ${path_ref} \
                                -p ${prefix}.golden.rpt.pickle \
                                -o ${prefix}.sim.${run_id}.bias
    else
        python3 ref_bi_context.py -s ${prefix}.${run_id}.sorted.het.bam \
                                  -v ${prefix}.het.vcf.gz \
                                  -f ${path_ref} \
                                  -p ${prefix}.golden.rpt.pickle \
                                  -o ${prefix}.sim.${run_id}.bias
    fi
fi
