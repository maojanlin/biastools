path_ref=$1
path_vcf=$2
path_out=$3
sample_id=$4
THR=$5
ALN=$6
ALN_IDX=$7
run_id=$8
prefix=${path_out}/${sample_id}

echo "[Biastools] Align sequences to the original reference"
if [[ ${ALN_IDX} == 'none' ]]; then
    ALN_IDX=${path_ref}
fi

if [[ ${ALN} == "bowtie2" ]]; then
    if [ ! -f ${ALN_IDX}.1.bt2 ]; then
        bowtie2-build ${path_ref} ${ALN_IDX}
    fi
    bowtie2 -p ${THR} -x ${ALN_IDX} --rg-id ${run_id}_hapA --rg SM:${sample_id} -1 ${prefix}.hapA_1.fq.gz -2 ${prefix}.hapA_2.fq.gz |\
        samtools sort -o ${prefix}.hapA.${run_id}.sorted.bam  
    bowtie2 -p ${THR} -x ${ALN_IDX} --rg-id ${run_id}_hapB --rg SM:${sample_id} -1 ${prefix}.hapB_1.fq.gz -2 ${prefix}.hapB_2.fq.gz |\
        samtools sort -o ${prefix}.hapB.${run_id}.sorted.bam 
    samtools merge -f ${prefix}.${run_id}.sorted.bam ${prefix}.hapA.${run_id}.sorted.bam ${prefix}.hapB.${run_id}.sorted.bam
elif [[ ${ALN} == "bwamem" ]]; then
    if [ ! -f ${ALN_IDX}.bwt ]; then
        bwa index ${path_ref} -p ${ALN_IDX}
    fi
    bwa mem -t ${THR} ${ALN_IDX} ${run_id}_hapA_1.fq.gz ${run_id}_hapA_2.fq.gz -R "@RG\tID:${run_id}_hapA\tSM:${sample_id}" |\
        samtools sort -@ ${THR} -o ${prefix}.hapA.${run_id}.sorted.bam -
    bwa mem -t ${THR} ${ALN_IDX} ${run_id}_hapB_1.fq.gz ${run_id}_hapB_2.fq.gz -R "@RG\tID:${run_id}_hapB\tSM:${sample_id}" |\
        samtools sort -@ ${THR} -o ${prefix}.hapB.${run_id}.sorted.bam -
    samtools merge -f ${prefix}.${run_id}.sorted.bam ${prefix}.hapA.${run_id}.sorted.bam ${prefix}.hapB.${run_id}.sorted.bam
fi

echo "[Biastools] Intersect the bam file and vcf file"
bedtools intersect -a ${prefix}.${run_id}.sorted.bam -b ${prefix}.het.vcf.gz | samtools view -bo ${prefix}.${run_id}.sorted.het.bam
samtools index ${prefix}.${run_id}.sorted.het.bam

