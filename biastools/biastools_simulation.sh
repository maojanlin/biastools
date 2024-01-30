path_ref=$1
path_vcf=$2
path_out=$3
sample_id=$4
THR=$5
coverage=$6
path_module=$7
prefix=${path_out}/${sample_id}

if [ ! -f "${path_ref}.fai" ]; then
    samtools faidx ${path_ref}
fi

bcftools norm -f ${path_ref} ${path_vcf} -m +any -Oz -o ${prefix}.normalized.vcf.gz
bcftools index ${prefix}.normalized.vcf.gz

echo "[Biastools] Generate haplotype consensus reference sequence"
bcftools consensus -f ${path_ref} -o ${prefix}.hapA.fa -H 1 ${prefix}.normalized.vcf.gz -c ${prefix}.ref2hapA.chain
bcftools consensus -f ${path_ref} -o ${prefix}.hapB.fa -H 2 ${prefix}.normalized.vcf.gz -c ${prefix}.ref2hapB.chain
samtools faidx ${prefix}.hapA.fa
samtools faidx ${prefix}.hapB.fa

echo "[Biastools] Calculate how many reads should be generated"
total_base=$(( $( cut -f2 ${path_ref}.fai | paste -s -d+ ) ))
read_num=$(expr ${total_base} / 151 / 4 \* ${coverage})
echo "generating ${read_num} 2x151 reads in each haplotype"

echo "[Biastools] Simulate sequences"
mason_simulator --illumina-read-length 151 --num-threads ${THR} -ir ${prefix}.hapA.fa -o ${prefix}.hapA_1.fq -or ${prefix}.hapA_2.fq -oa ${prefix}.gt.hapA.sam -n ${read_num}
mason_simulator --illumina-read-length 151 --num-threads ${THR} -ir ${prefix}.hapB.fa -o ${prefix}.hapB_1.fq -or ${prefix}.hapB_2.fq -oa ${prefix}.gt.hapB.sam -n ${read_num} --seed 9388
samtools sort -@ ${THR} ${prefix}.gt.hapA.sam > ${prefix}.gt.hapA.sorted.bam
samtools sort -@ ${THR} ${prefix}.gt.hapB.sam > ${prefix}.gt.hapB.sorted.bam
samtools index ${prefix}.gt.hapA.sorted.bam
samtools index ${prefix}.gt.hapB.sorted.bam
rm ${prefix}.gt.hapA.sam
rm ${prefix}.gt.hapB.sam

gzip -f ${prefix}.hapA_1.fq
gzip -f ${prefix}.hapA_2.fq
gzip -f ${prefix}.hapB_1.fq
gzip -f ${prefix}.hapB_2.fq

echo "[Biastools] Filter the heterozygous site in vcf file"
python3 ${path_module}filter_het_VCF.py -v ${prefix}.normalized.vcf.gz  -o ${prefix}.het.vcf.gz
tabix -p vcf ${prefix}.het.vcf.gz

echo "[Biastools] Generate golden distribution report"
python3 ${path_module}consensus_vcf_map_adaptive.py -v ${prefix}.het.vcf.gz \
    -c0 ${prefix}.ref2hapA.chain \
    -c1 ${prefix}.ref2hapB.chain \
    -f0 ${prefix}.hapA.fa \
    -f1 ${prefix}.hapB.fa \
    -s0 ${prefix}.gt.hapA.sorted.bam \
    -s1 ${prefix}.gt.hapB.sorted.bam \
    -o  ${prefix}.golden.rpt


