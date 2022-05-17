Work_dir="scanning_bias"
Range="chr17"

mkdir -p ${Work_dir}

Ref_fasta="/home/mlin77/data_blangme2/fasta/grch37/hg19.fa"
#BAM_file="/data/blangme2/naechyun/data_leviosam/illumina_wgs/hg002/bwa/HG002.novaseq.pcr-free.30x.bwa.grch37.bam"
BAM_file="/data/blangme2/naechyun/data_leviosam/illumina_wgs/hg002/experiments/chm13v2/chm13_hg19/HG002-bwa-chm13v2_hg19-final.bam"
sample_name="HG002.wgs.leviosam.sample"
#prefix="HG002.leviosam.chr17"
prefix="HG002_target.chm13_hg19.chr17.all"

echo "[BIASTOOLS] SAMPLE THE WHOLE ${BAM_file} as ${sample_name}.baseline ..." 
python3 sample_baseline.py -b ${BAM_file} -f ${Ref_fasta} -o ${sample_name}

echo "[BIASTOOLS] EXTRACT READS FOR ${Range}"
samtools view -h ${BAM_file} ${Range} -o ${prefix}.bam

echo "[BIASTOOLS] FORMAT VARIANT"
bcftools mpileup \
  --count-orphans \
  --annotate FORMAT/AD,FORMAT/DP\
  -f "${Ref_fasta}" \
  --min-BQ 0 \
  --min-MQ 0 \
  ${prefix}.bam \
  -o ${Work_dir}/${prefix}.mpileup

echo "[BIASTOOLS] SCANNING THE REGION"
python3 scanning_bias.py -g ${Work_dir}/${prefix}.mpileup -b ${sample_name}.baseline -wig -o ${prefix}.scanning

