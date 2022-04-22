Work_dir="scanning_bias"
Range="chr7:151832010-152133088"

mkdir ${Work_dir}

Ref_fasta="/home/mlin77/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
Bam_file_1="/scratch16/blangme2/naechyun/leviosam_exp/hg002_30x/hg19-kmt2c.bam"
Bam_file_2="/data/blangme2/naechyun/data_leviosam/illumina_wgs/hg002/experiments/chm13v2/chm13_hg19/chm13v2_hg19-kmt2c.bam"
prefix="KMT2C_chm13"
prefix="KMT2C_hg19"

echo "[BIASTOOLS] Format Variant"
bcftools mpileup \
  -A \
  --annotate FORMAT/AD,FORMAT/DP\
  -f "${Ref_fasta}" \
  --min-BQ 0 \
  --min-MQ 0 \
  ${Bam_file_1} \
  -o ${Work_dir}/${prefix}.mpileup

python3 scanning_bias.py -g ${Work_dir}/${prefix}.mpileup


#|\
#    bcftools call \
#      --multiallelic-caller \
#      --gvcf 1 \
#      --ploidy GRCh37 \
#      -Ob > ${Work_dir}/${prefix}.bcf

  #--redo-BAQ \
  #--per-sample-mF \
      #--variants-only \
  #--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \

