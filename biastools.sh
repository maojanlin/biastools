bgzip -c renamed_chr21_NA12878.vcf > renamed_chr21_NA12878.vcf.gz

bcftools index renamed_chr21_NA12878.vcf.gz
bcftools consensus -f GRCh38_chr21.fa -o GRCh38_chr21.hapA.fa -H 1 renamed_chr21_NA12878.vcf.gz
bcftools consensus -f GRCh38_chr21.fa -o GRCh38_chr21.hapB.fa -H 2 renamed_chr21_NA12878.vcf.gz
samtools faidx GRCh38_chr21.hapA.fa
samtools faidx GRCh38_chr21.hapB.fa

mason_simulator -ir GRCh38_chr21.hapA.fa -o hapA_1.fq -or hapA_2.fq -oa hapA.sam -n 3600000
mason_simulator -ir GRCh38_chr21.hapB.fa -o hapB_1.fq -or hapB_2.fq -oa hapB.sam -n 3600000

bowtie2-build GRCh38_chr21.fa ./index/GRCh38_chr21/chr21_index
bowtie2 -p 32 -x ./index/GRCh38_chr21/chr21_index -1 hapA_1.fq -2 hapA_2.fq -S hapA.bt2.sam
bowtie2 -p 32 -x ./index/GRCh38_chr21/chr21_index -1 hapB_1.fq -2 hapB_2.fq -S hapB.bt2.sam

samtools sort -@ 16 -o hapA.bt2.sorted.sam hapA.bt2.sam
samtools sort -@ 16 -o hapB.bt2.sorted.sam hapB.bt2.sam

samtools merge -r bt2.sorted.sam hapA.bt2.sorted.sam hapB.bt2.sorted.sam
samtools view -ho bt2.sorted.bam bt2.sorted.sam

bedtools intersect -a bt2.sorted.bam -b chr21_het.vcf | samtools view -h > bt2.sorted.het.sam

python3 ref_bi.py -s bt2.sorted.het.sam -v chr21_het.vcf -f GRCh38_chr21.fa -o bt.bias
