samtools faidx GRCh38_chr21.fa

bgzip -c renamed_chr21_NA12878.vcf > renamed_chr21_NA12878.vcf.gz
bcftools norm -f GRCh38_chr21.fa renamed_chr21_NA12878.vcf.gz -m +any -Oz -o normalized_chr21_NA12878.vcf.gz
bcftools index normalized_chr21_NA12878.vcf.gz

echo "[BIASTOOLS] Generate haplotype consensus reference sequence"
bcftools consensus -f GRCh38_chr21.fa -o GRCh38_chr21.hapA.fa -H 1 normalized_chr21_NA12878.vcf.gz
bcftools consensus -f GRCh38_chr21.fa -o GRCh38_chr21.hapB.fa -H 2 normalized_chr21_NA12878.vcf.gz
samtools faidx GRCh38_chr21.hapA.fa
samtools faidx GRCh38_chr21.hapB.fa

echo "[BIASTOOLS] Simulate sequences"
mason_simulator -ir GRCh38_chr21.hapA.fa -o hapA_1.fq -or hapA_2.fq -oa hapA.sam -n 3600000
mason_simulator -ir GRCh38_chr21.hapB.fa -o hapB_1.fq -or hapB_2.fq -oa hapB.sam -n 3600000

echo "[BIASTOOLS] Align sequences to the original reference"
bowtie2-build GRCh38_chr21.fa ./index/GRCh38_chr21/chr21_index
bowtie2 -p 32 -x ./index/GRCh38_chr21/chr21_index --rg-id hapA --rg SM:NA12878 -1 hapA_1.fq -2 hapA_2.fq -S hapA.bt2.sam
bowtie2 -p 32 -x ./index/GRCh38_chr21/chr21_index --rg-id hapB --rg SM:NA12878 -1 hapB_1.fq -2 hapB_2.fq -S hapB.bt2.sam

samtools sort -@ 16 hapA.bt2.sam -o hapA.bt2.sorted.bam
samtools sort -@ 16 hapB.bt2.sam -o hapB.bt2.sorted.bam

samtools merge -f bt2.sorted.bam hapA.bt2.sorted.bam hapB.bt2.sorted.bam

echo "[BIASTOOLS] Filter the heterozygous site in vcf file"
python3 filter_het_VCF.py -v normalized_chr21_NA12878.vcf.gz -o chr21_het.vcf.gz
tabix -p vcf chr21_het.vcf.gz

echo "[BIASTOOLS] Intersect the bam file and vcf file"
bedtools intersect -a bt2.sorted.bam -b chr21_het.vcf.gz | samtools view -bo bt2.sorted.het.bam
samtools index bt2.sorted.het.bam

echo "[BIASTOOLS] Test Octopus"
octopus -R GRCh38_chr21.fa -I bt2.sorted.het.bam -o chr21_octopus.vcf --bamout octopus.bt2.het.bam
python3 octopus_full_regrouper.py -p octopus.bt2.het.bam -s bt2.sorted.het.bam -o octopus.full.bt2.het.bam
samtools sort octopus.full.bt2.het.bam > octopus.bt2.sorted.het.bam

echo "[BIASTOOLS] Reference bias analysis"
#python3 ref_bi.py -s bt2.sorted.het.bam -v chr21_het.vcf.gz -f GRCh38_chr21.fa -o bt.bias
python3 ref_bi.py -s octopus.bt2.sorted.het.bam -v chr21_het.vcf.gz -f GRCh38_chr21.fa -o octopus.bt.bias
