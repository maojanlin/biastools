#!/bin/sh

bcftools view -e 'GT~"0|."' -G -o hapA.vcf renamed_chrv21_NA12878.vcf

bcftools view -e 'GT~".|0"' -G -o hapB.vcf renamed_chr21_NA12878.vcf

mason_simulator -iv hapA.vcf -ir GRCh38_chr21.fa -o hapA1.fq -or hapA2.fq -oa hapA.sam  -n 4800000

mason_simulator -iv hapB.vcf -ir GRCh38_chr21.fa -o hapB1.fq -or hapB2.fq -oa hapB.sam  -n 4800000

samtools sort hapA.sam -o sorted_hapA.sam

samtools sort hapB.sam -o sorted_hapB.sam

bcftools view -i 'GT="het"' -o chr21_het.vcf renamed_chr21_NA12878.vcf

python3 ref_bi.py -v chr21_het.vcf -s sorted_hapA.sam -f chr21_fasta.fa -o hapA_ref_bi.txt