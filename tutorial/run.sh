biastools --simulate -g grch38_chr20_part.fa -v HG002.chr20.part.vcf.gz -s HG002_part -r tutorial
biastools --align -a bowtie2 -g grch38_chr20_part.fa -v HG002.chr20.part.vcf.gz -s HG002_part -r tutorial
biastools --analyze -g grch38_chr20_part.fa -v HG002.chr20.part.vcf.gz -s HG002_part -r tutorial
biastools_scan --scan -g grch38_chr20_part.fa -s HG002_part -r tutorial -i out_dir/HG002_part.tutorial.sorted.bam
