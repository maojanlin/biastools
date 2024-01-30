# Tutorial: running biastools

In the tutorial, there are two initial files:
- ```grch38_chr20_part.fa```, which is the beginning 506,000 bases of the chr20 of GRCh38
- ```HG002.chr20.part.vcf.gz```, is the VCF file containing the first 1612 variants of HG002's chr20 called by Q100 project

After installation, the user can run the ```run.sh``` script, which simulate the reads from the reference genome and VCF file, and then align the simulated reads with Bowtie 2, 
and then analyze the alignment with context-aware assignment method. Finally, the biastools scan mode is used to scan the whole alignment bam file.

If the biastools is not installed, users can also directly called 
```
python3 biastools/biastools.py
```
and 
```
python3 biastools/biastools_scan.py
```
to run the procedure.
