
_Updated: Apr 17, 2025_
# Biastools: Measuring, visualizing and diagnosing reference bias

This github is originally forked from https://github.com/sheila12345/biastools

## Prerequisite programs
- samtools=v1.11
- bcftools=v1.9
- bedtools=v2.30.0
- gzip=v1.9
- tabix=v1.9
- bowtie2=v2.4.2
- bwa=v0.7.17
- mason_simulator=v2.0.9 (only for biastools --simulate)
- SeqAn=v2.4.0 (only for biastools --simulate)
  

## Installation
- [pip](https://pypi.org/project/biastools/)
```
pip install biastools
```
- [Github](https://github.com/maojanlin/biastools.git)
```
git clone https://github.com/maojanlin/biastools.git
cd biastools
```
Though optional, it is a good practice to install a virtual environment to manage the dependancies:

```
python -m venv venv
source venv/bin/activate
```
Now a virtual environment (named venv) is activated. Install biastools:

```
python setup.py install
```


## Usage

### Simulation, plotting, and analysis
```
$ biastools --simulate --align --analyze -o <work_dir> -g <ref.fa> -v <vcf> -s <sample_name> -r <run_id>
```

With the example command, biastools 
1. Simulates reads based on `<ref.fa>` and `<vcf>`, generating pair-end `.fq.gz` files for both haplotypes (`work_dir/sample_name.hap{A,B}_{1,2}.fq.gz`). 
2. Aligns the reads to the reference `<ref.fa>`, generating a BAM file with phasing information (`work_dir/sample_name.run_id.sorted.bam`).
3. Analyzes the BAM file with the context-aware assignment method, generating bias reports and plots.

#### Other aligners
Biastools supports [Bowtie 2](https://github.com/BenLangmead/bowtie2) and [bwa mem](https://github.com/lh3/bwa) aligners. BAM files from other aligners (named with `<work_dir/sample_name.run_id.sorted.bam>` and tagged with haplotype information) can be analyzed with

```
$ biastools --analyze -o <work_dir> -g <ref.fa> -v <vcf> -s <sample_name> -r <run_id>
```

#### Direct Analysis on Real sequence data
Biastools can also analyze real sequence data with the `--real` option using the context-aware assignment algorithm. The resulting plot does not include simulation information (`sample_id.real.indel_balance.pdf`).
```
$ biastools --analyze --real -t <thread> -o <work_dir> -g <ref.fa> -v <vcf> -s <sample_name> -r <run_id> \
                      --bam <path_to_target.bam>
```
Biastools first fetches the relevant alignments from the target BAM file, focusing only on heterozygous variant sites specified in the VCF file. These sites are then analyzed using a [context-aware algorithm](figures/context_aware.md). Finally, Biastools generates a bias report along with a bias-by-allele-length plot, both included in the output folder.


#### Combined Bias-by-allele-length plot
Multiple analysis results can be combined into a single Bias-by-allele-length plot. In biastools version 0.3.1, the default plotting module displays the 25th percentile, mean, and 75th percentile of the fraction of ALT alleles for variants stratified by allele length, using ticks to indicate the interquartile range and a central dot to mark the mean.

```
$ biastools --analyze -o <work_dir> -g <ref.fa> -v <vcf> -s <sample_name> -r <run_id> \
                      -lr file1.bias.all file2.bias.all file3.bias.all... \
                      -ld run_id1 run_id2 run_id3...
```

The output file `sample_name.combine.sim.indel_balance.pdf` plots the fraction of ALT alleles merged from the bias reports specified after the `-lr` option.  Users can use `-ld` option to specify the tool names, which will appear in the legend. To generate a combined plot using only real data bias reports (excluding simulation information), use the `--real` option.

An example of a combined bias-by-allele-length plot:
![multiple_indel_plot](figures/HG002.GIAB.4.2.1.demo.indel_balance.png?raw=true "multiple_indel_plot")


### Bias prediction from bias report
#### Real data
Biastools can predict if a variant is bias or not by:

```
$ biastools --predict -o <work_dir> -g <ref.fa> -v <vcf> -s <sample_name> -r <run_pd_id> -pr <path_to_bias_report>
```

With the example command, biastools
4. Generates two files: `sample_name.real.pd_id_bias.tsv` and `sample_name.real.pd_id_suspicious.tsv`. The `bias.tsv` report contains all sites predicted to be biased by the model. The `suspicious.tsv` file contains the sites which suspicious of lacking enough information from the VCF file. In another word, the reads align to the site shows different pattern to the haplotype indicated by the VCF file. 

#### Simulated guided prediction

```
$ biastools --predict -o <work_dir> -g <ref.fa> -v <vcf> -s <sample_name> -r <run_pd_id> \
                      -pr <path_to_bias_report> \
                      -ps <path_to_simulated_bias_report>
```

If the report of the sample based on simulated data is presented, biastools can generate cross prediction experiment result. In the experiment, the ground truth bias sites are based on simulation data.

### Scanning bias without vcf information
#### Scanning
```
$ biastools_scan --scan -o <work_dir> -g <ref.fa> -s <sample_name> -r <run_id> -i <path_to_target.bam>
```

Biastools transforms the `<path_to_target.bam>` into the mpileup format and generates baised and suspicious regions (`sample_name.run_id.bias.bed` and `sample_name.run_id.suspicious.bed`).


#### Compare two bam files with common baseline
```
$ biastools_scan --compare_bam -o <work_dir> -g <ref.fa> -s <sample_name> -r <run_id> \
                               -i  <path_to_target.bam> \
                               -i2 <path_to_second.bam> \
                               -m  <path_to_target.mpileup> \
                               -m2 <path_to_second.mpileup>
```
Biastools generates a common baseline from `path_to_target.bam` and `path_to_second.bam`, and uses the new common baseline to recalculate the bias regions based on the two mpileup files. The mpileup files can be generated by running **scanning** first, or directly run the **bcftools consensus**.



#### Directly compare two bias reports
User can also generate the comparison of the bias reports without a common baseline (not recommended):
```
$ biastools_scan --compare_rpt -o <work_dir> -s <sample_name> -r <run_id> \
                               -b1 <path_to_target_bias.bed> \
                               -b2 <path_to_improved_bias.bed> \
                               -l2 <path_to_improved_lowRd.bed>
```




