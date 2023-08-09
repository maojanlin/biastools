_Updated: Aug 6, 2023_
# Biastools: Measuring, visualizing and diagnosing reference bias

## Usage:

### Simulation, plotting, and analysis
```
$ Biastools.sh -123 -o <work_dir> -g <ref.fa> -v <vcf> -s 'sample_name' -r 'run_id'
```

With the example command, biastools will 
- [1] Simulate reads based on <ref.fa> and <vcf>, generating pair-end .fq.gz files for both haplotypes (work_dir/sample_name.hapA/B_1/2.fq.gz). 
- [2] Then align the reads to the reference <ref.fa>, generating bam file with phasing information (work_dir/sample_name.run_id.sorted.bam).
- [3] Analyze the bam file with context-aware asignment method, generating bias reports and plots.

#### Other aligners
Biastools support the alignment with `Bowtie 2` and `BWA MEM`. Additional alignment method can be performed on the simulated reads and feed in the analysis with command
```
$ Biastools.sh -3 -o <work_dir> -g <ref.fa> -v <vcf> -s 'sample_name' -r 'run_id'
```

Noted that the alignment file should be named with <work_dir/sample_name.run_id.sorted.bam> and tag with haplotype information.


#### Real data
Context-aware assignment can also analyze real data with `-R` option, but only generate the plot without simulation information (sample_id.real.indel_balance.pdf).
```
$ Biastools.sh -3 -o <work_dir> -g <ref.fa> -v <vcf> -s 'sample_name' -r 'run_id' -R
```


#### Multiple indel plots
Multiple analysis result can be combined into one single indel-balance plot.
```
$ Biastools.sh -3 -o <work_dir> -g <ref.fa> -v <vcf> -s 'sample_name' -r 'run_id' -l "file1.bias file2.bias file3.bias" -L "id1 id2 id3"
```


### Bias prediction on real data
```
$ Biastools.sh -4 -o <work_dir> -g <ref.fa> -v <vcf> -s 'sample_name' -r 'run_id' -R
```



