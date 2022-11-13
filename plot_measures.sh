bias_report="simulated_hg002/wgs_result_bwa/wgs.bwa.chr20.bias"
bias_report_real="simulated_hg002/wgs_result_bwa/real_data.bwa.chr20.bias"
bias_report_naive="simulated_hg002/wgs_result_bwa/wgs.bwa.naive.chr20.bias"
bias_het_vcf="simulated_hg002/chr20_het.vcf.gz"

# making the report of different categories of bias
python3 golden_graph_report.py   -mb ${bias_report}".SNP"
python3 golden_graph_report.py   -mb ${bias_report}".gap"
# plot the measures with NMB and NAB
python3 golden_graph.py          -mb ${bias_report}".SNP"
python3 golden_graph.py          -mb ${bias_report}".gap"
# generate the balance files without golden / show the mapQ and balance difference between bias sites if -SIM is specified
python3 depth_balance_grapher.py -mb ${bias_report}".SNP" -sim
python3 depth_balance_grapher.py -mb ${bias_report}".gap" -sim

if [ ! -z ${bias_report_naive+x} ]; then 
    # plotting the sites naive and two-stage assignment having more than 0.1 balance different
    python3 assignment_compare.py    -ar ${bias_report}".SNP" -nr   ${bias_report_naive}".SNP"
    # plotting the indel balance plot
    python3 indel_balance_plot.py    -ar ${bias_report} -nr ${bias_report_naive} -vcf ${bias_het_vcf} 
fi

if [ ! -z ${bias_report_real+x} ]; then 
    python3 compare_real.py                 -sr ${bias_report}".SNP" -rr ${bias_report_real}".SNP"
    python3 threshold_with_combine_score.py -sr ${bias_report}".SNP" -rr ${bias_report_real}".SNP" 
fi
