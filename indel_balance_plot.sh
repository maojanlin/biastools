#python3 indel_balance_plot.py \
#    -lr  bias_report_octopus/wgs.octopus.vg_pop5.chr20.bias \
#         bias_report_octopus/wgs.octopus.vg_major.chr20.bias \
#         bias_report_octopus/wgs.octopus.vg_linear.chr20.bias \
#         bias_report_octopus/wgs.octopus.giraffe_pop5.chr20.bias \
#         bias_report_octopus/wgs.octopus.giraffe_major.chr20.bias \
#         bias_report_octopus/wgs.octopus.giraffe_linear.chr20.bias \
#         bias_report_octopus/wgs.octopus.bt2.chr20.bias \
#         bias_report_octopus/wgs.octopus.bwa.chr20.bias \
#    -ln  "VG_pop5" \
#         "VG_major" \
#         "VG_linear" \
#         "Giraffe_pop5" \
#         "Giraffe_major" \
#         "Giraffe_linear" \
#         "bowtie2" \
#         "BWA MEM" \
#    -vcf ../biastools/simulated_hg002/chr20_het.vcf.gz \
#    -bd 20 \
#    -out vg_octopus_giraffe_map


#python3 indel_balance_plot.py \
#    -lr  bias_report/wgs.giraffe_1KGP001.chr20.bias \
#         bias_report/wgs.giraffe_pop5.chr20.bias \
#         bias_report/wgs.giraffe_major.chr20.bias \
#         bias_report/wgs.giraffe_linear.chr20.bias \
#         ../biastools/simulated_hg002/wgs_result/simulated.wgs.chr20.bias \
#         ../biastools/simulated_hg002/wgs_result_bwa/wgs.bwa.chr20.bias \
#    -ln  "Giraffe_1KGP001" \
#         "Giraffe_pop5" \
#         "Giraffe_major" \
#         "Giraffe_linear" \
#         "bowtie2" \
#         "BWA MEM" \
#    -vcf ../biastools/simulated_hg002/chr20_het.vcf.gz \
#    -bd 20 \
#    -out adaptive

python3 indel_balance_plot.py \
    -lr  bias_report_overlap/wgs.giraffe_1KGP001.chr20.overlap.bias \
         bias_report_overlap/wgs.giraffe_pop5.chr20.overlap.bias \
         bias_report_overlap/wgs.giraffe_major.chr20.overlap.bias \
         bias_report_overlap/wgs.giraffe_linear.chr20.overlap.bias \
         bias_report_overlap/wgs.bt2.chr20.overlap.bias \
         bias_report_overlap/wgs.bt2_local.chr20.overlap.bias \
         bias_report_overlap/wgs.bwa.chr20.overlap.bias \
    -ln  "Giraffe_1KGP001" \
         "Giraffe_pop5" \
         "Giraffe_major" \
         "Giraffe_linear" \
         "bowtie2" \
         "bowtie2-local" \
         "BWA MEM" \
    -vcf ../biastools/simulated_hg002/chr20_het.vcf.gz \
    -bd 20 \
    -map \
    -out overlap_simulation_map

#python3 indel_balance_plot.py \
#    -lr  real_bias_report_overlap/real_data.giraffe_1KGP001.chr20.overlap.bias \
#         real_bias_report_overlap/real_data.giraffe_pop5.chr20.overlap.bias \
#         real_bias_report_overlap/real_data.giraffe_major.chr20.overlap.bias \
#         real_bias_report_overlap/real_data.giraffe_linear.chr20.overlap.bias \
#         real_bias_report_overlap/real_data.bt2.chr20.overlap.bias \
#         real_bias_report_overlap/real_data.bt2_local.chr20.overlap.bias \
#         real_bias_report_overlap/real_data.bwa.chr20.overlap.bias \
#    -ln  "Giraffe_1KGP001" \
#         "Giraffe_pop5" \
#         "Giraffe_major" \
#         "Giraffe_linear" \
#         "bowtie2" \
#         "bowtie2-local" \
#         "BWA MEM" \
#    -vcf ../biastools/simulated_hg002/chr20_het.vcf.gz \
#    -bd 20 \
#    -real \
#    -out overlap_real


#python3 indel_balance_plot.py \
#    -lr  real_bias_report/real_data.giraffe_1KGP001.chr20.bias \
#         real_bias_report/real_data.giraffe_pop5.chr20.bias \
#         real_bias_report/real_data.giraffe_major.chr20.bias \
#         real_bias_report/real_data.giraffe_linear.chr20.bias \
#         real_bias_report/real_data.bt2.chr20.bias \
#         real_bias_report/real_data.bwa.chr20.bias \
#    -ln  "Giraffe_1KGP001" \
#         "Giraffe_pop5" \
#         "Giraffe_major" \
#         "Giraffe_linear" \
#         "bowtie2" \
#         "BWA MEM" \
#    -vcf ../biastools/simulated_hg002/chr20_het.vcf.gz \
#    -bd 20 \
#    -real \
#    -out giraffe_real
