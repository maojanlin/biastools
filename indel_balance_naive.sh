python3 indel_balance_naive.py -ar bias_report_overlap/wgs.bwa.chr20.overlap.bias \
                               -nr bias_report_naive/wgs.bwa.chr20.naive.bias \
                               -vcf chr20_het.vcf.gz \
                               -bd 20 \
                               -out naive_bwa

python3 indel_balance_naive.py -ar bias_report_overlap/wgs.bt2.chr20.overlap.bias \
                               -nr bias_report_naive/wgs.bt2.chr20.naive.bias \
                               -vcf chr20_het.vcf.gz \
                               -bd 20 \
                               -out naive_bt2

python3 indel_balance_naive.py -ar bias_report_overlap/wgs.bt2_local.chr20.overlap.bias \
                               -nr bias_report_naive/wgs.bt2_local.chr20.naive.bias \
                               -vcf chr20_het.vcf.gz \
                               -bd 20 \
                               -out naive_bt2_local

python3 indel_balance_naive.py -ar bias_report_overlap/wgs.giraffe_1KGP001.chr20.overlap.bias \
                               -nr bias_report_naive/wgs.giraffe_1KGP001.chr20.naive.bias \
                               -vcf chr20_het.vcf.gz \
                               -bd 20 \
                               -out naive_giraffe_1KGP001
