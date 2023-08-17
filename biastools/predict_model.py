import argparse
import numpy as np
import pandas as pd



def predict_bias(real_feature, miss_info, best_threshold, out_prefix):
    """
    quality score * balance score
    """
    real_feature['z_MAPQ'] = ((real_feature['AVG_MAPQ'] - 45) * -1).clip(lower=0)
    real_feature['combine_score'] = (real_feature['z_MAPQ']) * (real_feature['BALANCE']) #* (real_feature['BALANCE'])
    real_feature['plus_score']    = (real_feature['z_MAPQ']/45) + 1.5*real_feature['BALANCE']

    print(real_feature[real_feature['plus_score'] > best_threshold])
    real_feature[real_feature['plus_score'] > best_threshold].to_csv(out_prefix + "_bias.tsv", index=False, sep = "\t")




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-rr', '--real_report', help='the real data bias report')
    parser.add_argument('-thr', '--threshold',  help='the threshold for prediction model [1.5]', type=int, default=1.5)
    parser.add_argument('-out', '--out_prefix', help='the prefix for reports [predict]', type=str, default='predict')
    args = parser.parse_args()

    fn_real       = args.real_report
    best_th       = args.threshold
    out_prefix    = args.out_prefix
    
    df_real       = pd.read_csv(fn_real, sep='\t')

    # filter out the sites suspicious of imcomplete vcf information
    miss_info = (df_real['OTHER'] > df_real['NUM_READS'] * 0.9) + (df_real['OTHER'] > df_real['NUM_READS'] * 0.4) * \
                ( (df_real['REF'] == 0) + (df_real['ALT'] == 0 ))

    df_real[miss_info].to_csv(out_prefix + "_suspicious.tsv", index=False, sep = "\t")
    print("filtered number:", sum(miss_info))

    df_real_test  = df_real[~miss_info]
    predict_bias(df_real_test, miss_info, best_th, out_prefix)

    
