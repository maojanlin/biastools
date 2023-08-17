import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import ListedColormap
import pandas as pd
import random

from sklearn import datasets, metrics
from sklearn.metrics import roc_curve, precision_recall_curve, auc



def get_label(df_simulation):
    """
    sort and label the simulated data, real data
    """
    sp = pd.DataFrame()
    sp['Map_other']                  = list(df_simulation['MIS_MAP'])
    sp['Normalized Allelic Balance'] = list(df_simulation['BALANCE']-df_simulation['SIM_BALANCE']) # the average map_q score
    sp['Normalized Mapping Balance'] = list(df_simulation['MAP_BALANCE']-df_simulation['SIM_BALANCE']) # the average map_q score
    
    biased = (sp['Normalized Allelic Balance']**2 + sp['Normalized Mapping Balance']**2 > 0.01)
    b_loss = ((sp['Normalized Allelic Balance'] < sp['Normalized Mapping Balance']*2 + 0.1) * \
              (sp['Normalized Allelic Balance']*2 + 0.1 > sp['Normalized Mapping Balance']))
    b_flux = (sp['Normalized Allelic Balance'] > 0.1)*(sp['Map_other'] > 4)
    b_artifact = (sp['Normalized Allelic Balance'] > 0.1)*(sp['Map_other'] <= 4)

    sp['Category'] = biased*4
    sp['Category'] -= (biased * b_loss)*3
    sp['Category'] -= (biased * ~b_loss * b_flux)*2
    sp['Category'] -= (biased * ~b_loss * b_artifact)*1
    
    sp['binary_category'] = (sp['Category'] > 0)
    return sp

    
def print_accuracy(predict, label):
    print("Correct Num:",    np.sum(predict == label))
    TP = np.sum((predict == label) * (predict != 0))
    FP = np.sum((predict != label) * (predict != 0))
    FN = np.sum((predict != label) * (predict == 0))
    print("True Positive:",  TP)
    print("False Positive:", FP)
    print("False Negative:", FN)
    print("Precision:", TP/(TP+FP))
    print("Recall:", TP/(TP+FN))


def combine_score(sim_feature, sim_label, real_feature, real_label, miss_info, best_threshold, out_prefix):
    """
    quality score * balance score
    """
    sim_feature['label'] = sim_label
    sim_feature['z_MAPQ'] = ((sim_feature['AVG_MAPQ'] - 45) * -1).clip(lower=0)
    sim_feature['combine_score'] = (sim_feature['z_MAPQ']) * (sim_feature['BALANCE']) #* (sim_feature['BALANCE'])
    sim_feature['plus_score']    = (sim_feature['z_MAPQ']/45) + 1.5*sim_feature['BALANCE']
    sim_feature['mix_score']     = sim_feature['plus_score'] + sim_feature['combine_score'] / 20

    fpr_m, tpr_m, thresholds = metrics.roc_curve(sim_feature['label'], sim_feature['combine_score'], pos_label=True)
    fpr_p, tpr_p, thresholds = metrics.roc_curve(sim_feature['label'], sim_feature['plus_score'], pos_label=True)
    plt.plot(fpr_m, tpr_m, label="simulation_mul, auc="+str(round(auc(fpr_m,tpr_m),2)))
    plt.plot(fpr_p, tpr_p, label="simulation_add, auc="+str(round(auc(fpr_p,tpr_p),2)))

    real_feature['label'] = real_label
    real_feature['z_MAPQ'] = ((real_feature['AVG_MAPQ'] - 45) * -1).clip(lower=0)
    real_feature['combine_score'] = (real_feature['z_MAPQ']) * (real_feature['BALANCE']) #* (real_feature['BALANCE'])
    real_feature['plus_score']    = (real_feature['z_MAPQ']/45) + 1.5*real_feature['BALANCE']
    r_fpr_m, r_tpr_m, thresholds = metrics.roc_curve(real_feature['label'], real_feature['combine_score'], pos_label=True)
    r_fpr_p, r_tpr_p, thresholds = metrics.roc_curve(real_feature['label'], real_feature['plus_score'], pos_label=True)
    plt.plot(r_fpr_m, r_tpr_m, label="real_mul, auc="+str(round(auc(r_fpr_m, r_tpr_m),2)))
    plt.plot(r_fpr_p, r_tpr_p, label="real_add, auc="+str(round(auc(r_fpr_p, r_tpr_p),2)))

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()
    plt.savefig(out_prefix + "_ROC.pdf")
    plt.clf()


    precision, recall, thresholds = precision_recall_curve(sim_feature['label'], sim_feature['combine_score'])
    precision_p, recall_p, thresholds = precision_recall_curve(sim_feature['label'], sim_feature['plus_score'])
    r_precision, r_recall, thresholds = precision_recall_curve(real_feature['label'], real_feature['combine_score'])
    r_precision_p, r_recall_p, thresholds = precision_recall_curve(real_feature['label'], real_feature['plus_score'])
    
    plt.plot(recall,   precision,   label="simulation_mul, auc="+str(round(auc(recall, precision),2)))
    plt.plot(recall_p, precision_p, label="simulation_add, auc="+str(round(auc(recall_p, precision_p),2)))
    plt.plot(r_recall,   r_precision,   label="real_mul, auc="+str(round(auc(r_recall, r_precision),2)))
    plt.plot(r_recall_p, r_precision_p, label="real_add, auc="+str(round(auc(r_recall_p, r_precision_p),2)))
    
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.savefig(out_prefix + "_PRC.pdf")
    plt.clf()

    print("====== sim featue ========")
    print_accuracy(sim_feature['plus_score'] > best_threshold, sim_feature['label'])
    print("====== real featue ========")
    print_accuracy(real_feature['plus_score'] > best_threshold, real_feature['label'])
    print("======= overlap =========")
    print_accuracy(sim_feature[~miss_info]['plus_score'] > best_threshold, real_feature['plus_score'] > 1.5)
    print("sim label True", np.sum(sim_feature['label']))
    print("sim feature",  np.sum(sim_feature['plus_score'] > best_threshold))
    print("real feature", np.sum(real_feature['plus_score'] > best_threshold))
    FP = (sim_feature['plus_score'] > best_threshold)* ~(sim_feature['label'])
    FN = (sim_feature['plus_score'] <= best_threshold)* (sim_feature['label'])
    return FP, FN




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sr', '--simulation_report', help='the simulation bias report')
    parser.add_argument('-rr', '--real_report', help='the real data bias report')
    parser.add_argument('-thr', '--threshold',  help='the threshold for prediction model [1.5]', type=int, default=1.5)
    parser.add_argument('-out', '--out_prefix', help='the prefix for plottings [predict]', type=str, default='predict')
    args = parser.parse_args()

    fn_simulation = args.simulation_report
    fn_real       = args.real_report
    best_th       = args.threshold
    out_prefix    = args.out_prefix
    
    df_simulation = pd.read_csv(fn_simulation, sep='\t')
    df_real       = pd.read_csv(fn_real, sep='\t')

    sp_label = get_label(df_simulation)

    # filter out the sites suspicious of imcomplete vcf information
    miss_info = (df_real['OTHER'] > df_real['NUM_READS'] * 0.9) + (df_real['OTHER'] > df_real['NUM_READS'] * 0.4) * \
                ( (df_real['REF'] == 0) + (df_real['ALT'] == 0 ))
    no_info = df_simulation['AVG_MAPQ'].isnull()
    no_info += df_simulation['MAP_BALANCE'].isnull()
    no_info += df_simulation['BALANCE'].isnull()
    miss_info += no_info
    df_simulation = df_simulation[~no_info]
    sp_label      = sp_label[~no_info]
    print("filtered number:", sum(miss_info))

    df_real_test  = df_real[~miss_info]
    sp_real_label = sp_label[~miss_info]
    FP, FN = combine_score(df_simulation, sp_label.iloc[:, 4].values, df_real_test, sp_real_label.iloc[:, 4].values, \
                           miss_info, best_th, out_prefix)
    
    # print data of FP and FN
    with pd.option_context('display.max_rows', None):  # more options can be specified also
        print("False Positive:")
        print(df_simulation[FP])
        print("==================================================")
        print("False Negative:")
        print(df_simulation[FN])

    # plot false positive
    labels = ['Balanced', 'Bias (Loss)', 'Bias (Flux)', 'Bias (Local)', 'Outliers']
    idx_cat = set(sp_label[FP+FN]["Category"])
    labels = [labels[idx] for idx in sorted(idx_cat)]
    ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Allelic Balance",  hue = "Category", data = sp_label[FP+FN], \
            xlim=(-0.8,0.8), ylim=(-0.8,0.8), palette='Set2')
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.fig.suptitle("False Positive and False Negative")
    ax.fig.tight_layout()
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    plt.legend(h, labels, title="Category#", bbox_to_anchor=(0, 0), loc='lower left', borderaxespad=0.2)
    plt.savefig(out_prefix + '_FP_and_FN.pdf')
    plt.clf()
    
    labels = ['Balanced', 'Bias (Loss)', 'Bias (Flux)', 'Bias (Local)', 'Outliers']
    idx_cat = set(sp_label[~(FP+FN)]["Category"])
    labels = [labels[idx] for idx in sorted(idx_cat)]
    ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Allelic Balance",  hue = "Category", data = sp_label[~(FP+FN)], \
            xlim=(-0.8,0.8), ylim=(-0.8,0.8), palette='Set2')
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.fig.suptitle("True Positive and True Negative")
    ax.fig.tight_layout()
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    plt.legend(h, labels, title="Category#", bbox_to_anchor=(0, 0), loc='lower left', borderaxespad=0.2)
    plt.savefig(out_prefix + '_TP_and_TN.pdf')
    
    

