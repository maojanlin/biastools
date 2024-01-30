import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import math
import random
import numpy as np


colors = ["#bce4ff", "#8bd0fe", "#59bcfc", "#0099fc", "#0086dd", "#006bb1", "#004a7a", "#002740"]
colors = ["#f2dad5", "#e8bfc1", "#d9a4b2", "#c78ba6", "#aa719a", "#8b5b89", "#634271", "#3c2a4f"]

def map_mapq_to_size(mapq):
    if mapq >= 40:
        return 0
    elif mapq >= 30:
        return 1
    elif mapq >= 20:
        return 2
    elif mapq >= 10:
        return 3
    elif mapq >= 5:
        return 4
    elif mapq >= 3:
        return 5
    elif mapq >= 1:
        return 6
    return 7

labels = ['>40', '30~40', '20~30', '10~20', '5~10', '3~5', '1~3', '<1']

def map_color(
        var:float
        )-> int:
    """
    color_code = int(var/2)
    if color_code > 20:
        color_code = 20
    return color_code
    """
    if var > 0.5:
        return 0
    elif var > 0.3:
        return 1
    elif var > 0.1:
        return 2
    elif var > 0.05:
        return 3
    elif var > 0.01:
        return 4
    else:
        return 5

p_labels = ['>0.5', '0.3~0.5', '0.1~0.3', '0.05~0.1', '0.01~0.05', '<0.01']

def map_num_to_size(num):
    if num == 0:
        return 0
    elif num <= 3:
        return 1
    elif num <= 5:
        return 2
    elif num <= 10:
        return 3
    elif num <= 15:
        return 4
    elif num <= 20:
        return 5
    elif num <= 30:
        return 6
    return 7

n_labels = ['0', '1~3', '4~6', '7~10', '11~15', '16~20', '21~30', '>30']


def dist_origin(a, b):
    return math.dist((0,a) + (b,0))



def plot_golden(out_prefix, df_use):
    # Add columns
    mapQ   = list(df_use['AVG_MAPQ'])
    pValue = list(df_use['EVEN_P_VALUE'])
    
    sp = pd.DataFrame()
    sp['ALLELIC BALANCE']    = list(df_use['BALANCE'])
    sp['MAPPING BALANCE']    = list(df_use['MAP_BALANCE'])
    sp['SIMULATION BALANCE'] = list(df_use['SIM_BALANCE'])
    sp.head()
    
    mapped_mapQ = [map_mapq_to_size(q) for q in mapQ]
    mapped_p    = [map_color(p) for p in pValue]
    sp['Avg_MapQ_code'] = mapped_mapQ
    sp['Even_p_value']  = mapped_p
    sp['Assign_other']  = [map_num_to_size(n) for n in list(df_use['OTHER'])   ]
    sp['Map_other']     = [map_num_to_size(n) for n in list(df_use['MIS_MAP']) ]
    sp['MapQ'] = list(mapQ)

    #================== color map ====================
    set_mapQ_value = set(sp['Avg_MapQ_code'])
    color_mapQ = []
    for idx in sorted(set_mapQ_value):
        color_mapQ.append(colors[idx])
    
    set_misMap_value = set(sp['Map_other'])
    color_misMap = []
    for idx in sorted(set_misMap_value):
        color_misMap.append(colors[idx])
    
    #=========================== all merged plot ============================
    print("Ploting the Merged golden distribution Plot!")
    sp['Normalized Assignment Balance'] = list(df_use['BALANCE']-df_use['SIM_BALANCE']) # the average map_q score
    sp['Normalized Mapping Balance'] = list(df_use['MAP_BALANCE']-df_use['SIM_BALANCE']) # the average map_q score
    #ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Assignment Balance",  hue = "Avg_MapQ_code", data = sp, \
    #        xlim=(-0.6,0.6), ylim=(-0.6,0.6), palette=sns.color_palette(color_mapQ))
    ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Assignment Balance",  hue = "Map_other", data = sp, \
            xlim=(-0.8,0.8), ylim=(-0.8,0.8), palette=sns.color_palette(color_misMap))
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    #plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0, 0), loc='lower right', borderaxespad=0.2)
    plt.legend(h, n_labels, title="Mismapped Gain#", bbox_to_anchor=(0, 0), loc='lower right', borderaxespad=0.2)
    #plt.savefig(out_prefix + '.mismap.pdf')

    #print(df_use[sp['Normalized Assignment Balance']**2 + sp['Normalized Mapping Balance']**2 > 0.01])
    biased = (sp['Normalized Assignment Balance']**2 + sp['Normalized Mapping Balance']**2 > 0.01)
    b_loss = ((sp['Normalized Assignment Balance'] < sp['Normalized Mapping Balance']*2 + 0.1)*(sp['Normalized Assignment Balance']*2 + \
              0.1 > sp['Normalized Mapping Balance'])) + \
             ((sp['Normalized Assignment Balance'] + 0.1 > sp['Normalized Mapping Balance']*2)*(sp['Normalized Assignment Balance']*2 \
              < sp['Normalized Mapping Balance'] + 0.1))
    b_flux = (sp['Normalized Assignment Balance'] > 0.1)*(sp['Map_other'] >= 3) + \
             (sp['Normalized Assignment Balance'] < -0.1)*(sp['Map_other'] >= 3)
    b_artifact = (sp['Normalized Assignment Balance'] > 0.1)*(sp['Map_other'] < 3) + \
                 (sp['Normalized Assignment Balance'] < -0.1)*(sp['Map_other'] < 3)

    sp['Category'] = biased*4
    sp['Category'] -= (biased * b_loss)*3
    sp['Category'] -= (biased * ~b_loss * b_flux)*2
    sp['Category'] -= (biased * ~b_loss * b_artifact)*1
    labels = ['Balanced', 'Bias (Loss)', 'Bias (Flux)', 'Bias (Local)', 'Outliers']

    custom_palette = sns.color_palette('Set2')
    custom_palette = custom_palette[:4] + custom_palette[-1:]
    ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Assignment Balance",  hue = "Category", data = sp, \
            xlim=(-0.8,0.8), ylim=(-0.8,0.8), palette=custom_palette)
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    plt.legend(h, labels, title="Category#", bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0.2)
    plt.savefig(out_prefix + '.category.pdf')

    print("-------------------------------------------")
    print("Number of balanced:", sum(sp['Category'] == 0))
    print("Number of bias_loss:", sum(sp['Category'] == 1))
    print("Number of bias_flux:", sum(sp['Category'] == 2))
    print("Number of bias_local:", sum(sp['Category'] == 3))
    print("Number of outliers:", sum(sp['Category'] == 4))
    print("-------------------------------------------")

    df_use.loc[(sp['Category'] == 0).values, :].to_csv(out_prefix + '.balanced.tsv', index=False, sep="\t")
    df_use.loc[((sp['Category'] == 1)*(sp['Normalized Assignment Balance'] > 0)).values, :].to_csv(out_prefix + '.bias-loss.1.tsv', index=False, sep="\t")
    df_use.loc[((sp['Category'] == 1)*(sp['Normalized Assignment Balance'] < 0)).values, :].to_csv(out_prefix + '.bias-loss.2.tsv', index=False, sep="\t")
    df_use.loc[((sp['Category'] == 2)*(sp['Normalized Assignment Balance'] > 0)).values, :].to_csv(out_prefix + '.bias-flux.1.tsv', index=False, sep="\t")
    df_use.loc[((sp['Category'] == 2)*(sp['Normalized Assignment Balance'] < 0)).values, :].to_csv(out_prefix + '.bias-flux.2.tsv', index=False, sep="\t")
    df_use.loc[((sp['Category'] == 3)*(sp['Normalized Assignment Balance'] > 0)).values, :].to_csv(out_prefix + '.bias-local.1.tsv', index=False, sep="\t")
    df_use.loc[((sp['Category'] == 3)*(sp['Normalized Assignment Balance'] < 0)).values, :].to_csv(out_prefix + '.bias-local.2.tsv', index=False, sep="\t")
    df_use.loc[(sp['Category'] == 4).values, :].to_csv(out_prefix + '.bias-outlier.tsv', index=False, sep="\t")
    df_use.loc[(sp['Map_other'] > 4).values, :].to_csv(out_prefix + '.bias-mismap_gain.tsv', index=False, sep="\t")




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-mb', '--bias_report', help='bias report, must contain the golden information')
    parser.add_argument('-qt', '--quality_threshold', help='threshold that filtered the sites with avg_mapQ below the threshold', type=int)
    parser.add_argument('-out', '--output_prefix', help='the prefix for the output plots and report')
    args = parser.parse_args()
    
    fn_bias = args.bias_report
    mapQ_th = args.quality_threshold
    output_prefix = args.output_prefix
    if output_prefix == None:
        output_prefix = fn_bias

    df_use = pd.read_csv(fn_bias, sep='\t')
    if mapQ_th:
        df_use = df_use[df_use['AVG_MAPQ'] >= mapQ_th]
    df_use.head()

    plot_golden(output_prefix, df_use)

