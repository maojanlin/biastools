import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

import math
import numpy as np


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

labels = ['>0.5', '0.3~0.5', '0.1~0.3', '0.05~0.1', '0.01~0.05', '<0.01']


def filter_bias(
    fn_bias     :str
    ) -> None:
    """
    filter the variant site with unbalanced mapping
    """
    gold_idx = 11
    even_idx = 9
    read_idx = 10
    allele_idx = 2
    depth_idx = 7
    ref_idx   = 3
    alt_idx   = 4
    g_ref_idx = 12
    g_alt_idx = 13

    f_b = open(fn_bias, "r")
    headline = f_b.readline().strip()
    list_data = []
    for line in f_b:
        fields = line.strip().split()
        data = [fields[0], int(fields[1]), float(fields[2])] + [int(ele) for ele in fields[3:9]] + [float(ele) for ele in fields[9:-2]] + [int(ele) for ele in fields[-2:]]
        list_data.append(data)
    f_b.close()

    list_biased = []
    list_artifact = []
    for fields in list_data:
        if abs(fields[read_idx] - fields[gold_idx]) > 0.07:
            list_biased.append(fields)
        if abs(fields[read_idx] - fields[allele_idx]) > 0.07:
            list_artifact.append(fields)
    list_intersect = list(set([tuple(ele) for ele in list_biased]).intersection(set([tuple(ele) for ele in list_artifact])))
    print("Num bias sites:", len(list_biased))
    print("Num artifact sites:", len(list_artifact))
    print("Num of intersection:", len(list_intersect))
    
    # ============ output =============
    f = open("list_bias.rpt", 'w')
    for data in list_biased:
        for ele in data:
            f.write(str(ele) + ',')
        f.write('\n')
    f.close()
    
    f = open("list_artifact.rpt", 'w')
    for data in list_artifact:
        for ele in data:
            f.write(str(ele) + ',')
        f.write('\n')
    f.close()
    
    f = open("list_exclude.rpt", 'w')
    list_exclude = list(set([tuple(ele) for ele in list_data]) - (set([tuple(ele) for ele in list_artifact])) - (set([tuple(ele) for ele in list_biased])))
    for data in list_exclude:
        for ele in data:
            f.write(str(ele) + ',')
        f.write('\n')
    f.close()
    # ============ output =============
    
    avg_depth = sum([ele[depth_idx] for ele in list_data])/len(list_data)
    avg_fail  = sum([ele[depth_idx] - ele[ref_idx] - ele[alt_idx] for ele in list_data])/len(list_data)
    print("==============================================================================")
    print("Average Read depth:", avg_depth)
    print("Average Fail read:",  avg_fail)
    print("==============================================================================")
    print("chr, pos, ref, alt, depth, allele_blance, gold_dist, others, eveness")
    for ele in list_biased:
        print(ele[0], ele[1], str(ele[g_ref_idx])+'->'+str(ele[ref_idx]), str(ele[g_alt_idx])+'->'+str(ele[alt_idx]), ele[depth_idx], ele[ref_idx]/(ele[ref_idx] + ele[alt_idx]), \
                ele[gold_idx], ele[depth_idx]-ele[ref_idx]-ele[alt_idx], ele[even_idx], sep='\t')
    #print("==============================================================================")
    #for ele in list_artifact:
    #    print(ele[0], ele[1], ele[depth_idx], ele[ref_idx]/(ele[ref_idx] + ele[alt_idx]), ele[depth_idx]-ele[ref_idx]-ele[alt_idx])

    #=========================== standard ref_bias to read_distribute plot ============================
    all_data = pd.DataFrame()
    all_data['READ DEPTH'] = [ele[depth_idx] / avg_depth for ele in list_data]
    all_data['ALLELIC BALANCE'] = [ele[ref_idx]/(ele[ref_idx]+ele[alt_idx]) for ele in list_data]
    all_data['NUM FAIL']   = [ele[depth_idx] - ele[ref_idx] -ele[alt_idx] for ele in list_data]
    all_data['EVEN VAR']   = [map_color(ele[even_idx]) for ele in list_data]
    all_data.head()

    plt.clf()
    #ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "NUM FAIL", data = all_data)#hue="size", size="size", data=tips)
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "EVEN VAR", data = all_data)#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="P-value", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([0,1])
    plt.show()

    #=========================== standard ref_bias to read_distribute plot ============================
    biased_data = pd.DataFrame()
    biased_data['READ DEPTH'] = [ele[depth_idx] / avg_depth for ele in list_biased]
    biased_data['ALLELIC BALANCE'] = [ele[ref_idx]/(ele[ref_idx]+ele[alt_idx]) for ele in list_biased]
    biased_data['NUM FAIL']   = [ele[depth_idx] - ele[ref_idx] -ele[alt_idx] for ele in list_biased]
    biased_data['EVEN VAR']   = [map_color(ele[even_idx]) for ele in list_biased]
    biased_data.head()

    plt.clf()
    #ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "NUM FAIL", data = biased_data)#hue="size", size="size", data=tips)
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "EVEN VAR", data = biased_data)#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="P-value", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([0,1])
    plt.show()

    #=========================== standard ref_bias to read_distribute plot ============================
    artifact_data = pd.DataFrame()
    artifact_data['READ DEPTH'] = [ele[depth_idx] / avg_depth for ele in list_artifact]
    artifact_data['ALLELIC BALANCE'] = [ele[ref_idx]/(ele[ref_idx]+ele[alt_idx]) for ele in list_artifact]
    artifact_data['NUM FAIL']   = [ele[depth_idx] - ele[ref_idx] -ele[alt_idx] for ele in list_artifact]
    artifact_data['EVEN VAR']   = [map_color(ele[even_idx]) for ele in list_artifact]
    artifact_data.head()

    plt.clf()
    #ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "NUM FAIL", data = artifact_data)#hue="size", size="size", data=tips)
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "EVEN VAR", data = artifact_data)#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="P-value", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([0,1])
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-mb', '--merge_bias_report', help='merged bias report, should contain the golden')
    args = parser.parse_args()
    
    fn_bias = args.merge_bias_report
    filter_bias(fn_bias)
