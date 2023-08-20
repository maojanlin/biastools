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

def map_waste_to_color(value):
    return int(math.ceil(value*8))


def plot_golden(out_prefix, df_use):
    
    # Add columns
    df_use['WASTE_INFO']   = (df_use['OTHER'])/(df_use['NUM_READS']+0.01)
    mapQ   = list(df_use['AVG_MAPQ'])
    pValue = list(df_use['EVEN_P_VALUE'])
    
    sp = pd.DataFrame()
    sp['ASSIGNMENT BALANCE'] = list(df_use['BALANCE'])
    sp['MAPPING BALANCE']    = list(df_use['MAP_BALANCE'])
    sp['SIMULATION BALANCE'] = list(df_use['SIM_BALANCE'])
    sp.head()
    
    mapped_mapQ = [map_mapq_to_size(q) for q in mapQ]
    mapped_p    = [map_color(p) for p in pValue]
    waste_value = [map_waste_to_color(q) for q in list(df_use['WASTE_INFO'])]
    sp['Avg_MapQ_code'] = mapped_mapQ
    sp['Even_p_value']  = mapped_p
    sp['Waste_value']   = waste_value
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
    
    #=========================== standard ref_bias to read_distribute plot ============================
    print("Ploting the Standard Ref Bias Plot!")
    plt.clf()
    ax = sns.scatterplot(y="ASSIGNMENT BALANCE", x="MAPPING BALANCE",  hue = "Avg_MapQ_code", data = sp, palette=sns.color_palette(color_mapQ))
    #ax = sns.scatterplot(y="ASSIGNMENT BALANCE", x="MAPPING BALANCE",  hue = "Even_p_value", data = sp)#hue="size", size="size", data=tips)
    #ax = sns.scatterplot(y="ASSIGNMENT BALANCE", x="MAPPING BALANCE",  hue = "Waste_value", data = sp)#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,1])
    plt.ylim([0,1])
    
    FN_FIG = out_prefix + '.diff-assign2map_dot.pdf'
    plt.savefig(FN_FIG)
    
    #=========================== golden to read_distribute plot ============================
    print("Ploting the Golden distribution Plot!")
    plt.clf()
    ax = sns.scatterplot(x="SIMULATION BALANCE", y="MAPPING BALANCE",  hue = "Avg_MapQ_code", data = sp, palette=sns.color_palette(color_mapQ))#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,1])
    plt.ylim([0,1])
    
    FN_FIG = out_prefix + '.diff-sim2map_dot.pdf'
    plt.savefig(FN_FIG)
    
    #=========================== golden to ref_bias plot ============================
    plt.clf()
    ax = sns.scatterplot(x="SIMULATION BALANCE", y="ASSIGNMENT BALANCE",  hue = "Avg_MapQ_code", data = sp, palette=sns.color_palette(color_mapQ))#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,1])
    plt.ylim([0,1])
    
    FN_FIG = out_prefix + '.diff-sim2assign_dot.pdf'
    plt.savefig(FN_FIG)
    
    #=========================== all merged plot ============================
    print("Ploting the Merged golden distribution Plot!")
    plt.clf()
    sp['Normalized Assignment Balance'] = list(df_use['BALANCE']-df_use['SIM_BALANCE']) # the average map_q score
    sp['Normalized Mapping Balance'] = list(df_use['MAP_BALANCE']-df_use['SIM_BALANCE']) # the average map_q score
    ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Assignment Balance",  hue = "Avg_MapQ_code", data = sp, \
            xlim=(-0.8,0.8), ylim=(-0.8,0.8), palette=sns.color_palette(color_mapQ))
    #ax = sns.jointplot(x="Normalized Mapping Balance", y="Normalized Assignment Balance",  hue = "Map_other", data = sp, \
    #        xlim=(-0.8,0.8), ylim=(-0.8,0.8), palette=sns.color_palette(color_misMap))
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0, 0), loc='lower left', borderaxespad=0.2)
    #plt.legend(h, n_labels, title="Mismapped Gain#", bbox_to_anchor=(1,0), loc='lower right', borderaxespad=0.2)
    
    FN_FIG = out_prefix + '.category.MapQ.pdf'
    plt.savefig(FN_FIG)
    
    #======================= allelic difference plot =========================
    plt.clf()
    list_ref_diff = list(df_use['REF']-df_use['SIM_REF'])
    list_alt_diff = list(df_use['ALT']-df_use['SIM_ALT'])
    for idx in range(len(list_ref_diff)):
        list_ref_diff[idx] += random.uniform(-0.3, 0.3) # scatter plot
        list_alt_diff[idx] += random.uniform(-0.3, 0.3)
    sp['Ref# - Simulation Ref#'] = list_ref_diff
    sp['Alt# - Simulation Alt#'] = list_alt_diff
    
    #ax = sns.jointplot(x="Ref# - Simulation Ref#", y="Alt# - Simulation Alt#",  hue = "Even_p_value", data = sp, xlim=(-20,20), ylim=(-20,15))
    ax = sns.jointplot(x="Ref# - Simulation Ref#", y="Alt# - Simulation Alt#",  hue = "Avg_MapQ_code", data = sp, xlim=(-30,30), ylim=(-30,15), palette=sns.color_palette(color_mapQ))
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    #plt.legend(h, p_labels, title="Even P Value", bbox_to_anchor=(0,1), loc='upper right')
    plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(1,1), loc='upper right')
    
    FN_FIG = out_prefix + '.diff2-assign2sim.pdf'
    plt.savefig(FN_FIG)
    

    plt.clf()
    #ax = sns.jointplot(x="Ref# - Simulation Ref#", y="Alt# - Simulation Alt#",  hue = "Map_other", data = sp, xlim=(-30,30), ylim=(-30,15), palette=sns.color_palette(color_misMap))
    #ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
    #ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
    #ax.ax_joint.get_legend().remove()
    #h, l = ax.ax_joint.get_legend_handles_labels()
    #plt.legend(h, n_labels, title="Mismapped Gain#", bbox_to_anchor=(1,1), loc='upper right')
    #
    #FN_FIG = out_prefix + '-read_diff_allelic.mismap.pdf'
    #plt.savefig(FN_FIG)
    #====================== mapping difference plot =========================
    plt.clf()
    list_m_ref_diff = list(df_use['MAP_REF']-df_use['SIM_REF'])
    list_m_alt_diff = list(df_use['MAP_ALT']-df_use['SIM_ALT'])
    for idx in range(len(list_m_ref_diff)):
        list_m_ref_diff[idx] += random.uniform(-0.3, 0.3) # scatter plot
        list_m_alt_diff[idx] += random.uniform(-0.3, 0.3)
    sp['Mapping Ref# - Simulation Ref#'] = list_m_ref_diff
    sp['Mapping Alt# - Simulation Alt#'] = list_m_alt_diff
    
    #ax = sns.jointplot(x="Mapping Ref# - Simulation Ref#", y="Mapping Alt# - Simulation Alt#",  hue = "Even_p_value", data = sp, xlim=(-20,20), ylim=(-20,15))
    ax = sns.jointplot(x="Mapping Ref# - Simulation Ref#", y="Mapping Alt# - Simulation Alt#",  hue = "Avg_MapQ_code", data = sp, xlim=(-30,30), ylim=(-30,15), palette=sns.color_palette(color_mapQ))
    ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
    ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
    ax.ax_joint.get_legend().remove()
    h, l = ax.ax_joint.get_legend_handles_labels()
    #plt.legend(h, p_labels, title="Even P Value", bbox_to_anchor=(0,1), loc='upper right')
    plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(1,1), loc='upper right')
    
    FN_FIG = out_prefix + '.diff2-map2sim.pdf'
    plt.savefig(FN_FIG)
    
    
    plt.clf()
    #ax = sns.jointplot(x="Mapping Ref# - Simulation Ref#", y="Mapping Alt# - Simulation Alt#",  hue = "Map_other", data = sp, xlim=(-30,30), ylim=(-30,15), palette=sns.color_palette(color_misMap))
    #ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
    #ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
    #ax.ax_joint.get_legend().remove()
    #h, l = ax.ax_joint.get_legend_handles_labels()
    #plt.legend(h, n_labels, title="Mismapped Gain#", bbox_to_anchor=(1,1), loc='upper right')
    #
    #FN_FIG = out_prefix + '-read_diff_mapping.mismap.pdf'
    #plt.savefig(FN_FIG)
    #======================== read loss-gain plot ===========================
    plt.clf()
    array_m_ref_diff = -np.array(df_use['MAP_REF']-df_use['SIM_REF'])
    array_m_alt_diff = -np.array(df_use['MAP_ALT']-df_use['SIM_ALT'])
    list_read_loss = list(np.where(array_m_ref_diff < 0, 0, array_m_ref_diff) + np.where(array_m_alt_diff < 0, 0, array_m_alt_diff))
    list_read_gain = list(df_use["MIS_MAP"])
    for idx in range(len(list_m_ref_diff)):
        list_read_loss[idx] += random.uniform(0,0.5) # scatter plot
        list_read_gain[idx] += random.uniform(0,0.5) 
    sp["Loss of Read (Ref + Alt)"] = list_read_loss
    sp["Gain of Read"]             = list_read_gain

    #ax = sns.jointplot(x="Loss of Read (Ref + Alt)", y="Gain of Read", hue = "Avg_MapQ_code", data = sp, xlim=(0,30), ylim=(0,30), palette=sns.color_palette(color_mapQ))
    #ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
    #ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
    #ax.ax_joint.get_legend().remove()
    #h, l = ax.ax_joint.get_legend_handles_labels()
    #plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(1,1), loc='upper right')
    #
    #FN_FIG = out_prefix + '-loss_gain.pdf'
    #plt.savefig(FN_FIG)

    plt.close("all")
    sns.color_palette()
    ref_loss  = list(df_use['SIM_REF']-df_use['MAP_REF'])
    alt_loss  = list(df_use['SIM_ALT']-df_use['MAP_ALT'])
    read_gain = list(df_use["MIS_MAP"]) 
    hist_data = pd.DataFrame()
    hist_data['loss/gain in a variant'] = read_gain + ref_loss + alt_loss
    hist_data['category'] = ["MisMap gain"]*len(read_gain) + ["Ref loss"]*len(ref_loss) + ["Alt loss"]*len(alt_loss)

    bin_num = max(hist_data['loss/gain in a variant']) - min(hist_data['loss/gain in a variant'])
    plt.clf()
    ax = sns.displot(hist_data, x="loss/gain in a variant", bins=bin_num, hue="category", log_scale=(False,True), element="step")
    ax.set(ylabel="occurence")
    FN_FIG = out_prefix + '.loss_gain_occurence.pdf'
    plt.savefig(FN_FIG)
    
    #plt.clf()
    #ax = sns.displot(hist_data, x="loss/gain in a variant", bins=int(bin_num/3), hue="category", log_scale=(False,True), multiple="dodge")
    #ax.set(ylabel="occurence")
    #FN_FIG = out_prefix + '-loss_gain_occurence.dodge.pdf'
    #plt.savefig(FN_FIG)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-mb', '--bias_report', help='bias report, must contain the golden information')
    parser.add_argument('-qt', '--quality_threshold', help='threshold that filtered the sites with avg_mapQ below the threshold', type=int, default=0)
    parser.add_argument('-out', '--output_prefix', help='the prefix for the output plots and report')
    args = parser.parse_args()
    
    fn_bias = args.bias_report
    mapQ_th = args.quality_threshold
    output_prefix = args.output_prefix
    if output_prefix == None:
        output_prefix = fn_bias

    df_use = pd.read_csv(fn_bias, sep='\t')
    df_use = df_use[df_use['AVG_MAPQ'] >= mapQ_th]
    df_use.head()

    plot_golden(output_prefix, df_use)

