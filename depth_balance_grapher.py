import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import math
import random
import numpy as np


#============================= mapping function =================================
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

def map_waste_to_color(value):
    return int(math.ceil(value*8))



def plot_statistics(fn_bias, all_data):
    #============================= STEP 1 =================================
    fn_fig = fn_bias + 'allelic_balance_graphs.pdf'
    title = 'Allelic Balance Distribution' 
    
    rb = all_data['ALLELIC BALANCE']
    if 'NaN' in rb:
        print("true")
    #: show stats
    print (rb.quantile(q=[.01, .1, 0.25,0.5,0.75, .90, .99]))
    print("mean: ", np.mean(rb))
    
    plt.clf()
    n, bins, patches = plt.hist(rb, bins = list(np.linspace(0,1,26)))
    plt.ylim((1,40000))
    plt.xlim((0,1))
    plt.yscale('log')
    plt.xlabel('ALLELIC BALANCE')
    plt.ylabel('Counts')
    plt.title(title)
    plt.savefig(fn_fig)
    
    
    #============================= STEP 2 =================================
    fn_fig = fn_bias + '-allelic_balance_eighty.pdf'#'SRR622457-grch37-bt2.pdf'
    title = 'Reference Bias Distribution of HET sites with Reference Bias >= 0.80'
    
    df_eighty = all_data[all_data['ALLELIC BALANCE'] >= 0.80]
    df_eighty_refbi = df_eighty['ALLELIC BALANCE']
    #print(len(df_eighty_refbi))
    plt.clf()
    plt.hist(df_eighty_refbi)
    plt.yscale('log')
    plt.xlabel('ALLELIC BALANCE')
    plt.ylabel('Counts')
    plt.title(title)
    plt.savefig(fn_fig)
    


def plot_balance(
    fn_bias     :str,
    flag_sim    :bool
    ) -> None:
    """
    filter the variant site with unbalanced mapping
    """
    df_use = pd.read_csv(fn_bias, sep='\t')
    #df = pd.read_csv(FN_INPUT, sep='\t')
    #df_mid = df[df['ALT_COUNT'] + df['REF_COUNT']> 0]
    #df_use = df_use[df_use['NUM_READS'] >= 10]
    ####df_use = df_use[df_use['HET_SITE'].isin(dict_exclude['chr21'])]
    #df_use.loc[:,'THRESHOLD'] = 15
    df_use.head()
    
    rb = df_use['REFERENCE_BIAS']
    df_use['AVERAGE_MAPQ'] = df_use['SUM_MAPQ']/(df_use['NUM_READS']) # the average map_q score
    df_use['WASTE_INFO']   = (df_use['BOTH_COUNT']+df_use['NEITHER_COUNT'])/(df_use['NUM_READS']+0.01)
    mapQ = list(df_use['AVERAGE_MAPQ'])
    pValue = list(df_use['EVEN_P_VALUE'])
    mapped_mapQ = [map_mapq_to_size(q) for q in mapQ]
    mapped_p    = [map_color(p) for p in pValue]
    waste_value = [map_waste_to_color(q) for q in list(df_use['WASTE_INFO'])]
    
    avg_depth = np.mean(list(df_use['NUM_READS']))
    std_depth = np.std( list(df_use['NUM_READS']))
    avg_fail  = np.mean(list(df_use['BOTH_COUNT']+df_use['NEITHER_COUNT']))
    print("==============================================================================")
    print("Average Read depth:", avg_depth)
    print("Average Fail read:",  avg_fail)
    print("==============================================================================")
    
    all_data = pd.DataFrame()
    all_data['ALLELIC BALANCE']      = list(rb)
    all_data['READ DISTRIBUTION']    = list(df_use['READ_DISTRIBUTION'])
    all_data['x AVG READ DEPTH'] = list(df_use['NUM_READS']/avg_depth)
    all_data['AVERAGE MAPQ']         = list(df_use['AVERAGE_MAPQ'])
    
    all_data['Avg_MapQ_code'] = mapped_mapQ
    all_data['Even_p_value']  = mapped_p
    all_data['Waste_value']   = waste_value
    all_data['MapQ'] = list(mapQ)
    if flag_sim:
        all_data['GOLDEN_DISTRIBUTION'] = list(df_use['GOLDEN_DISTRIBUTION'])
        all_data['READ_DISTRIBUTION']   = list(df_use['READ_DISTRIBUTION'])
    all_data.head()

    #=========================== standard ref_bias to read_distribute plot ============================
    plt.clf()
    title = "Normalized Read Depth - Allelic Balance"
    #ax = sns.scatterplot(y="ALLELIC BALANCE", x="READ DEPTH",  hue = "NUM FAIL", data = all_data)#hue="size", size="size", data=tips)
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="x AVG READ DEPTH",  hue = "Even_p_value", data = all_data)#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, p_labels, title="P-value", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([-0.2,1])
    #plt.show()
    fn_fig = fn_bias + '-var2depth_EvenP.pdf'
    plt.title(title)
    plt.savefig(fn_fig)
    
    plt.clf()
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="x AVG READ DEPTH",  hue = "Avg_MapQ_code", data = all_data)#hue="size", size="size", data=tips)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="MapQ", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([-0.2,1])
    #plt.show()
    fn_fig = fn_bias + '-var2depth_MapQ.pdf'
    plt.title(title)
    plt.savefig(fn_fig)

    return all_data


def categorizing_simulated_data(all_data):
    biased_data    = all_data[abs(all_data['READ_DISTRIBUTION'] - all_data['GOLDEN_DISTRIBUTION']) > 0.07]
    artifact_data  = all_data[abs(all_data['ALLELIC BALANCE']   - all_data['READ_DISTRIBUTION'])   > 0.07]
    intersect_data = pd.merge(biased_data, artifact_data, how ='inner')

    print("==============================================================================")
    print("# of biased site:", len(biased_data))
    print("# of artifact site:", len(artifact_data))
    print("# of intersections:", len(intersect_data))
    print("==============================================================================")

    # print all data
    plt.clf()
    title = "All Sites"
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="x AVG READ DEPTH",  hue = "Avg_MapQ_code", data = all_data)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="MapQ", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([0,1])
    plt.title(title)
    plt.show()

    # print biased site
    plt.clf()
    title = "Biased Sites (Discordance between simulation and mapping)"
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="x AVG READ DEPTH",  hue = "Avg_MapQ_code", data = biased_data)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="MapQ", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([0,1])
    plt.title(title)
    plt.show()
    
    # print artifacts
    plt.clf()
    title = "Artifacts (Discordance between mapping and base-call)"
    ax = sns.scatterplot(y="ALLELIC BALANCE", x="x AVG READ DEPTH",  hue = "Avg_MapQ_code", data = artifact_data)
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, labels, title="MapQ", bbox_to_anchor=(0.87, 1), loc=2, borderaxespad=0., framealpha=1)
    plt.xlim([0,2.5])
    plt.ylim([0,1])
    plt.title(title)
    plt.show()

    



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-mb', '--merge_bias_report', help='merged bias report, should contain the golden')
    parser.add_argument('-sim', '--simulated_flag', action='store_true', help='specify to evoke categorizing plot')
    args = parser.parse_args()
    
    fn_bias  = args.merge_bias_report
    flag_sim = args.simulated_flag

    parsed_data = plot_balance(fn_bias, flag_sim)
    plot_statistics(fn_bias, parsed_data)
    if flag_sim:
        categorizing_simulated_data(parsed_data)


