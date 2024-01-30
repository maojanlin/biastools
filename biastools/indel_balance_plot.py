import argparse
import seaborn as sns
import matplotlib.pyplot as plt

import math
import numpy as np
import pysam



def read_bias_report(fn_bias_report):
    list_bias_SNP = []
    list_bias_gap = []
    f = open(fn_bias_report, 'r')
    header = f.readline()
    for line in f:
        fields = line.split()
        if fields[-1] == '.':
            list_bias_gap.append(fields)
        else:
            list_bias_SNP.append(fields)
    f.close()
    return list_bias_SNP, list_bias_gap


def calculate_SNP_balance(assign_SNP, flag_real):
    """
    Return for simulated read:
    [[simulate_balance], [map_balance], [assign_balance]]
    Return for real read:
    [assign_balance]
    """
    if flag_real:
        record = [[float(fields[5]) for fields in assign_SNP]]
    else:
        record = [[],[],[]]
        for idx in range(len(assign_SNP)):
            record[0].append(float(assign_SNP[idx][14]))
            record[1].append(float(assign_SNP[idx][10]))
            record[2].append(float(assign_SNP[idx][5]))
    return record


def calculate_gap_balance(assign_gap, f_vcf, len_bd, get_idx):
    list_insert = [ [] for _ in range(len_bd) ]
    list_delete = [ [] for _ in range(len_bd) ]
    for idx in range(len(assign_gap)):
        ref_name  = assign_gap[idx][0]
        var_start = int(assign_gap[idx][1])
        var_segment = f_vcf.fetch(contig=ref_name, start=var_start-1, stop=var_start+1) # get exactly the variant at the site
        for var in var_segment:
            if var.start+1 != var_start:
                continue
            len_ref = len(var.ref)
            if len(var.alts) == 1:
                len_alt = len(var.alts[0])
            else:
                hap = var.samples[0]['GT']
                if hap[0] != 0:
                    len_alt = len(var.alts[hap[0]-1])
                else:
                    len_alt = len(var.alts[hap[1]-1])
            
            if len_ref > len_alt: # deletion
                diff = min(len_ref - len_alt -1, len_bd-1)
                record = float(assign_gap[idx][get_idx])
                list_delete[diff].append(record)
            else: # 0 and insertions
                diff = min(len_alt - len_ref -1, len_bd-1)
                record = float(assign_gap[idx][get_idx])
                list_insert[diff].append(record)
    return list_insert, list_delete


def addlabels(x,y,len_bd):
    for i in range(len(x)):
        plt.text(i-len_bd, y[i], y[i], ha = 'center')


def plot_balance(balance_delete, balance_SNP, balance_insert, output_name, len_bd, list_incidents, list_plot_name):
    len_plot = len(list_plot_name)
    balance_list = [np.zeros(2*len_bd+1) for idx in range(len_plot)]
    balance_1st  = [np.zeros(2*len_bd+1) for idx in range(len_plot)]
    balance_3rd  = [np.zeros(2*len_bd+1) for idx in range(len_plot)]
    
    for idy, list_delete in enumerate(balance_delete):
        for idx in range(len_bd):
            list_balance = np.array(list_delete[idx])
            if len(list_balance) > 1:
                median = np.median(list_balance[~np.isnan(list_balance)])
                balance_list[idy][len_bd-1-idx] = median
                balance_1st [idy][len_bd-1-idx] = median - np.quantile(list_balance[~np.isnan(list_balance)], 0.25)
                balance_3rd [idy][len_bd-1-idx] = np.quantile(list_balance[~np.isnan(list_balance)], 0.75) - median
            else:
                balance_list[idy][len_bd-1-idx] = np.nan
                balance_1st [idy][len_bd-1-idx] = np.nan
                balance_3rd [idy][len_bd-1-idx] = np.nan
    
    for idy, list_balance in enumerate(np.array(balance_SNP)):
        median = np.median(list_balance[~np.isnan(list_balance)])
        balance_list[idy][len_bd] = median
        balance_1st [idy][len_bd] = median - np.quantile(list_balance[~np.isnan(list_balance)], 0.25)
        balance_3rd [idy][len_bd] = np.quantile(list_balance[~np.isnan(list_balance)], 0.75) - median
    
    for idy, list_insert in enumerate(balance_insert):
        for idx in range(len_bd):
            list_balance = np.array(list_insert[idx])
            if len(list_balance) > 1:
                median = np.median(list_balance[~np.isnan(list_balance)])
                balance_list[idy][len_bd+1+idx] = median
                balance_1st [idy][len_bd+1+idx] = median - np.quantile(list_balance[~np.isnan(list_balance)], 0.25)
                balance_3rd [idy][len_bd+1+idx] = np.quantile(list_balance[~np.isnan(list_balance)], 0.75) - median
            else:
                balance_list[idy][idx+len_bd+1] = np.nan
                balance_1st[idy][idx+len_bd+1] = np.nan
                balance_3rd[idy][idx+len_bd+1] = np.nan
    """
    for idx in range(len(balance_list[0])):
        print(idx, balance_list[0][idx], balance_1st[0][idx], balance_3rd[0][idx])
    print(balance_list[0])
    print(balance_list[1])
    print(balance_list[2])
    print(balance_list[3])
    """
    t = list(range(-len_bd,len_bd+1))
    f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3,1]})
    f.set_size_inches(15,8)
    #f.figsize = (15,13)
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    #colors = colors[:2] + colors[5:]
    #colors = colors[1:]
    for idx, name in enumerate(list_plot_name):
        a0.errorbar(t, (1-balance_list[idx]), yerr=(balance_3rd[idx], balance_1st[idx]), capsize=3, fmt='-o', label=name, color=colors[idx])
    
    a0.legend()
    #a0.set_ylim([0.3, 0.7])
    a0.axhline(y=0.5, color='gray', linestyle='dashdot', linewidth=0.9)
    a0.set(ylabel='Fraction of alternate allele')
    
    a1.set(xlabel='Insertion (+) or deletion (-) length')
    a1.set(ylabel='# of variants')
    a1.bar(t, list_incidents, align='center', width=0.5, log=True)
    a1.set_ylim([1, max(list_incidents)*5])
    addlabels(t, list_incidents, len_bd)
    #a1.set_yscale('log')
    a1.grid(axis='y', color='gray', linestyle='dashdot', linewidth=0.6)
    plt.savefig(output_name + '.indel_balance.pdf')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-lr', '--list_report', nargs='+', required=True, help='the list of assignment bias report')
    parser.add_argument('-ln', '--list_name',   nargs='+', required=True, help='the second bias report')
    parser.add_argument('-vcf', '--vcf_report', help='the vcf report for the bias report regions')
    parser.add_argument('-bd', '--boundary', type=int, default=40, help='the boundary indel lengths extend from 0')
    parser.add_argument('-map', '--flag_mapping', action='store_true', help='show the mapping rather than local result')
    parser.add_argument('-real', '--flag_real', action='store_true', help='specify if the report contains no simulation information')
    parser.add_argument('-out', '--output_name', help="output file name")
    args = parser.parse_args()

    list_report = args.list_report
    list_name   = args.list_name
    fn_vcf = args.vcf_report
    boundary = args.boundary
    flag_map = args.flag_mapping
    flag_real = args.flag_real
    output_name = args.output_name
    if output_name == None:
        output_name = list_name[0]

    assert len(list_report) == len(list_name), "Number of bias_report and bias names are different."
    
    f_vcf = pysam.VariantFile(fn_vcf)
    # read the bias report
    list_bias_report = []
    for fn_assign_report in list_report:
        assign_report = read_bias_report(fn_assign_report)
        list_bias_report.append(assign_report)
    
    # fetch the SNP balance information
    list_balance_SNP = []
    for assign_SNP, assign_gap in list_bias_report:
        balance_SNP = calculate_SNP_balance(assign_SNP, flag_real)
        list_balance_SNP.append(balance_SNP)

    if flag_real: # no simulation of mapping information provided
        list_plot_name = [name + '(real)' for name in list_name]
        
        # fetch the gap balance information
        list_balance_delete = []
        list_balance_insert = []
        for assign_SNP, assign_gap in list_bias_report:
            balance_insert, balance_delete = calculate_gap_balance(assign_gap, f_vcf, boundary, 5)
            list_balance_insert.append(balance_insert)
            list_balance_delete.append(balance_delete)

        balance_SNP    = list_balance_SNP
        balance_delete = list_balance_delete
        balance_insert = list_balance_insert
    else: # to plot the simulated reads, the first entry is the simulated balance information, then we can choose map or local_assignment
        flag_choice = 2
        gap_choice  = 5
        list_plot_name = ["simulated"]
        if flag_map:
            flag_choice = 1
            gap_choice  = 10
            list_plot_name += [name + '(map)' for name in list_name]
        else:
            list_plot_name += [name + '(assign)' for name in list_name]
    
        # fetch the gap balance information
        balance_insert, balance_delete = calculate_gap_balance(list_bias_report[0][1], f_vcf, boundary, 14) # getting the simulated information
        list_balance_delete = [balance_delete]
        list_balance_insert = [balance_insert]
        for assign_SNP, assign_gap in list_bias_report:
            balance_insert, balance_delete = calculate_gap_balance(assign_gap, f_vcf, boundary, gap_choice)
            list_balance_insert.append(balance_insert)
            list_balance_delete.append(balance_delete)

        balance_SNP    = [list_balance_SNP[0][0]] + [balance[flag_choice] for balance in list_balance_SNP]
        balance_delete = list_balance_delete
        balance_insert = list_balance_insert
    

    # get the incident numbers of the indels
    list_incidents = [len(balance) for balance in list_balance_delete[0]][::-1] + [len(list_balance_SNP[0][0])] + [len(balance) for balance in list_balance_insert[0]]
    
    plot_balance(balance_delete, balance_SNP, balance_insert, output_name, boundary, list_incidents, list_plot_name)

