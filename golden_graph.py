import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import math
import random
import numpy as np
import pysam

#============================= STEP 0 DATA INPUT =================================
#FN_INPUT = 'merge.tribus.bt.bias.gap'
#f_exclusion = "GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed"
#f_e = open(f_exclusion, 'r')
#dict_exclude = {}
#for line in f_e:
#    fields = line.strip().split()
#    if dict_exclude.get(fields[0]):
#        dict_exclude[fields[0]] += list(range(int(fields[1]), int(fields[2])))
#    else:
#        dict_exclude[fields[0]] = list(range(int(fields[1]), int(fields[2])))

FN_INPUT = 'merge.compact.bt.bias.SNP'
FN_INPUT = 'merge.outrim.all.bt.bias.SNP'
FN_INPUT = 'merge.chrX.bias.SNP'
FN_INPUT = 'merge.chr21.naive.bias.SNP'
FN_INPUT = 'merge.chm13_grch38.chr21.bias.SNP'
FN_INPUT = 'merge.chr21.bias.SNP'
FN_INPUT = 'merge.refflow-10.chr21.bias.SNP'
FN_INPUT = 'merge.major.chr21.bias.SNP'

df_use = pd.read_csv(FN_INPUT, sep='\t')
#df = pd.read_csv(FN_INPUT, sep='\t')
#df_mid = df[df['ALT_COUNT'] + df['REF_COUNT']> 0]
#df_use = df_use[df_use['NUM_READS'] >= 10]
####df_use = df_use[df_use['HET_SITE'].isin(dict_exclude['chr21'])]
#df_use.loc[:,'THRESHOLD'] = 15
df_use.head()

rb = df_use['REFERENCE_BIAS']
if 'NaN' in rb:
    print("true")
print(rb)


#============================= STEP 3 =================================
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


# Add columns
df_use['AVERAGE_MAPQ'] = df_use['SUM_MAPQ']/(df_use['NUM_READS']) # the average map_q score
df_use['WASTE_INFO']   = (df_use['BOTH_COUNT']+df_use['NEITHER_COUNT'])/df_use['NUM_READS']
mapQ = list(df_use['AVERAGE_MAPQ'])
pValue = list(df_use['EVEN_P_VALUE'])

sp = pd.DataFrame()
sp['VARIANT BIAS']        = list(rb)
sp['READ DISTRIBUTION']   = list(df_use['READ_DISTRIBUTION'])
sp['GOLDEN DISTRIBUTION'] = list(df_use['GOLDEN_DISTRIBUTION'])
sp.head()

mapped_mapQ = [map_mapq_to_size(q) for q in mapQ]
mapped_p    = [map_color(p) for p in pValue]
waste_value = [map_waste_to_color(q) for q in list(df_use['WASTE_INFO'])]
sp['Avg_MapQ_code'] = mapped_mapQ
sp['Even_p_value']  = mapped_p
sp['Waste_value']   = waste_value
sp['MapQ'] = list(mapQ)

#=========================== standard ref_bias to read_distribute plot ============================
plt.clf()
ax = sns.scatterplot(y="VARIANT BIAS", x="READ DISTRIBUTION",  hue = "Avg_MapQ_code", data = sp)#hue="size", size="size", data=tips)
#ax = sns.scatterplot(y="VARIANT BIAS", x="READ DISTRIBUTION",  hue = "Waste_value", data = sp)#hue="size", size="size", data=tips)
h, l = ax.get_legend_handles_labels()
plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., framealpha=1)
plt.xlim([0,1])
plt.ylim([0,1])

FN_FIG = FN_INPUT + '-var2read_dot.pdf'
plt.savefig(FN_FIG)
print("Standard Ref Bias Plot!")

#=========================== golden to read_distribute plot ============================
plt.clf()
ax = sns.scatterplot(x="GOLDEN DISTRIBUTION", y="READ DISTRIBUTION",  hue = "Avg_MapQ_code", data = sp)#hue="size", size="size", data=tips)
h, l = ax.get_legend_handles_labels()
plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., framealpha=1)
plt.xlim([0,1])
plt.ylim([0,1])

FN_FIG = FN_INPUT + '-gold2read_dot.pdf'
plt.savefig(FN_FIG)

#=========================== golden to ref_bias plot ============================
plt.clf()
ax = sns.scatterplot(x="GOLDEN DISTRIBUTION", y="VARIANT BIAS",  hue = "Avg_MapQ_code", data = sp)#hue="size", size="size", data=tips)
h, l = ax.get_legend_handles_labels()
plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., framealpha=1)
plt.xlim([0,1])
plt.ylim([0,1])

FN_FIG = FN_INPUT + '-gold2var_dot.pdf'
plt.savefig(FN_FIG)

#=========================== all merged plot ============================
plt.clf()
sp['Normalized Variant Bias']      = list(df_use['REFERENCE_BIAS']-df_use['GOLDEN_DISTRIBUTION']) # the average map_q score
sp['Normalized Read Distribution'] = list(df_use['READ_DISTRIBUTION']-df_use['GOLDEN_DISTRIBUTION']) # the average map_q score
ax = sns.jointplot(x="Normalized Read Distribution", y="Normalized Variant Bias",  hue = "Avg_MapQ_code", data = sp, xlim=(-0.6,0.6), ylim=(-0.6,0.6))
ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
ax.ax_joint.get_legend().remove()
h, l = ax.ax_joint.get_legend_handles_labels()
plt.legend(h, labels, title="Avg MapQ", bbox_to_anchor=(0, 0), loc='lower right', borderaxespad=0.2)

FN_FIG = FN_INPUT + '-merged_golden_dot.pdf'
plt.savefig(FN_FIG)

#======================= allelic difference plot =========================
plt.clf()
list_ref_diff = list(df_use['REF_COUNT']-df_use['GOLDEN_REF_COUNT'])
list_alt_diff = list(df_use['ALT_COUNT']-df_use['GOLDEN_ALT_COUNT'])
for idx in range(len(list_ref_diff)):
    list_ref_diff[idx] += random.uniform(-0.3, 0.3) # scatter plot
    list_alt_diff[idx] += random.uniform(-0.3, 0.3)
sp['Ref read diff'] = list_ref_diff
sp['Alt read diff'] = list_alt_diff

ax = sns.jointplot(x="Ref read diff", y="Alt read diff",  hue = "Even_p_value", data = sp, xlim=(-20,20), ylim=(-20,15))
ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
ax.ax_joint.get_legend().remove()
h, l = ax.ax_joint.get_legend_handles_labels()
plt.legend(h, p_labels, title="Even P Value", bbox_to_anchor=(0,1), loc='upper right')

FN_FIG = FN_INPUT + '-read_diff_allelic.pdf'
plt.savefig(FN_FIG)

#====================== mapping difference plot =========================
plt.clf()
list_m_ref_diff = list(df_use['NUM_READS']*df_use['READ_DISTRIBUTION']-df_use['GOLDEN_REF_COUNT'])
list_m_alt_diff = list(df_use['NUM_READS']*(1-df_use['READ_DISTRIBUTION'])-df_use['GOLDEN_ALT_COUNT'])
for idx in range(len(list_m_ref_diff)):
    list_m_ref_diff[idx] += random.uniform(-0.3, 0.3) # scatter plot
    list_m_alt_diff[idx] += random.uniform(-0.3, 0.3)
sp['Ref read diff'] = list_m_ref_diff
sp['Alt read diff'] = list_m_alt_diff

ax = sns.jointplot(x="Ref read diff", y="Alt read diff",  hue = "Even_p_value", data = sp, xlim=(-20,20), ylim=(-20,15))
ax.ax_joint.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.1)
ax.ax_joint.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.1)
ax.ax_joint.get_legend().remove()
h, l = ax.ax_joint.get_legend_handles_labels()
plt.legend(h, p_labels, title="Even P Value", bbox_to_anchor=(0,1), loc='upper right')

FN_FIG = FN_INPUT + '-read_diff_mapping.pdf'
plt.savefig(FN_FIG)



#=========================== all merged plot ============================
"""
plt.clf()
sp['Variant - Read Bias'] = list(df_use['REFERENCE_BIAS']-df_use['READ_DISTRIBUTION']) # the average map_q score
sp['Read - Golden Bias']  = list(df_use['READ_DISTRIBUTION']-df_use['GOLDEN_DISTRIBUTION']) # the average map_q score
ax = sns.scatterplot(x="Read - Golden Bias", y="Variant - Read Bias",  hue = "Avg_MapQ_code", data = sp)#hue="size", size="size", data=tips)
#ax = sns.scatterplot(x="Read - Golden Bias", y="Variant - Read Bias",  hue = "Waste_value", data = sp)#hue="size", size="size", data=tips)
ax.axhline(y=0, color='gray', linestyle='dashdot', linewidth=0.2)
ax.axvline(x=0, color='gray', linestyle='dashdot', linewidth=0.2)
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
plt.xlim([-0.6,0.6])
plt.ylim([-0.6,0.6])

FN_FIG = FN_INPUT + '-merged_read_dot.pdf'
plt.savefig(FN_FIG)
"""



