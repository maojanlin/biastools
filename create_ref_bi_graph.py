import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

#============================= STEP 1 =================================
#print(plt.style.available)
TITLE = 'Reference Bias Distribution' # of NA12878 Reads \n Aligned to GRCh37 Using Bowtie2'
FN_FIG = 'reference_bias_graphs.pdf'#'SRR622457-grch37-bt2.pdf'
FN_INPUT = './bt.bias'

df = pd.read_csv(FN_INPUT, sep='\t')
df_mid = df[df['ALT_COUNT'] + df['REF_COUNT']> 0]
#df_almost = df_mid[df_mid['REF COUNT'] > 0]
df_use = df_mid[df_mid['NUM_READS'] >= 15]
df_use.loc[:,'THRESHOLD'] = 15
df_use.head()


#============================= STEP 2 =================================
rb = df_use['REFERENCE_BIAS']
if 'NaN' in rb:
    print("true")
print(rb)
#: show stats
print (rb.quantile(q=[.01, .1, 0.25,0.5,0.75, .90, .99]))
print("mean: ", np.mean(rb))

plt.clf()
n, bins, patches = plt.hist(rb, bins = list(np.linspace(0,1,26)))
plt.ylim((1,40000))
plt.xlim((0,1))
plt.yscale('log')
plt.xlabel('Reference Bias')
plt.ylabel('Counts')
plt.title(TITLE)
plt.savefig(FN_FIG)


#============================= STEP 3 =================================
FN_FIG = 'reference_bias_eighty.pdf'#'SRR622457-grch37-bt2.pdf'
df_eighty = df_use[df_use['REFERENCE_BIAS'] >= 0.80]
df_eighty_refbi = df_eighty['REFERENCE_BIAS']
#print(len(df_eighty_refbi))

title = 'Reference Bias Distribution of HET sites with Reference Bias >= 0.80'
plt.clf()
plt.hist(df_eighty_refbi)
plt.yscale('log')
plt.xlabel('Reference Bias')
plt.ylabel('Counts')
plt.title(title)
plt.savefig(FN_FIG)



#============================= STEP 4 =================================
FN_FIG = 'reference_bias_dot.pdf'#'SRR622457-grch37-bt2.pdf'
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
sp = pd.DataFrame()

# Add columns
df_eighty['AVERAGE_MAPQ'] = df_eighty['SUM_MAPQ']/(df_eighty['NUM_READS'])
mapQ = list(df_eighty['AVERAGE_MAPQ'])

sp['REFERENCE BIAS'] = list(df_eighty_refbi)
sp['READ DISTRIBUTION'] = list(df_eighty['READ_DISTRIBUTION'])
sp.head()

mapped_mapQ = [map_mapq_to_size(q) for q in mapQ]
sp['Ambiguity'] = mapped_mapQ
sp['MapQ'] = list(mapQ)
#sp['# Variants'] = list(size_in['# Variants']+1)

plt.clf()
#ax = sns.scatterplot(x="REFERENCE BIAS", y="Read Distribution",  hue = "Ambiguity", size = "# Variants", data = sp)#hue="size", size="size", data=tips)
ax = sns.scatterplot(x="REFERENCE BIAS", y="READ DISTRIBUTION",  hue = "Ambiguity", data = sp)#hue="size", size="size", data=tips)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(FN_FIG)





# In[5]:

'''
print (min(rb))
print (max(rb))

import math 
count = 0
for i, r in enumerate(rb):
    if math.isnan(r):
        print("count: ", count)
        print (df_use.iloc[i])
    count += 1
print(count)


# In[6]:


bplot = sns.boxplot(y='REFERENCE BIAS', x='THRESHOLD', 
                 data=df_use, 
                 width=0.7,
                 palette="colorblind").set_title("Allelic Bias of SRR622457 reads aligned to GRCh37 using Bowtie2 With At least 15x Coverage")

sns.set(style="whitegrid")


# In[7]:


df_ninety = df_use[df_use['REFERENCE BIAS'] >= 0.90]
#print(df_ninety)
df_ninety_refbi = df_ninety['REFERENCE BIAS']
#print(df_ninety_refbi)

title = 'Reference Bias Distribution of HET sites with Reference Bias >= 0.90'
plt.clf()
plt.hist(df_ninety_refbi)
plt.yscale('log')
plt.xlabel('Reference Bias')
plt.ylabel('Counts')
plt.title(title)


# In[11]:


FN_INPUT = '/Users/sheilaiyer/Documents/JHU/text_files/bt2_extra_bias.txt'#'/net/langmead-bigmem-ib.bluecrab.cluster/storage/sheila/faster_complete-full-bt2-refbias.txt'
bias = pd.read_csv(FN_INPUT, sep='\t')

gb = bias['ALT BIAS']
ob = bias['OTHER BIAS']


# In[14]:


import matplotlib.pyplot as plt

title = 'Gap Bias for Bowtie2 Real Data'
plt.clf()
plt.hist(gb)
plt.yscale('log')
plt.xlabel('Gap Bias')
plt.ylabel('Counts')
plt.title(title)


# In[15]:


title = 'Other Bias for Bowtie2 Real Data' 
plt.clf()
plt.hist(ob)
plt.yscale('log')
plt.xlabel('Other Bias')
plt.ylabel('Counts')
plt.title(title)


# In[20]:


print ('#gap>0.2 (%)', sum(gb>0.2)/float(len(gb))*100)
print ('avg #gap', sum(gb)/float(len(gb)))
print ('#other>0.2 (%)', sum(ob>0.2)/float(len(ob))*100)
print ('avg #other', sum(ob)/float(len(ob)))


# In[ ]:
'''



