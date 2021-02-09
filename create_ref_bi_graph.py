#!/usr/bin/env python
# coding: utf-8

# In[9]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#styling: https://matplotlib.org/users/style_sheets.html
print(plt.style.available)
#FN_INPUT = '/Users/sheilaiyer/Documents/JHU/text_files/SRR622457_chr21_bt2_python.txt'#'/net/langmead-bigmem-ib.bluecrab.cluster/storage/sheila/faster_complete-full-bt2-refbias.txt'
TITLE = 'Reference Bias Distribution' # of NA12878 Reads \n Aligned to GRCh37 Using Bowtie2'
FN_FIG = 'reference_bias_graphs.pdf'#'SRR622457-grch37-bt2.pdf'
FN_INPUT = '/Users/sheilaiyer/resources/biastools/hapA_ref_bi.txt'


df = pd.read_csv(FN_INPUT, sep='\t')
df_mid = df[df['ALT_COUNT'] + df['REF_COUNT']> 0]
#df_almost = df_mid[df_mid['REF COUNT'] > 0]
df_use = df_mid[df_mid['NUM_READS'] >= 15]
df_use.loc[:,'THRESHOLD'] = 15
df_use.head()


# In[10]:


import numpy as np

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



