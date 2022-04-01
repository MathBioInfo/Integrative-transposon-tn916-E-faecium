#!/usr/bin/env python
# coding: utf-8

# In[165]:


# to run this code we need some data
import pandas as pd
import numpy as np
from Bio import SeqIO, SeqRecord, Seq

import glob
from pathlib import Path
import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from numpy import sin
from numpy import sqrt
from scipy.optimize import curve_fit
dfL = pd.read_csv('processed_ICE_lengths.csv', sep=',') #read the processed lengths
#get the Tn916 names
ICEname = []
for i in range(len(dfL)):
    name = dfL['sseqid'][i].split('|')[-4] # example of 'sseqid' is ICEberg|329|Tn6079|GenBank|GU951538|462..28872 
    ICEname.append(name)
dfL['Icename']=ICEname # insert a new column into the dataframe with heading as "Icename" and entries from "ICEname"
dfL1= dfL[['isolate','Icename','length']] # define a new dataframe which consist of these three columns only
inpath = 'ICEDB/' #define the path
records = [] # read all the txt files from BLAST and look for ORF-1 .. 24. also record them
orf_map = {'ORF-{0}'.format(i+1): 0 for i in range(24)}
for file in Path(inpath).glob("*.txt"):
    # Get the empty ORF counter
    result_dict = orf_map.copy()
    
    # the variable 'file' will be the full path to a fasta file
    # We will save just the file name as a variable called "stem"
    stem = str(file).replace(inpath, '').replace('.txt', '').split('_')[-3]
    df = pd.read_csv(file, sep='\t', header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    
    # this has the ORFS
    v = df['sseqid'].tolist()
    v1 = df['pident'].tolist()
    zip_iterator = zip(v, v1)
    dicti = dict(zip_iterator)
    for orf in v:
        if dicti[orf]> 70:
            result_dict[orf] = 1
        
    v2 = df['qseqid'].unique()[0]
    
    result_dict['isolate'] = stem
    result_dict['Icename'] = str(v2)    
    #df1 = pd.DataFrame (v, columns = ['ORFs'])
    #df2=df1.drop_duplicates()
    records.append(result_dict)
r = pd.DataFrame.from_records(records)
#reorder the columns.
columns = ['isolate',"Icename"] + ["ORF-{}".format(i+1) for i in range(24)]
r = r[columns] # r now is a dataframe consisiting of teh record of all ORFs
joined = r.merge(dfL1, on=['isolate','Icename']) # join the two data frames by 'isolate' and 'Icename'
xis = ['ORF-1']
Int = ['ORF-2', 'ORF-2', 'ORF-3']
reco_cols = ['ORF-1','ORF-2','ORF-3','ORF-4']
regu_cols = ['ORF-5','ORF-6','ORF-7','ORF-8','ORF-9','ORF-10','ORF-12']
assc_cols = ['ORF-11']
cong_cols = ['ORF-13','ORF-14','ORF-15','ORF-16','ORF-17','ORF-18','ORF-19','ORF-20','ORF-21','ORF-22','ORF-23','ORF-24']
xissum = joined[xis].sum(axis = 1)
intsum = joined[Int].sum(axis = 1)
congsum = joined[cong_cols].sum(axis = 1)
recosum = joined[reco_cols].sum(axis = 1)
regsum = joined[regu_cols].sum(axis = 1)
accsum = joined[assc_cols].sum(axis=1)
orf = ['ORF-1', 'ORF-2', 'ORF-3', 'ORF-4', 'ORF-5', 'ORF-6', 'ORF-7', 'ORF-8', 'ORF-9', 'ORF-10', 'ORF-11', 'ORF-12', 'ORF-13', 'ORF-14', 'ORF-15', 'ORF-16', 'ORF-17', 'ORF-18', 'ORF-19', 'ORF-20', 'ORF-21', 'ORF-22', 'ORF-23', 'ORF-24']
allsum = joined[orf].sum(axis = 1)
joined['sumorf'] = allsum
# drop elements where length is greater than 5000 and number of ORfs are less than < 3.
# this will result in 2023 Tn916 elements in our data
joined.drop(joined[(joined['length'] > 5000) & (joined['sumorf'] < 4)].index, inplace = True)

# reindef the DataFrame
joined = joined.reset_index(drop=True)
# plot of length vs number of orfs
fig = plt.figure(figsize = (18, 5))
plt.plot(joined['length'], joined['sumorf'], 'o')
plt.legend(['Tn916 lengths vs number of ORFs'])
plt.xlabel('lenght', fontsize=16)
plt.ylabel('Number of ORFs', fontsize=16)
plt.grid()
plt.show()

plt.figure(figsize=(18,5))
plt.hist(joined['length'], bins = 50, color = 'm')
plt.legend(['The distribution of Tn916'])
plt.xlabel('length', fontsize=16)
plt.ylabel('Numbers', fontsize=16)
plt.grid()
plt.show()

sumorf = joined[orf].sum(axis = 0)
orf_count = dict(zip(orf, sumorf)) #create a dictionary of ORFs
ORF = list(orf_count.keys())
COUNT = list(orf_count.values())
fig = plt.figure(figsize = (18, 5))
fig.suptitle("ORFs Count", fontsize=16)
fig.text(0.5, 0.04, 'ORFs', ha='center', fontsize=16)
fig.text(0.09, 0.5, 'Number', va='center', rotation='vertical', fontsize=16)
plt.bar(ORF, COUNT, color ='maroon', width = 0.5)
plt.show()

# The enrichment of different clas of genes
Recombination = recosum/allsum
Excision=xissum/allsum
Integration = intsum/allsum
Conjugation = congsum/allsum
Regulation = regsum/allsum
tetM = accsum/allsum
# define the true objective function
def objective(x, a, b, c,  f):
	return (a * x) + (b * x**2) + (c * x**3)  + f

fig = plt.figure(figsize = (18, 5))
fig.suptitle("The enrichment of different class of genes", fontsize=16)
fig.text(0.5, 0.04, 'ORFs', ha='center', fontsize=16)
fig.text(0.09, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=16)

plt.subplot(2,3,1)

x, ycon = allsum, Conjugation
# curve fit
popt, _ = curve_fit(objective, x, ycon)
# summarize the parameter values
a, b, c,   f = popt
# plot input vs output
plt.scatter(x, ycon, color = 'g')
# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(min(x), max(x)+1, 1)
# calculate the output for the range
y_con = objective(x_line, a, b, c, f)
# create a line plot for the mapping function
plt.plot(x_line, y_con, '-', color='blue', linewidth=2)
plt.legend(['Best Fit', 'Conjugation Data'])
plt.grid()

plt.subplot(2,3,2)
x, yRec = allsum, Recombination
popt, _ = curve_fit(objective, x, yRec)
a, b, c,  f = popt
plt.scatter(x, yRec, color = 'r')
x_line = np.arange(min(x), max(x)+1, 1)
y_Rec = objective(x_line, a, b, c, f)
plt.plot(x_line, y_Rec, '-', color='blue', linewidth=2)
plt.legend(['Best Fit','Recombination Data'])
plt.grid()

plt.subplot(2,3,3)
x, yExc = allsum, Excision
popt, _ = curve_fit(objective, x, yExc)
a, b, c, f = popt
plt.scatter(x, yExc, color = 'c')
x_line = np.arange(min(x), max(x)+1, 1)
y_Exc = objective(x_line, a, b, c, f)
plt.plot(x_line, y_Exc, '-', color='blue', linewidth=2)
plt.legend(['Best Fit','Excision Data'])
plt.grid()

plt.subplot(2,3,4)
x, yReg = allsum, Regulation
popt, _ = curve_fit(objective, x, yReg)
a, b, c, f = popt
plt.scatter(x, yReg, color = 'm')
x_line = np.arange(min(x), max(x)+1, 1)
y_Reg = objective(x_line, a, b, c,f)
plt.plot(x_line, y_Reg, '-', color='blue', linewidth=2)
plt.legend(['Best Fit','Regulation Data'])
plt.grid()

plt.subplot(2,3,5)
x, ytetM = allsum, tetM
popt, _ = curve_fit(objective, x, ytetM)
a, b, c, f = popt
plt.scatter(x, ytetM, color = 'y')
x_line = np.arange(min(x), max(x)+1, 1)
y_tetM = objective(x_line, a, b, c, f)
plt.plot(x_line, y_tetM, '-', color='blue', linewidth=2)
plt.legend(['Best Fit','Accessory Data'])
plt.grid

plt.subplot(2,3,6)
x, ytetM = allsum, Integration
popt, _ = curve_fit(objective, x, ytetM)
a, b, c, f = popt
plt.scatter(x, ytetM, color = 'k')
x_line = np.arange(min(x), max(x)+1, 1)
y_tetM = objective(x_line, a, b, c, f)
plt.plot(x_line, y_tetM, '-', color='blue', linewidth=2)
plt.legend(['Best Fit','Integration Data'])
plt.grid
plt.show()
status=[None]*len(joined)
for i in range(len(joined)):
    if congsum.iloc[i] > 5 and xissum.iloc[i] > 0 and intsum.iloc[i] >1:
        status[i]='functional'
    else:
        status[i]='nonfunctional'
joined['Status'] = status # inser a new column 'status'
#joined.to_csv('pre_abs_ORFs.csv', index=False) #save as a csv
Functional = joined[joined["Status"] == 'functional']
#Functional.to_csv('functional_25_tn916.csv', index=False) # save to csv
Non_Functional = joined[joined["Status"] == 'nonfunctional']
#Non_Functional.to_csv('nonfunctional_25_tn916.csv', index=False) #save to csv


# In[ ]:




