import sys
import os
import pandas as pd
from pathlib import Path
result = []
root = "../ICEtest/"

for file in Path(root).glob("*.txt"):
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    # split(char) splits a string into a list using char as a delimiter
    # list[-1] returns the last entry in a list
    stem = str(file).replace('.txt', '').split('/')[-1]
    rslt_df = df.loc[(df['length']>100)]
    dataff = rslt_df.sort_values(by='length',ascending=False)
    dataff.drop_duplicates(subset ="qstart", keep = "first",inplace = True)
    dataff1 = dataff.sort_values(by='length',ascending=False)
    dataff1.drop_duplicates(subset ="qend", keep = "first",inplace = True)
    # creat a directiry
    #os. mkdir(stem)
    #store dataff1 to the created directroy
    #subroot = stem
    #subdir = os.path.join(root, subroot)
    #dataff1.to_csv(os.path.join(subdir, stem+".csv"), index=False)
    
    final = dataff1.groupby('sseqid')['length'].sum().sort_values(ascending=False).to_frame(name = 'length').reset_index()
    final1 = final[(final.length >5040)]
        
    
    final1['isolate'] = stem
    
    result.append(final1)

all_result = pd.concat(result)
all_result.to_csv('processed_ICE_lengths.csv', index=False)    
sys.exit(0)

