#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Extracting ICE sequences
import pandas as pd
import numpy as np
from Bio import SeqIO, SeqRecord, Seq

import glob
from pathlib import Path
df = pd.read_csv("processed_ICE_lengths.csv") # read the cvs file 
for i, row in df.iterrows():                  # get different columns from the file
    ICEname = str(row['sseqid']).split('|')[-4]
    num = str(i)
    Iso = str(row['isolate'])
    save_path = '/home/amjad/Desktop/Dsktop/test/amjad_presence_absence/ICEamrTEST' # path to savee file
    name_to_save = Iso+'_' + num+'.fasta'
    completeName = os.path.join(save_path, name_to_save)
    out_file = open(completeName, "w") # open file
    
    length_name = Iso+".csv" # read individual csv file, containing information for the individual ICE
    length_file = "/home/amjad/Desktop/Dsktop/test/amjad_presence_absence/BALST_results/"+length_name
    df1 = pd.read_csv(length_file)
    df1['Isolate'] = Iso
    for j, rowq in df1.iterrows():
        ICEcontig = str(j)
        contig = rowq['qseqid']
        ICEname1 = str(rowq['sseqid']).split('|')[-4]
        length = str(rowq['length'])
        #print(contig)
        #print(ICEname1)
        start = rowq['qstart']
        end = rowq['qend']
        name = Iso+".fasta"
        filename = "/home/amjad/NewFASTA_Assemblies_Efaecium/"+name
        
        for seq in SeqIO.parse(filename,"fasta"):
            i =1
            if (seq.id == str(contig)) & (ICEname == ICEname1):
                
                #out_file.write(">" + ICEname+"_"+str(contig)+ "_"+str(start)+"_"+str(end)+"\n")
                out_file.write(">"+"\n")
                
                out_file.write( str(seq.seq[start-1:end]) + "\n")
                
                
                
            
    #break
out_file.close()

