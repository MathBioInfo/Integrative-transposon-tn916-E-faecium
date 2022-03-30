#!/usr/bin/env python
# coding: utf-8

# In[58]:


import pandas as pd
import random
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import time
deltat= 0.1
rD = 0.2
rE = 1
rS = 0.1
rF = 0.1
rD = rD*deltat
rE = rE*deltat
rC = rC*deltat
rS = rS*deltat
rF = rF*deltat
Ubound = 5000 #upper bound on ICE population
Ntime = 100
timeiteration = Ntime/deltat
Re_enter = pd.DataFrame()
Congplasmid = pd.DataFrame()
Plasmid = pd.DataFrame()
allsum = pd.DataFrame()
Degraded = pd.DataFrame()
totDegraded = pd.DataFrame()
cong = pd.DataFrame()
reco = pd.DataFrame()
reg = pd.DataFrame()
asc = pd.DataFrame()
GENEMEAN = pd.DataFrame()
GenVSice = pd.DataFrame()
GenVexciable = pd.DataFrame() # to keep track of ice capable of conjugation
GenVnonxeciable = pd.DataFrame() # to kepp tyrack of ICE uncapable fo conjugation
GenVcongplasmid = pd.DataFrame() # to keep track of conjugative plasmids
rows, cols = (24, 1000)
ICE = [[1 for i in range(rows)] for j in range(cols)]
df = pd.DataFrame(ICE)
orf_map = {'ORF-{0}'.format(i+1): 1 for i in range(24)}
df1 = pd.DataFrame([orf_map])
head = df1.columns
df.columns=head
ICEcomplete = pd.DataFrame(np.array([1]*24))
ICEcomplete = ICEcomplete.transpose()
ICEcomplete.columns = head
xis = ['ORF-1']
Int = ['ORF-2', 'ORF-2', 'ORF-3']
reco_cols = ['ORF-1','ORF-2','ORF-3','ORF-4']
regu_cols = ['ORF-5','ORF-6','ORF-7','ORF-8','ORF-9','ORF-10','ORF-12']
assc_cols = ['ORF-11']
cong_cols = ['ORF-13','ORF-14','ORF-15','ORF-16','ORF-17','ORF-18','ORF-19','ORF-20','ORF-21','ORF-22','ORF-23','ORF-24']
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
totdegraded = 0 ###################
for t in tqdm(np.arange(Ntime/deltat)):# main time loop
    df = df.reset_index(drop=True)
    for i in range(len(df)-1):
        for j in range(rows):
            rand = random.uniform(0, 1)
            if rD>rand:
                df.iloc[i,j] = 0
    allsum = df.sum(axis = 1)
    xissum = df[xis].sum(axis = 1)
    intsum = df[Int].sum(axis = 1)
    congsum = df[cong_cols].sum(axis = 1)
    recosum = df[reco_cols].sum(axis = 1)
    regsum = df[regu_cols].sum(axis = 1)
    conjable = 0 # excisable Tn916 elements
    for i in range(len(df)-1):
        if xissum.iloc[i] > 0: # xis is important for the excission of Tn916 element
            conjable += 1 # excisable elements
    exable = pd.DataFrame([{'Gen':t, 'Number':conjable}])
    GenVexciable = GenVexciable.append(exable)  # excisable
    ###############################################################
    nonconjable = len(df)-conjable
    nonexciable = pd.DataFrame([{'Gen':t, 'Number':nonconjable}])
    GenVnonxeciable = GenVnonxeciable.append(nonexciable) 
    #################################################################
    if conjable == 0: # use this if we are using closed system (rF = 0)
        print('"There is no more excisable ICE in the population at generation ',t)
        break
    re_enter = [] # keep track of ICE that can renter into genome
    congplasmid =[] # keep track oof ICE not capable of entering into genome
    for i in range (len(congsum)):
        rand = random.uniform(0, 1)
        if congsum.iloc[i] > 5 and xissum.iloc[i] >0 and intsum.iloc[i] >1 and rand < rE: # if conjsuative and xis machinery is present
            re=df.iloc[i].copy()
            re_enter.append(re)
            #print(congsum.iloc[i], rand, rE, i)
            # if an ice is having conjugative machinary intact, xis present and Int absent then they may have life as plasmid
        elif congsum.iloc[i] > 5 and xissum.iloc[i] >0  and intsum.iloc[i] == 0 and rand < rE:
            plas = df.iloc[i].copy()
            congplasmid.append(plas)
            #print(recosum.iloc[i], rand, rE, i)
            df.iloc[i] = 0
    re_enter = pd.DataFrame(re_enter)
   
    
    df = df.append(re_enter) #append to the total population ICEs that were capable to reenter into the genome
    congplasmid = pd.DataFrame(congplasmid)
    allsum = df.sum(axis = 1)
    for i in range (len(allsum)):
        if allsum.iloc[i] ==0:
            df.drop([i], inplace=True)
            totdegraded += 1 # keep track of in sll generations
    ############################################
    for i in range(len(df)-1): # new ICEs can also enter into the population
        rand = random.uniform(0,1)
        if rand < rF:
            df = df.append(ICEcomplete) 
    ##############################################################
    Re_enter = Re_enter.append(re_enter) # total re_entered ICEs
    Congplasmid = Congplasmid.append(congplasmid) # total plasmid type? 
    totdegr = pd.DataFrame([{'Gen':t, 'Number':totdegraded}])
    totDegraded = totDegraded.append(totdegr)  # these total ICEs lost from the population not just due to degradation
    # degradation just to mutation
    degr = pd.DataFrame([{'Gen':t, 'Number':totdegraded-len(Congplasmid)}])
    Degraded = Degraded.append(degr)
    #################################################
    plasmid = pd.DataFrame([{'Gen':t, 'Number':len(Congplasmid)}])
    Plasmid = Plasmid.append(plasmid) # to plot Cojugative plasmids
    ####################################################################
    ############################################# Selection
    sel = []
    asc = df.iloc[:,10] # tetM gene
    for i in range (len(asc)):
        rand = random.uniform(0, 1)
        if rand < rS*asc.iloc[i]:
            newsel = df.iloc[i].copy()
            sel.append(newsel)
    sel = pd.DataFrame(sel)
    # Reproduction and population regulations, we are usingt Wright_Fishaer kind of model with slight modafication
    Rep = []
    Wright_Fisher =  pd.DataFrame()
    frac = 1 - len(df)/(10*Ubound + len(sel))
    for i in range(len(df)):
        rand = random.uniform(0, 1)
        if rand < frac:
            next_gen = df.iloc[i].copy()
            Rep.append(next_gen)
    Rep = pd.DataFrame(Rep)
    Wright_Fisher = sel.append(Rep)
#     print('Frac_lose is :', frac_lose, 'Randon is :', rand)
#     print('The length of origional df is :', len(df))
#     print('The length of new Rep is :', len(Rep))
#     print('the length of sel is:' , len(sel))
    ############################################# New pole of ICE for next generation
    df = Wright_Fisher
 #   print('The length of wf is :', len(Wright_Fisher))
    if len(Wright_Fisher) == 0:
        print('ICE population got extinct. Please restart the simulations')
        break    
#    step = t*deltat
#     if step%5 == 0:
#         gmean = df.mean(axis=0)
#         gmean = pd.DataFrame(gmean)
#         gmean = gmean.transpose()
#         gmean.insert(loc=0, column='Time', value=step)
#         GENEMEAN = GENEMEAN.append(gmean)
    gmean = df.mean(axis=0)
    gmean = pd.DataFrame(gmean)
    gmean = gmean.transpose()
    gmean.insert(loc=0, column='Time', value = t)
    GENEMEAN = GENEMEAN.append(gmean) 
        
    genvsice = pd.DataFrame([{"Gen": t, "Copies": len(df)}])
    GenVSice = GenVSice.append(genvsice)
#     print('Frac_lose is :', frac_lose, 'Randon is :', rand)
#     print('The length of origional df is :', len(df))
#     print('The length of new Rep is :', len(Rep))
#     print('the length of sel is:' , len(sel))
    ############################################# New pole of ICE for next generation
    df = Wright_Fisher
#    print('The length of wf is :', len(Wright_Fisher))
    #break
    time.sleep(0.1)
GENEMEAN.to_csv(r'genemean.csv',index = False)
df.to_csv(r'ICEPresAbsTry.csv', index = False)
GenVSice.to_csv(r'GenVSice.csv', index = False)
GenVexciable.to_csv(r'GenVexciable.csv',index = False) # ICEs capale of conjugation
Degraded.to_csv(r'Degraded.csv',index = False) # Degraded ICE
summean = pd.DataFrame()
summean['reco'] = GENEMEAN[reco_cols].sum(axis=1)/4
summean['regu'] = GENEMEAN[regu_cols].sum(axis=1)/7
summean['assc'] = GENEMEAN[assc_cols]/1
summean['cong'] = GENEMEAN[cong_cols].sum(axis=1)/12
plt.figure(figsize=(10, 5),dpi =100)
plt.subplot(3,1,1)
plt.plot(GENEMEAN['Time'], summean['reco'], color = 'blue')
plt.plot(GENEMEAN['Time'], summean['regu'], color = 'green')
plt.plot(GENEMEAN['Time'], summean['assc'], color = 'red')
plt.plot(GENEMEAN['Time'], summean['cong'], color = 'yellow')
plt.legend(['Recombination', 'Regulation', 'Accessory',  'Conjugation'])
#plt.xlabel('Time')
plt.ylabel('Ave: # of genes')
plt.title('Pinit ='+str(cols)+',rD ='+str(rD)+',rE = '+str(rE)+',rS='+str(rS)+',rF='+str(rF))
plt.grid()
plt.subplot(3, 1, 2)
plt.plot(GenVSice['Gen'], GenVSice['Copies'], color = 'blue')
plt.plot(GenVexciable['Gen'], GenVexciable['Number'], linestyle = 'dotted', color = 'red')
plt.plot(GenVnonxeciable['Gen'], GenVnonxeciable['Number'], linestyle = 'dotted',  color = 'green')
plt.legend(['Total # of ICE', '# of excisable', '# of non_excisable'])
#plt.xlabel('Generation')
plt.ylabel('Number')
plt.grid()
plt.subplot(3,1,3)
plt.plot(totDegraded['Gen'], totDegraded['Number'], color = 'red')
plt.plot(Degraded['Gen'], Degraded['Number'],  linestyle='-.', color = 'blue')
plt.plot(Plasmid['Gen'], Plasmid['Number'],  linestyle='-.', color = 'green')
plt.legend(['# of ICE lost from the pop:', '# of degraded ICE','# conjugative plasmids'])
plt.xlabel('Generation')
plt.ylabel('Number')
plt.grid()
#plt.savefig('figure.pdf')
plt.show()
    


# In[ ]:




