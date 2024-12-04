### modules
import numpy as np
import pandas as pd
from scipy import constants
import functions as iscore

### set drug
antag = ['CTLA4']   # multiple antagonists can be considered
# Nivolumab (anti-PD1)          'PDCD1'
# Ipilimumab (anti-CTLA4)       'CTLA4'
# Durvalmab (anti-PDL1)         'CD274'
# Relatlimab (anti-LAG3)        'LAG3'
# Tiragolumab (anti-TIGIT)      'TIGIT'

### read data
data1 = pd.read_csv("./dat/data_ExpressionMatrix_Pre.csv", index_col=0)
#data1 = pd.read_csv("./dat/data_ExpressionMatrix_During.csv", index_col=0)
#data1 = pd.read_csv("./dat/data_ExpressionMatrix_JustFinish.csv", index_col=0)

### reference data
data2 = pd.read_csv("./dat/data_proteomics_reference.csv", index_col=0)
df1 = pd.read_csv("./dat/data_Ligand_Receptor_KD.csv")

### cell types considered
l_celltype = ['B.memory', 'B.naive', 'B.plasma', 'Treg.memory', 'Treg.naive', 'T4.naive', 'T4.CM', 'T4.EM', 'T4.EMRA', 'Th1', 'Th17', 'Th2', 'T8.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'mDC', 'pDC', 'Mono.classical', 'Mono.intermediate', 'Mono.nonclassical', 'NK.bright', 'NK.dim', 'Basophil', 'Eosinophil', 'Neutrophil']
state1 = '_activated'
state2 = '_steady-state'

l1 = []
for st1 in l_celltype:
    if st1 in list(data1.columns):
        l1.append(st1)
#print(l1)

l2 = []
for i1, st1 in enumerate(l1):
    for st2 in l1[i1:]:
        #print(st1, st2)
        if (st1 == st2):
            if (st1+state1 not in data2.columns):
                l2.append([st1+state2, st2+state2, 1])
            else:
                l2.append([st1+state2, st2+state2, 0.25])
                l2.append([st1+state2, st2+state1, 0.25])
                l2.append([st1+state1, st2+state2, 0.25])
                l2.append([st1+state1, st2+state1, 0.25])
        elif (st1 != st2) and (st1+state1 not in data2.columns):
            if (st2+state1 not in data2.columns):
                l2.append([st1+state2, st2+state2, 1])
            else:
                l2.append([st1+state2, st2+state2, 0.5])
                l2.append([st1+state2, st2+state1, 0.5])
        elif (st1 != st2) and (st1+state1 in data2.columns):
            if (st2+state1 not in data2.columns):
                l2.append([st1+state2, st2+state2, 0.5])
                l2.append([st1+state1, st2+state2, 0.5])
            else:
                l2.append([st1+state2, st2+state2, 0.25])
                l2.append([st1+state2, st2+state1, 0.25])
                l2.append([st1+state1, st2+state2, 0.25])
                l2.append([st1+state1, st2+state1, 0.25])
#print(l2)

### calculate interaction score for each cell type pair
l3 = []
l4 = []
for row1 in l2:
    #print(row1)
    st1 = row1[0]
    st2 = row1[1]
    x1 = iscore.interaction_score(st1, st2, data1, data2, df1, antag) * row1[2]
    #print( x1 )
    l3.append([st1, st2])
    l4.append(x1)
df_out = pd.concat( [pd.DataFrame(l3), pd.DataFrame(l4)], axis=1 )
clm = ['celltype1', 'celltype2', 'x1', 'y1']
df_out.columns = clm
# x1: total score woDrug, y1: score for target woDrug

#print(df_out)
#df_out.to_csv('result_IS_x.csv')

df1 = pd.DataFrame(l4, columns=['x1', 'y1'])
df1['score2'] = df1['x1'] * df1['y1']
#print(df1['score2'])

### treatment score
print('Treatment Score with drug for', antag)
print(df1['score2'].sum())

