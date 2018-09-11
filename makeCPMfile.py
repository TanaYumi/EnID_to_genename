#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:36:11 2018

@author: tanakayumiko
"""

# In[43]:


import pandas as pd
import numpy as np


# In[44]:
filepath = "directory/of/htseq/output"

#%%
#concat 15 counter files
Data = pd.DataFrame()
for sample in range(1,12):
    
    for il in range(1,12):
        #remove mushikui
            print('concat : sample%d_IL%d'%(sample,il))
            #read files
            filename = "%s/data/HTSeq/countdata-IL%02d-%d.txt"%(filepath,il,sample)
            data = pd.read_table(filename, sep='\t', header=None)
            #set columns names
            data.columns = ['Ensembl', '%02d-%d'%(il,sample)]
            #set index
            data = data.set_index('Ensembl')
            #concat
            Data = pd.concat([Data,data], axis=1)
            print('IL%02d-%d : %d indices'%(il,sample,len(data)))
Data = Data[:-5]
print(Data.head())
print('Data : %d genes'%len(Data.index.drop_duplicates()))


# In[45]:


#read Ensembl ID and gene name
EnID = pd.read_csv("%s/data/mart_export_human.csv"%filepath)
EnID = EnID.set_index('Gene stable ID')
print('raw EnID: len %d'%len(EnID))


# In[46]:


#remove doublet
print('Gene name duplicated : %d'%len(EnID[EnID.duplicated()]))
print('EnID duplicated : %d'%len(EnID[EnID.index.duplicated()]))
EnID = EnID.ix[EnID.index.drop_duplicates()]
print('EnID len : %d'%len(EnID))


# In[47]:


#ser EnID to DataFrame
EnData = pd.concat([EnID, Data], axis=1)
print('EnID and Genename, EnData\n', EnData.head())


# In[48]:


#check duplicated 'Gene name'
print('duplicated Gene name')
print(EnData[EnData['Gene name'].duplicated()].sort_values(by='Gene name').head())


# In[49]:


#drop NaN data
print('before drop Nan : %d genes'%len(EnData))
print('Nan genes : %d genes'%EnData.isnull().any(axis=1).sum())
EnData_dropna = EnData.dropna()
print('drop Nan : %d genes'%len(EnData_dropna))


# In[50]:


#save data
savefile = '%s/work/EnID_GeneName_count.csv'%filepath
EnData_dropna.sort_values(by='Gene name').to_csv(savefile)


# In[51]:


#EnData_dropnaの中で、Gene nameが重複しているものを表示
print('duplicated EnData_dropna : %d genes'%len(EnData_dropna[EnData_dropna.duplicated('Gene name')]))
print(EnData_dropna[EnData_dropna.duplicated('Gene name')].head())


# In[52]:


#Gene name をindexに指定する
EnData_g = EnData_dropna.set_index('Gene name')

print('EnData_g : %d genes'%len(EnData_g))
print(EnData_g.head())
EnData_g.ix['CD3E']


# In[53]:


#Gene nameの重複した行のみ取り出す
duplicated_GeneName_index = EnData_g[EnData_g.index.duplicated()].index.drop_duplicates()
print('duplicated EnID : %d genes'%len(EnData_g.ix[duplicated_GeneName_index]))
print('duplicated Gene name : %d genes'%len(duplicated_GeneName_index))
print(EnData_g.ix[duplicated_GeneName_index].head())
w=EnData_g.ix[duplicated_GeneName_index]
w.ix['Metazoa_SRP']


# In[54]:


#Gene nameの重複した行どうしで和を取る
#まずcolumnsに重複した遺伝子を追加していき、後で転換する
dub = pd.DataFrame(index = EnData_g.columns)
for i in duplicated_GeneName_index.drop_duplicates():
    dub['%s'%i] = EnData_g.ix[i].sum(axis=0)
dub_T = dub.T.sort_index()
print('index is unique : ', dub_T.index.is_unique)
print('duplicated genes dub_T : %d genes'%len(dub_T))
print(dub_T.head())
dub_T.ix['Metazoa_SRP']


# In[55]:


print('EnData_g length : %d genes'%len(EnData_g.index))
print('drop duplicates from EnData_g.index : %d genes'%len(EnData_g.index.drop_duplicates()))


# In[56]:


#重複遺伝子を除いたDataFrameをまず作成し、あとから重複遺伝子の和のDataFrameを結合する
print(EnData_g.ix['CD3E'])
unique_index = EnData_g.index.drop_duplicates()

unique = EnData_g.ix[unique_index]
unique = unique.drop(dub_T.index)
print(unique.ix['CD3E'])
Data = pd.concat([unique,dub_T])
print('index of Data is unique : ', Data.index.is_unique)
print('concat Data : %d'%len(Data))
print(Data)
print(Data.ix['CD3E'])
print(Data.ix['Metazoa_SRP'])


# In[57]:


#Nan, duplicates除去データ保存
Data.to_csv('%s/work/dn_dd_count.csv'%filepath)
#Nan, duplicates, allzero除去データ保存
az_Data = Data[Data.sum(axis=1)!=0]
print(az_Data.head())
print('remove all_zero genes : %d'%len(az_Data))
az_Data.to_csv('%s/work/dn_dd_az_count.csv'%filepath)
az_Data


