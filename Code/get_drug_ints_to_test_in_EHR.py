# written to reformat the outputs from the drug interaction
# summary table to look at EHR examples of drug interactions
# written 12-11-19, JLW

import pickle, os, csv
import pandas as pd

# load data
f = '../data/cotherapy/supp4_summary_drug_interactions.xlsx'
dfs = pd.read_excel(f, sheetname=None)

# merge all sheets into a single sheet
sumdf = pd.concat([df.assign(dme=n) for n,df in dfs.items()])
new_cols = ['Predicted Cotherapy','Relationship','dme','Drugs with labeled DME','DME pathway intermediate','Intermediate node count','PMID','Notes']
sumdf = sumdf[new_cols]
sumdf.to_csv(open('../data/cotherapy/summary_data_for_EHR.csv','w'),index=False)

# specifically select for DME combinations that would exacerbate a DME
toi = ['Activates','Aggravates','Induces'] # interaction types of interest
actdf = sumdf.loc[sumdf['Relationship'].isin(toi)]
split_data = []
for (c,r,d,dlist,pgenes,p,n,ext) in actdf.itertuples(index=False):
	d_split = list(set([d for a in dlist.split('|') for d in a.split(',')]))
	for dname in d_split:
		split_data.append([c,r,d,dname,pgenes,p,n,ext])
actdf_split = pd.DataFrame(split_data, columns = actdf.columns)
actdf_split.to_csv(open('../data/cotherapy/activate_summary_data_for_EHR.csv','w'),index=False)

# repeat process to look for combinations that prevent a DME
toi = ['Prevcents','Prevents','Prevents?','Inhibits']
predf = sumdf.loc[sumdf['Relationship'].isin(toi)]
split_data = []
for (c,r,d,dlist,pgenes,p,n,ext) in predf.itertuples(index=False):
	d_split = list(set([d for a in dlist.split('|') for d in a.split(',')]))
	for dname in d_split:
		split_data.append([c,r,d,dname,pgenes,p,n,ext])
predf_split = pd.DataFrame(split_data, columns = predf.columns)
predf_split.to_csv(open('../data/cotherapy/prevent_summary_data_for_EHR.csv','w'),index=False) 

