# written to reformat the outputs from the drug interaction
# summary table to look at EHR examples of drug interactions
# written 12-11-19, JLW

import pickle, os, csv
import pandas as pd

f = '../data/cotherapy/supp4_summary_drug_interactions.xlsx'

dfs = pd.read_excel(f, sheetname=None)
sumdf = pd.concat([df.assign(dme=n) for n,df in dfs.items()])
new_cols = ['Predicted Cotherapy','Relationship','dme','Drugs with labeled DME','DME pathway intermediate','Intermediate node count','PMID','Notes']

sumdf = sumdf[new_cols]

sumdf.to_csv(open('../data/cotherapy/summary_data_for_EHR.csv','w'),index=False)

toi = ['Activates','Aggravates','Induces'] # interaction types of interest

actdf = sumdf.loc[sumdf['Relationship'].isin(toi)]
#actdf.explode('Drugs with labeled DME')
#
#pd.DataFrame([
#    [c, s, e, n]
#    for c, S, E, n in df.itertuples(index=False)
#    for s, e in zip(S.split(','), E.split(','))
#], columns=df.columns)

split_data = []
for (c,r,d,dlist,pgenes,p,n,ext) in actdf.itertuples(index=False):
	d_split = list(set([d for a in dlist.split('|') for d in a.split(',')]))
	for dname in d_split:
		split_data.append([c,r,d,dname,pgenes,p,n,ext])
actdf_split = pd.DataFrame(split_data, columns = actdf.columns)
actdf_split.to_csv(open('../data/cotherapy/activate_summary_data_for_EHR.csv','w'),index=False)

