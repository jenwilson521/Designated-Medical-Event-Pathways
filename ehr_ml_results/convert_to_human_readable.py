# written to create output tables from the ATC and ICD10 codes
# written 5-19-21 JLW

import pickle,os,csv,json,math
from collections import defaultdict
import pandas as pd
import numpy as np

# load ATC codes
db_df_2 = pd.read_excel('../data/Drugbank050120.xlsx')
db_df_2 = db_df_2.replace(np.nan, 'NONE', regex=True)
db2atc = defaultdict(list)
atc2db = defaultdict(list)
db2name = pickle.load(open('../rscs/drugbankid_to_name.pkl','rb'))
db2name['DB14018'] = 'Bromotheophylline'
for (dbid,atc_entry) in zip(db_df_2['drugbank_id'].to_list(),db_df_2['atc_codes'].to_list()):
	if atc_entry=='NONE' or (type(atc_entry)==type(0.1) and math.isnan(atc_entry)):
		continue
	if '|' in atc_entry:
		for atcc in atc_entry.split('|'):
			db2atc[dbid].append(atcc)
			atc2db[atcc].append(dbid)
	else:
			db2atc[dbid] = [atc_entry]
			atc2db[atc_entry].append(dbid)

# remove ATC codes that map to more than one drug as we didn't want to consider combo products
atc2db_singles = dict([x for x in atc2db.items() if len(x[1])==1])

# manually curated ICD10 codes
cont_dmes_icd10 = dict([x for x in zip(['Delirium', 'Edema', 'Hypertension', 'Myopathy', 'Pancreatitis', 'Pneumonia', 'Sepsis'],['F05','R60','I10','G72,I42','K86,K85','J18,J15,J12,J16','A40,A41,B37.7'])])
icd10_to_dmes = dict([x for x in zip(['F05','R60','I10','G72,I42','K86,K85','J18,J15,J12,J16','A40,A41,B37.7'],['Delirium', 'Edema', 'Hypertension', 'Myopathy', 'Pancreatitis', 'Pneumonia', 'Sepsis'])])

# fcomp = 'target_comp_result.txt'
fcombo = 'aggregated_result.txt'
fcombo_lines = [l.strip() for l in open(fcombo,'r').readlines()]
out_lines = []
for row in fcombo_lines:        
	row_dic = json.loads(row)
	icd = row_dic["icd10codes"]
	[net_drugs,non_net_drugs] = row_dic["atc_codes"]
	combo_drug = list(set(row_dic['required_atc']))[0]
	net_drugs_atc = [atc2db_singles[x][0] for x in net_drugs if x in atc2db_singles]
	net_drugs_name = list(set([db2name[x] for x in net_drugs_atc]))

	non_net_drugs_atc = [atc2db_singles[x][0] for x in non_net_drugs if x in atc2db_singles]
	non_net_drugs_name = list(set([db2name[x] for x in non_net_drugs_atc]))
	combo_drug_name = db2name[atc2db[combo_drug][0]]
	dme_name = icd10_to_dmes[icd[0]]
	row_line = {'ExpNum':int(row_dic['exp_count'][0]),'DME':dme_name}
	row_line['Network Drugs'] = ','.join(net_drugs_name)
	row_line['Non Network Drugs'] = ','.join(non_net_drugs_name)

	row_line['Combo Drug'] = combo_drug_name
	out_lines.append(row_line)

converted_df = pd.DataFrame(out_lines).sort_values(by=['ExpNum'])
converted_df.to_excel('aggregated_result_drugnames.xlsx',index=False)

# repeat the process with the original experiment numbers for all 58 examples
f = 'aggregated_combinations_for_stride_ml.json'
f_lines = [l.strip() for l in open(f,'r').readlines()]
out_lines = []
for row in f_lines:
	row_dic = json.loads(row)
	icd = row_dic["icd10codes"]
	[net_drugs,non_net_drugs] = row_dic["atc_codes"]
	combo_drug = list(set(row_dic['required_atc']))[0]
	net_drugs_atc = [atc2db_singles[x][0] for x in net_drugs[0].split(',') if x in atc2db_singles]
	net_drugs_name = list(set([db2name[x] for x in net_drugs_atc]))
	non_net_drugs_atc = [atc2db_singles[x][0] for x in non_net_drugs[0].split(',') if x in atc2db_singles]
	non_net_drugs_name = list(set([db2name[x] for x in non_net_drugs_atc]))
	combo_drug_name = [db2name[atc2db_singles[cd][0]] for cd in combo_drug.split(',') if cd in atc2db_singles][0]
	dme_name = icd10_to_dmes[icd[0]]
	row_line = {'ExpNum':int(row_dic['exp_count'][0]),'DME':dme_name}
	row_line['Network Drugs'] = ','.join(net_drugs_name)
	row_line['Non Network Drugs'] = ','.join(non_net_drugs_name)

	row_line['Combo Drug'] = combo_drug_name
	out_lines.append(row_line)

converted_df = pd.DataFrame(out_lines).sort_values(by=['ExpNum'])
converted_df.to_excel('aggregated_combinations_for_stride_ml_drugnames.xlsx',index=False)

#        row_dic = json.loads(row)
#        if "prop_p" in row_dic:
#                exp_num = "exp"+row_dic["exp_count"][0]
#                prop_p = row_dic["prop_p"]
#                prop_coef = row_dic["prop_coef"]
#                combo_HRs[exp_num] = (prop_p,prop_coef)

#	# add p-values for where there were enough patients to measure HR
#	if "prop_p" in row_dic:
#		row_line['HR'] = math.exp(row_dic["prop_p"])
#		row_line['P-value'] = "{:.3e}".format(row_dic["prop_coef"])
#	else:
#		row_line['HR'] = 'not detected'
#
#
