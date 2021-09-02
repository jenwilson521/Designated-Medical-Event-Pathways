# written 13-7-20 JLW to do a more 
# thorough characterization of predicted
# combinations and benchmark with TWOSIDES
# data

import pickle,os,csv, math
import pandas as pd
import numpy as np
from collections import defaultdict

# drug-DME from labels
drug_dmes = pickle.load(open('../data/true_positives_dbid.pkl','rb')) # value = (DBID, phen, phen_genes) # e.g., ('DB01223', 'Sepsis', 'CXCL8,POMC...'),()] # Some drugs are matched to multiple rele "phens"
drugs_matched_to_dmes = defaultdict(set)
dmes_to_drugs = defaultdict(set)
for (dme,drug_list) in drug_dmes.items():
	for (dbid,dphen,dgenes) in drug_list:
		drugs_matched_to_dmes[dbid].add(dme)
		dmes_to_drugs[dme].add(dbid)

co_f = '../data/cotherapy/potential_co_therapies.xlsx'
xl = pd.ExcelFile(co_f)
all_dmes = xl.sheet_names

keep_cols = ['Drug','Drug Target found in DME pathway']

# Save all combinations, keep track of the network proteins for class definition
all_combos = []
for sn in all_dmes:
	df = pd.read_excel(co_f,sheet_name =sn)
	df = df[keep_cols]
	labeled_drugs = drug_dmes[sn]
	for (dname,phen,phen_genes) in labeled_drugs:
		phen_genes_list = phen_genes.split(',')
		for (dcombo,net_prot) in zip(df['Drug'].to_list(),df['Drug Target found in DME pathway'].to_list()):
			if net_prot in phen_genes_list:
				all_combos.append((dname,dcombo,net_prot,sn,phen))

# select out unique combinations by DME category instead of specific DME phenotype
combos_DME_cat = [(dname,dcombo,sn) for (dname,dcombo,net_prot,sn,phen) in all_combos] # some "phens" map to the same "sn"/DME
network_classes = set([net_prot for (dname,dcombo,net_prot,sn,phen) in all_combos]) # 172 network proteins
combos_DMEs = set([sn for (dname,dcombo,net_prot,sn,phen) in all_combos]) # 12 DMEs
unique_combos = set(combos_DME_cat) # 18988, e.g., ('DB00652', 'Loperamide', 'Hypertension'), ('DB00434', 'Sr-123781a', 'Edema')

# now remove pairs if the drugs share protein targets
dint = pickle.load(open('../rscs/drug_intome_targets.pkl','rb')) # DrugBank data from original pathways analysis
def determine_overlap(d1,d2):
	if d1 in dint and d2 in dint:
		d1targs = set(dint[d1])
		d2targs = set(dint[d2])
		if len(d1targs.intersection(d2targs)) >=1:
			overlap = True
		else:
			overlap = False
	else: # if drug targets are unknown, return True to remove them from the rest of the analysis
		overlap = True
	return overlap

non_overlapping = [(dname,dcombo,sn) for (dname,dcombo,sn) in unique_combos if not determine_overlap(dname,dcombo)] # 6098 combos, also removed those without DB targs
just_drugs = list(set([(dname,dcombo) for (dname,dcombo,sn) in non_overlapping])) # 5246

# convert the DBIDs to common names for TWOSIDES analysis
db2name = pickle.load(open('../rscs/drugbankid_to_name.pkl','rb'))
non_overlapping_names = [(db2name[dname],dcombo,sn) for (dname,dcombo,sn) in non_overlapping]
pickle.dump(non_overlapping_names,open('../char_data/ppi_predicted_non_overlapping_names_combos_with_dmes.pkl','wb'))
pred_combs = [sorted([da.lower(),db.lower()]) for (da,db,sn) in non_overlapping_names] #6098

# load TWOSIDES and save unique combinations
tsf = '../TWOSIDES/TWOSIDES.csv'
tsdf = pd.read_csv(tsf)
tsds_combos = zip(tsdf['drug_1_concept_name'].to_list(),tsdf['drug_2_concept_name'].to_list())
tsides_all_conditions = set(tsdf['condition_concept_name'].to_list()) #
tsides_joinName = [sorted([da.lower(),db.lower()]) for (da,db) in tsds_combos] # 42,920,391 drug-drug-AR combos
tisdes_joinName_dbid = [(name2dbid[da],name2dbid[db]) for (da,db) in tsides_joinName if da in name2dbid and db in name2dbid] # 35,750,160
tisdes_joinName_non_overlapping = [(da,db) for (da,db) in tisdes_joinName_dbid if not determine_overlap(da,db)] # 14,380,015 non-overlapping
unique_total_tsides = list(set([(a,b) for (a,b) in tsides_joinName])) # 211,990 drug-drug combinations

# load DrugBank synomyns for comparison of drug combos
db_df = pd.read_csv('../data/drugbank_vocabulary.csv')
dsyns = db_df['Synonyms'].to_list()
drug_syns = defaultdict(list)
for dentry in dsyns:
	if type(dentry)==type('string'): 
		dlist = dentry.split(' | ')
		if len(dlist) >1:
			for ds in dlist[1:]:
				drug_syns[dlist[0].lower()].append(ds.lower())

# also load ATC codes
db_df_2 = pd.read_excel('../data/Drugbank050120.xlsx')
db_df_2 = db_df_2.replace(np.nan, 'NONE', regex=True)
db2atc = defaultdict(list)
for (dbid,atc_entry) in zip(db_df_2['drugbank_id'].to_list(),db_df_2['atc_codes'].to_list()):
	if atc_entry=='NONE' or (type(atc_entry)==type(0.1) and math.isnan(atc_entry)):
		continue
	if '|' in atc_entry:
		for atcc in atc_entry.split('|'):
			db2atc[dbid].append(atcc)
	else:
		db2atc[dbid] = [atc_entry]

# recreate drugname <-> DBID mappings to be all lower case
db2name = {} 
name2dbid = {}
for (db_no,dname) in zip(db_df_2['drugbank_id'].to_list(),db_df_2['name'].to_list()):
	db2name[db_no] = dname.lower()
	name2dbid[dname.lower()] =db_no 

# first look at combinations that are recorded in TWOSIDES regardless of the actual side effect
#pred_in_2sides = [] 
#for [dca,dcb] in pred_combs: 
#	if sorted([dca,dcb]) in tsides_joinName: # check if the sorted list in TWOSIDES
#		pred_in_2sides.append([dca,dcb])
#	else:
#		if dca in drug_syns:
#			all_a = [dca] + drug_syns[dca]
#		else:
#			all_a = [dca]
#		if dcb in drug_syns:
#			all_b = [dcb] + drug_syns[dcb]
#		else:
#			all_b = [dcb]
#		if len(all_a) > 1 or len(all_b) > 1: # if synonyms found, check all combination of synonyms
#			for dcasyn in all_a:
#				for dcbsyn in all_b:
#					if sorted([dcasyn,dcbsyn]) in tsides_joinName:
#						pred_in_2sides.append([dca,dcb])
#		else:
#			continue
#
#
#pickle.dump(pred_in_2sides,open('drug_combos_in_TWOSIDES.pkl','wb'))
pred_in_2sides = pickle.load(open('../char_data/drug_combos_in_TWOSIDES.pkl','rb')) #457, the original names, not the synonyms, these do not include drugs with shared proteins
pred_in_2_sides_set = set([tuple(sorted(skey)) for skey in pred_in_2sides]) #405
# pred_in_2_sides_set_dbids = [(name2dbid[da],name2dbid[db]) for (da,db) in pred_in_2_sides_set]
# pred_in_2_sides_non_overlapping_drugs = [(da,db) for (da,db) in pred_in_2_sides_set_dbids if not determine_overlap(da,db)] # [('goserelin', 'loperamide'), ('abatacept', 'celecoxib'),... # 404 drugs

# keep TWOSIDES drugs with names or synonyms in DrugBank
all_drug_names = set()
for (dn,dlist) in drug_syns.items():
	all_drug_names.add(dn.lower())
	for dls in dlist:
		all_drug_names.add(dls.lower())
for dn in db2name.values():
	all_drug_names.add(dn.lower())

# limit TWOSIDES data to drugs that are in my dataset to expedite searching
tsides_short = tsdf[((tsdf['drug_1_concept_name'].str.lower().isin(all_drug_names)) &(tsdf['drug_2_concept_name'].str.lower().isin(all_drug_names)))] # (40675109, 13)
pickle.dump(tsides_short, open('../char_data/tsides_short.pkl','wb'))
tsides_short = pickle.load(open('../char_data/tsides_short.pkl','rb'))

all_tsides_combos = set([tuple(sorted((da.lower(),db.lower()))) for (da,db) in zip(tsides_short['drug_1_concept_name'].to_list(),tsides_short['drug_2_concept_name'].to_list())]) # 193960 combinations
pickle.dump(all_tsides_combos,open('../char_data/all_twosides_combos_in_pathfx.pkl','wb'))
tsides_abbrev = [x for x in zip(tsides_short['drug_1_concept_name'].to_list(),tsides_short['drug_2_concept_name'].to_list(),tsides_short['condition_concept_name'].to_list())] # 40675109 drug-drug-ARs 
#... or just all of the combos that overlap with the drugnames that I studied

pre_ddi_in_twosides = [(da.lower(),db.lower(),sn.lower()) for (da,db,sn) in non_overlapping_names if sorted([da.lower(),db.lower()]) in pred_in_2sides] # 456 
pre_ddi_in_twosides_set = set([(tuple(sorted([da,db])),sn) for (da,db,sn) in pre_ddi_in_twosides]) # 456
# expand predictions to include synonyms that could have matched case reports in TWOSIDES
pre_ddic_with_syns = {}
syns_to_orig_combo = {}
for (da,db,sn) in pre_ddi_in_twosides: 
	orig_tup = tuple(sorted([da,db]))
	pre_ddic_with_syns[orig_tup] = sn
	if da in drug_syns:
		all_a = [da] + drug_syns[da]
	else:
		all_a = [da]
	if db in drug_syns:
		all_b = [db] + drug_syns[db]
	else:
		all_b = [db]
	if len(all_a) >1 or len(all_b) >1:
		for das in all_a:
			for dbs in all_b:
				syn_tup = tuple(sorted([das,dbs]))
				pre_ddic_with_syns[syn_tup] = sn
				syns_to_orig_combo[syn_tup] = orig_tup

# search drug-drug-DME sets for overlap in abbreviated TWOSIDES data
def word_match(s1,s2):
	wmatch = False
	if s1.lower() in s2.lower() or s2.lower() in s1.lower():
		wmatch = True
	else:
		s1_words = [s.lower() for s in s1.split(' ')]
		s2_words = [s.lower() for s in s2.split(' ')]
		if len(set(s1_words).intersection(set(s2_words))) >= 1:
			wmatch = True
	return wmatch
	
# explore matches with just TWOSIDES data, just exploratory analysis
ddi_wth_2sides_evidence = [] # 368
explore_dic = defaultdict(list)
tisdes_sort_key = [tuple(sorted([da,db])) for (da,db,dme_name) in tsides_abbrev]
for (da,db,dme_name) in tsides_abbrev:
	sort_key = tuple(sorted([da.lower(),db.lower()]))
	if sort_key in pre_ddic_with_syns:
		pred_dme = pre_ddic_with_syns[sort_key]
		# if dme_name.lower() in pred_dme.lower() or pred_dme.lower() in dme_name.lower():
		if word_match(dme_name.lower(),pred_dme.lower()):
			if sort_key in syns_to_orig_combo: # convert back to orginal drug combos if synomyns are used
				save_key = syns_to_orig_combo[sort_key]
			else:
				save_key = sort_key
			ddi_wth_2sides_evidence.append((save_key,pred_dme))
		explore_dic[(sort_key)].append(dme_name) # just look at how many side effects are associated with a combination

# method to split apart combinations
def split_drug_str(dstring):
	all_drugs = set()
	if '|' in dstring and ',' in dstring:
		all_drugs = set([d.lower() for sub_str in dstring.split('|') for d in sub_str.split(',')])
	elif '|' not in dstring and ',' in dstring:
		all_drugs = set([d.lower() for d in dstring.split(',')])
	return all_drugs

def convert_to_atc_str(dset):
	all_drugs = set()
	for dname in dset:
		dbid = name2dbid[dname.lower()]
		atc_codes = db2atc[dbid] # this is a list, but may be empty because it's a default dict
		all_drugs = all_drugs.union(set(atc_codes)) 
#		if atc_codes==[]:
#			all_drugs.add("FIX:"+dname)
#		else:
#			all_drugs = all_drugs.union(set(atc_codes)) 
	return ','.join(list(all_drugs))

def return_non_overlapping(cdrug,cset):
	comparator_drugs = set()
	comb_targs = set(dint[name2dbid[cdrug]])
	for cmpd in cset:
		cmp_targs = set(dint[name2dbid[cmpd]])	
		if len(comb_targs.intersection(cmp_targs))<1:
			comparator_drugs.add(cmpd)
	return comparator_drugs

def expand_pairs_with_syns(drug_name):
	if drug_name in drug_syns:
		all_names = [drug_name] + drug_syns[drug_name]
	else:
		all_names = [drug_name]
	return all_names
	
def check_drug_pair_match(npair,pair_list): # need a separate method  to compare all synonyms
	nmatch = False
	(da,db) = npair
	all_a = expand_pairs_with_syns(da)
	all_b = expand_pairs_with_syns(db)
	new_pairs = []
	for da1 in all_a:
		for db1 in all_b:
			sort_key = tuple(sorted((da1,db1)))
			new_pairs.append(sort_key)
	if len(set(new_pairs).intersection(set(pair_list))) > 1:
		nmatch = True
	return nmatch	

# aggregate combinations based on network genes
# combinations need to share a network gene for this analysis
tps = pickle.load(open('../data/true_positives_dbid.pkl','rb')) # all drugs that had a network connection to the DME
net_dmes = [k for k in tps.keys()]
all_tsds_dmes = list(set([dme for (dkey,dme) in ddi_wth_2sides_evidence]))
cont_dmes = [ndme for ndme in net_dmes if ndme.lower() in all_tsds_dmes] # only pursue the overlapping ones
cont_dmes_icd10 = dict([x for x in zip(['Delirium', 'Edema', 'Hypertension', 'Myopathy', 'Pancreatitis', 'Pneumonia', 'Sepsis'],['F05','R60','I10','G72,I42','K86,K85','J18,J15,J12,J16','A40,A41,B37.7'])])
# clean up the all_combos object to group combos by mechanism
all_combo_clean_and_sort = []
for (dbid,d2name,nmech,ndme,dme_phen) in all_combos:
	db_name = db2name[dbid].lower()
	d2name = d2name.lower()
	sort_key = tuple(sorted((db_name,d2name)))
	all_combo_clean_and_sort.append((sort_key,nmech,ndme.lower()))


outf = open('aggregated_combinations_for_stride_ml.json','w')
outf2 = open('just_target_comp_for_stride_ml.json','w')
outf3 = open('aggregated_combinations_human_readable.json','w')
out_line = '{"icd10codes":["ICD10CODES"],"atc_codes": [["TARGETDRUGS"],["COMPARATORDRUGS"]],"required_atc": ["COMBODRUG"],"exp_count":["EXPCOUNT"]}'
exp_count = 0
net_mech_rows = []
for ndme in cont_dmes:
	if ndme == 'Proteinuria':
		continue
	out_icd10 = cont_dmes_icd10[ndme]
	all_tps = set([db2name[dbid] for (dbid,dphen,genelist) in tps[ndme]])
	ddis = list(set([tuple(sorted((da,db))) for ((da,db),dme) in ddi_wth_2sides_evidence if ndme.lower()==dme.lower()])) # there will be repeats if drug synonyms matched TWOSIDES combos
	combos_with_mechs = [(skey,nmech,netdme) for (skey,nmech,netdme) in all_combo_clean_and_sort if skey in ddis] 
	all_mechs = list(set([nmech for (skey,nmech,netdme) in combos_with_mechs]))
	for smech in all_mechs: # separate analysis based on network protein
		drug_combos = [skey for (skey,nmech,netdme) in combos_with_mechs if nmech==smech]
		all_considered_drugs = set([dn for skey in drug_combos for dn in skey])
		combo_drugs = all_considered_drugs.difference(all_tps) # back-out which are the combo drugs
		for cdrug in combo_drugs:
			out_req_str = convert_to_atc_str({cdrug})
			smech_dic = {'DownstreamProtein':smech,'DME':ndme,'ComboATC':out_req_str}

			pairs_with_combo = [skey for skey in drug_combos if skey[0]==cdrug or skey[1]==cdrug]
			target_set = set([dn for skey in pairs_with_combo for dn in skey if dn!=cdrug])
			out_targ_str = convert_to_atc_str(target_set)
			smech_dic['TargetDrugs'] = out_targ_str
			
			# find comparators by selecting TPs w/o nmech that are also co-prescribed in TWOSIDES
			pos_net_no_mech = set([db2name[dbid] for (dbid,dphen,genelist) in tps[ndme] if smech not in genelist]).difference(target_set)
			comp_combos = set([pnnm for pnnm in pos_net_no_mech if tuple(sorted([cdrug,pnnm])) in all_tsides_combos])
			out_comp_str = convert_to_atc_str(comp_combos)
			smech_dic['ComparatorDrugs'] = out_comp_str
			net_mech_rows.append(smech_dic)

			exp_count+=1
			final_out_line = out_line.replace('ICD10CODES',out_icd10).replace('TARGETDRUGS',out_targ_str).replace('COMPARATORDRUGS',out_comp_str).replace('COMBODRUG',out_req_str).replace('EXPCOUNT',str(exp_count))
			n=outf.write(final_out_line+'\n') # with the combo drug
			no_combo_out_line = final_out_line.replace(out_req_str,'') # without the combo drug, the baseline hazard ratio
			n=outf2.write(no_combo_out_line+'\n')
			readable_line = out_line.replace('ICD10CODES',out_icd10+' '+ndme).replace('TARGETDRUGS',','.join(list(target_set))).replace('COMPARATORDRUGS',','.join(list(comp_combos))).replace('COMBODRUG',cdrug).replace('EXPCOUNT',str(exp_count))
			n=outf3.write(readable_line+'\n')
			
outf.close()
outf2.close()
outf3.close()
net_mechs_df = pd.DataFrame(net_mech_rows)
net_mechs_df.to_excel('../char_data/network_mechanisms_for_ehr_ml.xlsx',index=False)

### Literature comparisons
# further subset drug-drug-AR pairs based on literature-derived directionality
co_ther_file = '../data/cotherapy/supp4_summary_drug_interactions_noSharedTargets.xlsx'
xl = pd.ExcelFile(co_ther_file)
# all_rels = set()
ddi_rels = ['Aggravates','Activates','Induces']
co_ther = ['Inhibits','Prevents','Prevcents']
lit_agg = set()
lit_inhib = set()
for dme_name in xl.sheet_names:
	df = xl.parse(sheet_name=dme_name)
	df_agg = df[df['Relationship'].isin(ddi_rels)]
	for (da,db) in zip(df_agg['Drugs with labeled DME'],df_agg['Predicted Cotherapy']):
		sort_key = tuple(sorted([da.lower(),db.lower()]))
		lit_agg.add((sort_key,dme_name.lower()))
	df_inh = df[df['Relationship'].isin(co_ther)]
	for (da,db) in zip(df_inh['Drugs with labeled DME'],df_inh['Predicted Cotherapy']):
		sort_key = tuple(sorted([da.lower(),db.lower()]))
		lit_inhib.add((sort_key,dme_name.lower()))
#	drels = df['Relationship'].to_list()
#	all_rels = all_rels.union(set(drels))

# expand drug names because of formatting
lit_agg_format = set()
for ((dat,dbt),dme_nm) in lit_agg:
	if '|' in dat or ',' in dat:
		drug_splits = [xx for x in dat.split('|') for xx in x.split(',')]
		for dasplit in drug_splits:
			sort_key = tuple(sorted([dasplit,dbt]))
			lit_agg_format.add((sort_key,dme_nm))
	elif '|' in dbt or ',' in dbt:
		drug_splits = [xx for x in dbt.split('|') for xx in x.split(',')]
		for dbsplit in drug_splits:
			sort_key = tuple(sorted([dat,dbsplit]))
			lit_agg_format.add((sort_key,dme_nm))
	else:
		sort_key = tuple(sorted([dat,dbt]))
		lit_agg_format.add((sort_key,dme_nm))

lit_agg_dic = dict(lit_agg_format)
pre_ddi_lit_agg = list(set([(da,db,dname) for (da,db,dname) in pre_ddi_in_twosides if tuple(sorted([da,db])) in lit_agg_dic]))
pre_ddi_lit_agg_syns_draft = defaultdict(set)
for (da,db,dname) in pre_ddi_lit_agg:
	pre_ddi_lit_agg_syns_draft[tuple(sorted([da,db]))].add(dname)
###pre_ddi_lit_agg_syns_draft = dict([(tuple(sorted([da,db])),dname) for (da,db,dname) in pre_ddi_lit_agg]) # seed dictionary with direct names
pre_ddi_lit_agg_syns = defaultdict(set) 
for ((da,db),dset) in pre_ddi_lit_agg_syns_draft.items():
	if da in drug_syns:
		all_a = [da] + drug_syns[da]
	else:
		all_a = [da]
	if db in drug_syns:
		all_b = [db] + drug_syns[db]
	else:
		all_b = [db]
	if len(all_a) >1 or len(all_b) >1:
		for das in all_a:
			for dbs in all_b:
				syn_tup = tuple(sorted([das,dbs]))
				pre_ddi_lit_agg_syns[syn_tup] = dset
pre_ddi_lit_agg_syns.update(pre_ddi_lit_agg_syns_draft)
lit_agg_and_tsides = set() # ends up being only 8 combinations
for (da,db,dme_name) in tsides_abbrev:
	sort_key = tuple(sorted([da.lower(),db.lower()]))
	if sort_key in pre_ddi_lit_agg_syns:
		dme_set = pre_ddi_lit_agg_syns[sort_key]
		for pred_name in dme_set:
			if word_match(dme_name.lower(),pred_dme.lower()):
				if sort_key in syns_to_orig_combo:
					save_key = syns_to_orig_combo[sort_key]
				else:
					save_key = sort_key
				lit_agg_and_tsides.add((sort_key,pred_dme))

### Development ###
# pre_in_2sides = set(pred_combs).intersection(set(tsides_joinName))

# convert to ATC codes
#outf = open('lit_and_twosides.txt','w')
#outline = '{"icd10codes":["I10"], "atc_codes":[[DRUGA],[DRUGB]],"required_atc":[COMBODRUG]}\n' # since they are all hypertension, use the same ITC Code I10
#
#for ((da,db),dme) in lit_agg_and_tsides:
#	da_dbid = name2dbid[da.lower()]
#	db_dbid = name2dbid[db.lower()]
#	da_atcs = db2atc[da_dbid]
#	db_atcs = db2atc[db_dbid] 
#	for datc in da_atcs:
#		for dbtc in db_atcs:
#			n_out = outline.replace('DRUGA',datc).replace('DRUGB',datc).replace('COMBODRUG,dbtc) ### FIX THIS TO FIT STRIDE_ML format
#			outf.write(n_out) 
#
#	#### DIDN'T FINISH YET


