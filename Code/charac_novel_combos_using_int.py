# written to count DDIs using intermediate genes
# as well as terminal genes
# written 10-19-21 JLW

import pickle,os,csv, math
import pandas as pd
import numpy as np
from collections import defaultdict

# load drugbank synonyms for later
db_df = pd.read_csv('../data/drugbank_vocabulary.csv')
db_df = db_df.fillna('none')
dbid_to_all_names = defaultdict(list)
syns_to_dbid = {}
for (dbid,cn,syn_str) in zip(db_df['DrugBank ID'].to_list(),db_df['Common name'].to_list(),db_df['Synonyms'].to_list()):
	dbid_to_all_names[dbid].append(cn.lower())
	syns_to_dbid[cn.lower()] = dbid
	if syn_str == 'none':
		continue
	for dsyn in syn_str.split(' | '):
		dbid_to_all_names[dbid].append(dsyn.lower())
		syns_to_dbid[dsyn.lower()] = dbid

# Get intermediate nodes and in which networks they are found, per DME
ncf = 'supp2_DME_merged_node_counts.xlsx' 
nc_xl = pd.ExcelFile(ncf)

# get predicted co-therpies based on if their targets are in the merged DME
# networks
co_f = '../data/cotherapy/potential_co_therapies.xlsx'
co_xl = pd.ExcelFile(co_f) # note: these sheetnames have _ and no spaces

# method for splitting all drugs and returning a single set
def split_drugs(dstr):
	drugs_phens = dstr.split(';')
	(drugs,phens) = zip(*[tuple(x.split(':')) for x in drugs_phens])
	return set(drugs)

# count combinations when intermediate genes along the shortest path are used to predict combinations
all_net_classes = set()
all_combinations = set()
combos_with_dmes = set()
for dme in nc_xl.sheet_names:
	# dme = 'Thrombocytopenia'
	nc_df = nc_xl.parse(dme)
	nc_df = nc_df.set_index('Node Name').drop(['Unnamed: 0'],axis=1)
	nc_df['All drugs'] = nc_df['Included in these drug pathways'].apply(lambda x: split_drugs(x))
	nc_df['Unique Drug Count'] = nc_df['All drugs'].apply(lambda x: len(x))
	targs2type = dict(zip(nc_df.index,nc_df['Node Type'].to_list()))
	targs2dmedrugs = dict(zip(nc_df.index,nc_df['All drugs'].to_list()))
	co_df = co_xl.parse(dme.replace(' ','_')) # some genes are listed in multiple rows because there are multiple drugs that target this gene
	co_df = co_df.set_index('Drug Target found in DME pathway').drop(['Unnamed: 0'],axis=1)
	co_df['Node Type'] = co_df.index.map(targs2type)
	co_df['All drugs'] = co_df.index.map(targs2dmedrugs)
	for net_node in co_df.index:
		all_net_classes.add((dme,net_node))
	dme_net_classes = set(co_df.index)
	print(dme)
	print('no. of net classes: '+str(len(dme_net_classes)))
	for (net_node,cdrug,dme_drug_set) in zip(co_df.index,co_df['Drug'].to_list(),co_df['All drugs'].to_list()):
		for dme_drug in dme_drug_set:
			combo_key = '__'.join(sorted([cdrug,dme_drug]))
			all_combinations.add(combo_key)
			combos_with_dmes.add((combo_key,dme))



# how many DMEs had classifications
net_class_dmes = list(set([d for (d,g) in all_net_classes]))

# original drugbank binding information
dint = pickle.load(open('../rscs/drug_intome_targets.pkl','rb'))

# method to get targets, use synonyms if needed
def get_targets(dn):
	if dn in dint:
		dtargs = dint[dn]
	elif dn.lower() in syns_to_dbid and syns_to_dbid[dn.lower()] in dint:
		dbid = syns_to_dbid[dn.lower()]
		dtargs = dint[dbid]
	else:
		dtargs = []
	return dtargs

# method to look up drugs and exclude those with any shared targets
def overlap_targets(dcombo):
	# dcombo = 'Brompheniramine__Sitaxentan'
	shard = True
	[d1,d2] = dcombo.split('__')
	d1targs = get_targets(d1)
	d2targs = get_targets(d2)
	if len(d1targs) > 1 and len(d2targs) >1: # must have documented targets and then non-overlapping. Skip drugs without targets
		ovlp = set(d1targs).intersection(set(d2targs))
		if len(ovlp) <1:
			shard = False
	return shard

# subset the data based on no shared proteins
combos_not_shared = [x for x in all_combinations if not overlap_targets(x)] # 11904 drug-drug combinations
combos_with_dmes_not_shared = [x for x in combos_with_dmes if x[0] in combos_not_shared] # 19741 drug-drug-ARs

# filtered TWOSIDES data based on drugs in this analysis
tsides_short = pickle.load(open('../char_data/tsides_short.pkl','rb'))
ts_drugs = list(set(list(tsides_short.drug_1_concept_name) + list(tsides_short.drug_2_concept_name)))

# convert both datasets to DBIDs to facilitate matching and comparison
ts_in_dbid = [d for d in ts_drugs if d.lower() in syns_to_dbid] # 1722
ts_all_combos = set()
ts_combos_ars = set()
for (d1,d2,ar) in zip(tsides_short.drug_1_concept_name,tsides_short.drug_2_concept_name,tsides_short.condition_concept_name):
	d1_id = syns_to_dbid[d1.lower()]
	d2_id = syns_to_dbid[d2.lower()]
	dkey = '__'.join(sorted([d1_id,d2_id]))
	ts_all_combos.add(dkey)
	ts_combos_ars.add((dkey,ar))

# some drug names are in dint but not in syns_to_dbid?
all_pfx_drugs = list(set([d for dc in combos_not_shared for d in dc.split('__')]))
pfx_not_in_syns = [d for d in all_pfx_drugs if d.lower() not in syns_to_dbid] # 4 drugs
syns_to_dbid['fondaparinux sodium'] = syns_to_dbid['fondaparinux']
syns_to_dbid['ipratropium bromide'] = syns_to_dbid['ipratropium']
syns_to_dbid['menthol'] = syns_to_dbid['menthol natural']
syns_to_dbid['isoetarine'] = 'DB00221'

# convert network-class predictions to DBIDs
combos_not_shared_dbid = set()
combos_with_dmes_not_shared_dbid = set()
for (c,parr) in combos_with_dmes_not_shared:
	[d1,d2] = c.split('__')
	d1_db = syns_to_dbid[d1.lower()]
	d2_db = syns_to_dbid[d2.lower()]
	dkey = '__'.join(sorted([d1_db,d2_db]))
	combos_not_shared_dbid.add(dkey)
	combos_with_dmes_not_shared_dbid.add((dkey,parr))

# method for looking at synonymous adverse reactions
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

# compare combos in TWOSIDES
pfx_in_tsides = combos_not_shared_dbid.intersection(ts_all_combos) # how many combos are prescribed in the real world = 964
pfx_in_tsides_with_ars = list(set([x for x in combos_with_dmes_not_shared_dbid if x[0] in pfx_in_tsides])) # 1565
ts_combos_ars_dic = defaultdict(list)
for (dk,tsarr) in ts_combos_ars:
	ts_combos_ars_dic[dk].append(tsarr) 

# look if the ARs are similar between datasets
pfx_in_tsides_matched_ar = set() # 50.22364217252397 -> 786/1565
for (dk,parr) in pfx_in_tsides_with_ars:
	for tsarr in ts_combos_ars_dic[dk]:
		if word_match(parr,tsarr):
			pfx_in_tsides_matched_ar.add((dk,parr))



