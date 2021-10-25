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
all_dmes_net_classes = pd.DataFrame()
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
	co_df['AR'] = dme
	all_dmes_net_classes = pd.concat([all_dmes_net_classes,co_df]) # save this later for creating supplement file
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
ts_prr_data = []
for (d1,d2,ar,prr) in zip(tsides_short.drug_1_concept_name,tsides_short.drug_2_concept_name,tsides_short.condition_concept_name,tsides_short.PRR):
	d1_id = syns_to_dbid[d1.lower()]
	d2_id = syns_to_dbid[d2.lower()]
	dkey = '__'.join(sorted([d1_id,d2_id]))
	ts_all_combos.add(dkey)
	ts_combos_ars.add((dkey,ar))
	row_data = {'DrugKey':dkey,'drug_1_concept_name':d1,'drug_2_concept_name':d2,'condition_concept_name':ar,'PRR':prr}
	ts_prr_data.append(row_data)

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
pfx_in_tsides_matched_ar = set() # 50.22364217252397 <- 786/1565
tsides_full_keys = set()
for (dk,parr) in pfx_in_tsides_with_ars:
	for tsarr in ts_combos_ars_dic[dk]:
		if word_match(parr,tsarr):
			pfx_in_tsides_matched_ar.add((dk,parr))
			full_key = dk+'__'+tsarr
			tsides_full_keys.add(full_key)


# save all PRRs for supplement
ts_prr_df = pd.DataFrame(ts_prr_data) # (40675109, 13)
ts_prr_df = ts_prr_df[ts_prr_df['DrugKey'].isin(pfx_in_tsides)] # (213908, 5)
ts_prr_df['fullkey'] = ts_prr_df['DrugKey']+'__'+ts_prr_df['condition_concept_name']
ts_prr_df = ts_prr_df[ts_prr_df['fullkey'].isin(tsides_full_keys)] # (1742, 6)
ts_prr_df.to_excel('../TWOSIDES/all_pfx_pred_tsides_prr.xlsx',index=False)

# keep track of all network genes for faster lookup in final table
true_positives_dbid = pickle.load(open('../data/true_positives_dbid.pkl','rb'))
all_pos_dbs = list(set([x[0] for (ar,dlist) in true_positives_dbid.items() for x in dlist]))
nfdir = '../data/all_drugbank_network_files/'
nfs = [f for f in os.listdir(nfdir)]
nf_dic = dict([(f.split('_')[0],f) for f in nfs if 'merged_neighborhood_.txt' in f])
db_to_netgenes = defaultdict(set)
for (dbid,nf) in nf_dic.items():
	net_data = [l.strip().split('\t') for l in open(os.path.join(nfdir,nf)).readlines()]
	if net_data:
		(prota,protb,scr) = zip(*net_data)
		all_prot = prota+protb
		db_to_netgenes[dbid] =set(all_prot)

# write all predicted drug-drug-AR-network classes that overlap with TWOSIDES for supplement.
# add non-network, non-overlapping drug list too
all_class_pred = [] #non-overlapping targets, have an outcome in TWOSIDES
predicted_data = [(nnode,dc,nt,addrugs,ar) for (nnode,dc,nt,addrugs,ar) in zip(all_dmes_net_classes.index,all_dmes_net_classes['Drug'].to_list(), all_dmes_net_classes['Node Type'].to_list(), all_dmes_net_classes['All drugs'].to_list(), all_dmes_net_classes.AR)]
for (nnode,dc,nt,addrugs,ar) in predicted_data:
#	print((nnode,dc,nt,addrugs,ar))
	if dc.lower() not in syns_to_dbid:
		continue
	dc_dbid = syns_to_dbid[dc.lower()]
	db1 = syns_to_dbid[dc.lower()]
	keep_addrugs = set() 
	for adrug in addrugs:
		db2 = syns_to_dbid[adrug.lower()]
		dkey = '__'.join(sorted([db1,db2]))
		if (dkey,ar) in pfx_in_tsides_matched_ar:
#			print([dkey,ar])
			keep_addrugs.add(adrug)
	if len(keep_addrugs) > 0:
		keep_addrugs_dbid = set([syns_to_dbid[x.lower()] for x in keep_addrugs])
		all_true_positives_db = set([x[0] for x in true_positives_dbid[ar]])
		non_class_tps = all_true_positives_db.difference(keep_addrugs_dbid)
		non_class_non_overlap = set([dbid for dbid in non_class_tps if not overlap_targets('__'.join([dc_dbid,dbid]))]) 
		nonClass_nonOv_noNnode = [dbid for dbid in non_class_non_overlap if nnode not in db_to_netgenes[dbid]] # remove drugs with nnode in their network; some have the gene, but wasn't on shortest path between drug target and AR-assoc gene
		common_names = [dbid_to_all_names[dbid][0] for dbid in nonClass_nonOv_noNnode]
		row_data = {'AR':ar,'Predicted Combo':dc,'NetworkClassGene':nnode,'NodeType':nt,'AR-associated drugs':','.join(sorted(keep_addrugs)),'NonNetworkClass,AR-associated Drugs':','.join(sorted(common_names))} 
		all_class_pred.append(row_data)

final_class_pred = pd.DataFrame(all_class_pred)
final_class_pred.to_excel('all_SP_drug_class_predictions.xlsx',index=False)
		
