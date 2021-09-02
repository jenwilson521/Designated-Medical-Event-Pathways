# written to compare sensitivity of 
# PPI to random choices from tsides data
# written 6-9-21 JLW

import pickle, os, csv, random
import pandas as pd
from collections import defaultdict

# these are drug-drug combos from twosides where a drug or synonym is also in DrugBank, 193,960 combinations
tsides_short = pickle.load(open('../char_data/tsides_short.pkl','rb'))
all_tsides_combos = pickle.load(open('../char_data/all_twosides_combos_in_pathfx.pkl','rb')) # e.g., [('etanercept', 'omeprazole'), ('enalapril', 'petrolatum'), ('oxcarbazepine', 'vitamin e')]
all_conditions = list(set(tsides_short['condition_concept_name'].to_list()))

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

# mapping from synonyms to DBID
db_df = pd.read_csv('../data/drugbank_vocabulary.csv')
db_df = db_df.fillna('none')
all_syns_to_dbid = {}
for (dbid,cname,syns) in zip(db_df['DrugBank ID'].to_list(),db_df['Common name'].to_list(),db_df['Synonyms'].to_list()):
	all_syns_to_dbid[cname.lower()] = dbid
	if syns=='none':
		continue
	if ' | ' in syns:
		all_syns = syns.split(' | ')
		for asy in all_syns:
			all_syns_to_dbid[asy.lower()] = dbid
	else:
		all_syns_to_dbid[syns.lower()] = dbid

# method to look if targets overlap
dint = pickle.load(open('../rscs/drug_intome_targets.pkl','rb')) # DrugBank data from original pathways analysis
def is_overlapping(da,db):
	atargs = set(dint[da])
	btargs = set(dint[db])
	if len(atargs.intersection(btargs)) >= 1:
		return True
	else:
		return False

# map back to DrugBank ID
all_tsides_combos_dbid = [(all_syns_to_dbid[da],all_syns_to_dbid[db]) for (da,db) in all_tsides_combos]
# drug combos in tsides with drugbank entry and none-overlapping targets
tsides_dbid_non_ov = [(da,db) for (da,db) in all_tsides_combos_dbid if da in dint and db in dint and not is_overlapping(da,db)] # 86,292 combinations

# PPI predicted drugs with non-overlapping protein targets and dmes
non_overlapping = pickle.load(open('../char_data/ppi_predicted_non_overlapping_combos_with_dmes.pkl','rb')) # e.g., [('Aripiprazole', 'Drotrecogin alfa', 'Hyperlipidemia'), ('Droxidopa', 'Suramin', 'Hypertension'), 
ppi_dmes = list(set([sn for (da,db,sn) in non_overlapping]))
keep_tsides_conditions = [connm for connm in all_conditions for sn in ppi_dmes if word_match(connm,sn)] # 258 relevant conditions

tsides_rel_cond = tsides_short[tsides_short['condition_concept_name'].isin(keep_tsides_conditions)] # shape: (1480768, 13)
[da_all,db_all,cond_all] = [tsides_rel_cond['drug_1_concept_name'].to_list(),tsides_rel_cond['drug_2_concept_name'].to_list(),tsides_rel_cond['condition_concept_name'].to_list()]
tsides_real_answers = dict([(tuple(sorted([da,db])),cnm) for (da,db,cnm) in zip(da_all,db_all,cond_all)])

# now shuffle and ranomly pick with replacement to preserve reporting rates for certain conditions and measure sensitivity
def shuf(l):
	return random.sample(l,len(l))

def sample_rands(a,b,c):
	shuf_dic = dict([(tuple(sorted([a1,b1])),c1) for (a1,b1,c1) in zip(shuf(a),shuf(b),shuf(c))])
	return shuf_dic

def measure_sensitivity(dic1,dic2):
	tp = set(dic1.items()).intersection(set(dic2.items()))
	sens = float(len(tp))/len(dic2)
	return sens
	
rdir = '../char_data/shuffle_tsides/'
if not os.path.exists(rdir):
	os.makedirs(rdir)

all_sens = []
total_rand = 100
for i in range(100):
	if i%10 == 0:
		print(str(i)+' / '+str(total_rand))
	shuf_dic = sample_rands(da_all,db_all,cond_all)
	pickle.dump(shuf_dic,open(os.path.join(rdir,'tsides_shuf_dic_'+str(i)+'.pkl'),'wb'))
	sens = measure_sensitivity(shuf_dic,tsides_real_answers)
	# print(sens)
	all_sens.append(sens)

pickle.dump(all_sens,open(os.path.join(rdir,'tsides_suf_all_sens_values.pkl'),'wb'))

# plot and test for value
all_sens = pickle.load(open(os.path.join(rdir,'tsides_suf_all_sens_values.pkl'),'rb'))
