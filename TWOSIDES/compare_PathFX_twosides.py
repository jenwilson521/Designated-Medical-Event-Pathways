# written to look at TWOSIDES data
# for drug interactions predictions
# written 11-20-19 JLW

import pickle,os,csv
from collections import defaultdict

# co-therapies from PathFX analysis
co_dir = '../data/cotherapy/'
co_f = os.path.join(co_dir,'activate_summary_data_for_EHR.csv')
#co_f = os.path.join(co_dir,'prevent_summary_data_for_EHR.csv')
dR = csv.DictReader(open(co_f,'r'))
pred_ddi = defaultdict(set)
for r in dR:
	pddi = r['Predicted Cotherapy'].lower()
	orig_drug = r['Drugs with labeled DME'].lower()
	dme = r['dme'].lower()
	pred_ddi[(orig_drug,pddi)].add(dme.lower())
print('number of interactions: '+str(len(pred_ddi)))


f = 'TWOSIDES.csv'
dR = csv.DictReader(open(f,'r'))
#outf2 = open('pathfs_vs_TWOSIDES_strict_compare_prevent.txt','w')
outf2 = open('pathfs_vs_TWOSIDES_strict_compare_activate.txt','w')
#hdr = ['drug1','drug2','DME','TWOSIDES condition name','PRR','\n']
#outf.write('\t'.join(hdr))
hdr = ['drug1','drug2','DME','TWOSIDES condition name','PRR','PRR error','mean reporting freq','\n']
outf2.write('\t'.join(hdr))
phens_to_map = set()
for r in dR:
	c1 = r['drug_1_concept_name'].lower()
	c2 = r['drug_2_concept_name'].lower()
	cn = r['condition_concept_name'].lower()
	prr = r['PRR']
	prr_err = r['PRR_error']
	mrf = r['mean_reporting_frequency']
	if (c1,c2) in pred_ddi or (c2,c1) in pred_ddi:
		if (c1,c2) in pred_ddi:
			dme = pred_ddi[(c1,c2)]
		else:
			dme = pred_ddi[(c2,c1)]
#		print((c1,c2,'dme: ',dme,'Twosides',cn))
#		outdata = [c1,c2,'|'.join(dme),cn,str(prr),'\n']
#		outf.write('\t'.join(outdata))
		for d in dme: # iterate through the set
			if d in cn:	
				outdata = [c1,c2,d,cn,str(prr),prr_err,mrf,'\n']
				outf2.write('\t'.join(outdata))

		else:
			for ddme in dme:
				phens_to_map.add((cn,ddme))
#outf3 = open('phenotypes_not_direct_matches.txt','w')
#for (cn,dme) in phens_to_map:
#	outf3.write('\t'.join([cn,dme,'\n']))

outf2.close()
#outf.close()

## development code
#all_cond = set()
#phen = 'myocardial infarction'
#phen_data = defaultdict(list) 
#phen_out = phen.replace(' ','_')
#pickle.dump(phen_data,open(phen_out+'.pkl','wb'))
#
## explore scores
#prrs = [float(x[2]) for x in phen_data['Acute myocardial infarction']]
#myo_fix = [(a.capitalize(),b.capitalize(),float(c)) for (a,b,c) in phen_data['Acute myocardial infarction']]
#rev_sort = sorted(myo_fix,key=lambda x:x[2],reverse=True)
#
#test_a = 'Dronabinol'
#test_b_all = 'Olanzapine,Sumatriptan,Lurasidone,Trimipramine,Zolmitriptan,Pindolol,Naloxone,Bromocriptine,Clozapine,Mirtazapine,Memantine,Maprotiline,Almotriptan,Pergolide,Amitriptyline,Pasireotide,Gabapentin,Trazodone,Pilocarpine,Naratriptan,Progesterone,Maraviroc,Octreotide,Amoxapine,Rizatriptan,Frovatriptan,Atropine,Pramipexole,Desipramine,Risperidone,Propranolol,Tinzaparin,Paroxetine,Dihydroergotamine,Aripiprazole,Pentazocine,Metoclopramide,Eletriptan,Menthol,Nicardipine,Perphenazine,Ropinirole,Ketoprofen,Dexamethasone,Nortriptyline,Adenosine,Apomorphine,Tropicamide,Disopyramide,Nalbuphine,Tramadol,Fentanyl,Diphenhydramine,Imipramine|Olanzapine,Sumatriptan,Trimipramine,Zolmitriptan,Naloxone,Bromocriptine,Clozapine,Maprotiline,Almotriptan,Pergolide,Amitriptyline,Naratriptan,Rizatriptan,Frovatriptan,Atropine,Pramipexole,Desipramine,Risperidone,Paroxetine,Aripiprazole,Pentazocine,Dihydroergotamine,Eletriptan,Nicardipine,Ropinirole,Nortriptyline,Apomorphine,Tropicamide,Nalbuphine,Tramadol,Fentanyl,Imipramine'
#
#test_b_list = list(set([x for short in test_b_all.split('|') for x in short.split(',')]))
#
#	
#check_ints = [x for x in rev_sort for b in test_b_list if test_a in x and b in x]
