# written to look at ORs before and after
# adding drug combinations
# written 3-17-21 JLW

import pickle,os,csv,json,matplotlib,math
import numpy as np
matplotlib.use("AGG")
import matplotlib.pyplot as plt

fcomp = 'target_comp_result.txt'
fcombo = 'aggregated_result.txt'

no_combo_HRs = {}
fcomp_lines = [l.strip() for l in open(fcomp,'r').readlines()]
for row in fcomp_lines:
	row_dic = json.loads(row)
	if "prop_p" in row_dic:
		exp_num = "exp"+row_dic["exp_count"][0]
		prop_p = row_dic["prop_p"]
		prop_coef = row_dic["prop_coef"]
		no_combo_HRs[exp_num] = (prop_p,prop_coef)

combo_HRs = {}
fcombo_lines = [l.strip() for l in open(fcombo,'r').readlines()]
for row in fcombo_lines:
	row_dic = json.loads(row)
	if "prop_p" in row_dic:
		exp_num = "exp"+row_dic["exp_count"][0]
		prop_p = row_dic["prop_p"]
		prop_coef = row_dic["prop_coef"]
		combo_HRs[exp_num] = (prop_p,prop_coef)


writer = pd.ExcelWriter('all_HR_pvalues.xlsx',engine='xlsxwriter')
nchrs = [[epn,math.exp(pp),"{:.3e}".format(cpv)] for (epn,(cpv,pp)) in no_combo_HRs.items()]
nchrs_df = pd.DataFrame(nchrs, columns = ['ExpNum','HR','P-value'])
nchrs_df.to_excel(writer,sheet_name = 'No Combos',index=False)
chrs = [[epn,math.exp(pp),"{:.3e}".format(cpv)] for (epn,(cpv,pp)) in combo_HRs.items()]
chrs_df = pd.DataFrame(chrs, columns = ['ExpNum','HR','P-value'])
chrs_df.to_excel(writer,sheet_name = 'Combos', index=False)
writer.save()


# method for labeling rectangles
def autolabel(rects):
	for rect in rects:
		height = rect.get_height()
		if height>0:
			ax.annotate('{:.2f}'.format(height),xy=(rect.get_x() + rect.get_width() / 2, height), xytext=(0, 3),textcoords="offset points",ha='center', va='bottom')
		if height<0:
			ax.annotate('{:.2f}'.format(height),xy=(rect.get_x() + rect.get_width() / 2, height), xytext=(0, -15),textcoords="offset points",ha='center', va='bottom')

# report table of p-values
hdr = ["Exp Number","Between Class HR","Between Class HR with combo","\n"]
outlines = "\t".join(hdr)
for (exp_num,(ppval,pcoef)) in no_combo_HRs.items():
	if exp_num in combo_HRs:
		cppval = combo_HRs[exp_num][0]
	else:
		cppval = 1.0
	out_data = [exp_num,"{:.3e}".format(ppval),"{:.3e}".format(cppval),"\n"]
	outlines = outlines + "\t".join(out_data)
print(outlines)
	

# plot coefficients
fig,ax = plt.subplots()
labels = [k for k in no_combo_HRs.keys()]
noCombos = []
combos = []
for k in labels:
	baseline = math.exp(no_combo_HRs[k][1])
	if k in combo_HRs:
		comboHR = math.exp(combo_HRs[k][1])
	else:
		comboHR = 0
	noCombos.append(baseline)
	combos.append(comboHR)

x = np.arange(len(labels))
width=0.35
r1 = ax.bar(x-width/2,noCombos,width,label='baseline',color='grey')
r2 = ax.bar(x+width/2,combos,width,label='with combo drugs',color='red')
ax.axhline(y=1.0, color='grey', linestyle=':')
ax.set_xticks(x)
ax.set_xticklabels(labels,rotation=90)
ax.legend()
ax.set_ylabel('cox coefficient')
ax.set_title('cox coefficient w/ and w/o combo drugs')
right_side = ax.spines["right"]
right_side.set_visible(False)
top_side = ax.spines["top"]
top_side.set_visible(False)
plt.subplots_adjust(bottom=0.25)
plt.savefig('stride_ml_cox_coef.png',format='png')

# plot differences
cc_diffs = [(exn,(cc-nc)) for (exn,nc,cc) in zip(labels,noCombos,combos) if cc!=0]
fig,ax = plt.subplots()
xlabels = [e for (e,ccd) in cc_diffs]
pos_c = [ccd if ccd>0 else 0 for (e,ccd) in cc_diffs]
neg_c = [ccd if ccd<0 else 0 for (e,ccd) in cc_diffs]
x = np.arange(len(xlabels))
r1 = ax.bar(x-width/2,pos_c,width,label='Increased Cox Coeff',color='red')
r2 = ax.bar(x+width/2,neg_c,width,label='Decreased Cox Coeff',color='dodgerblue')
ax.axhline(y=1.0, color='k')
ax.set_xticks(x)
ax.set_xticklabels(xlabels,rotation=90)
#ax.legend(loc='lower right')
ax.set_ylabel('cox coefficient')
ax.set_title('changes in cox coefficient')
ax.set_ylim([-0.5,0.5])
autolabel(r1)
autolabel(r2)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.subplots_adjust(bottom=0.25)
plt.savefig('stride_ml_changes_in_cox_coeff.png',format='png')

# plot the HR ratio
cc_fc = [(exn,(float(cc)/nc)) for (exn,nc,cc) in zip(labels,noCombos,combos) if cc!=0]
fig,ax = plt.subplots()
xlabels = [e for (e,ccd) in cc_diffs]
pos_c = [ccd if ccd>1.0 else 0 for (e,ccd) in cc_fc]
neg_c = [ccd if ccd<=1.0 else 0 for (e,ccd) in cc_fc]
x = np.arange(len(xlabels))
r1 = ax.bar(x-width/2,pos_c,width,label='Increased Cox Coeff',color='red')
r2 = ax.bar(x+width/2,neg_c,width,label='Decreased Cox Coeff',color='dodgerblue')
ax.axhline(y=1.0, color='k')
ax.set_xticks(x)
ax.set_xticklabels(xlabels,rotation=90)
ax.set_ylabel('HR Ratio')
ax.set_title('Fold change in cox coefficient')
ax.set_ylim([0,1.5])
autolabel(r1)
autolabel(r2)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.subplots_adjust(bottom=0.25)
plt.savefig('stride_ml_hr_ratios.png',format='png')
