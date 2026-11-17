"""Publication figures drawn exclusively from exported survey results."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];F=P/'figures'
NAVY='#153D64';BLUE='#78AED3';TEAL='#229C91';PALE='#A9D6D0';GRAY='#DCE2E8';DARK='#53616C'
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.titlesize':15,'axes.titleweight':'bold','axes.labelsize':13,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#4B5156','xtick.color':'#50565B','ytick.color':'#50565B','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white','axes.unicode_minus':True})
profiles=pd.read_csv(P/'results/profiles.csv');adj=pd.read_csv(P/'results/adjusted_contrasts.csv');sensitivity=pd.read_csv(P/'results/sensitivities.csv')
def val(year,group,measure,data=profiles):return data[(data.year==year)&data.group.eq(group)&data.measure.eq(measure)].iloc[0]
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(F/(name+'.'+ext),dpi=600,bbox_inches='tight',pad_inches=.14)
 plt.close(fig)
def clean(ax):ax.set_axisbelow(True);ax.grid(axis='x',color='#E5E8EB',lw=.8);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
def ci(ax,x,y,lo,hi,color,marker='o',filled=True):ax.errorbar(x,y,xerr=[[x-lo],[hi-x]],fmt=marker,color=color,mfc=color if filled else 'white',mec=color,ms=7.5,lw=2.1,capsize=0,zorder=3)
def f1():
 fig,axs=plt.subplots(1,2,figsize=(14,6.2),gridspec_kw={'width_ratios':[1,1.18]});fig.subplots_adjust(wspace=.35,bottom=.23)
 a=axs[0];names=['Dilation','Photography','Either'];colors=[BLUE,TEAL,NAVY]
 for i,(n,c) in enumerate(zip(names,colors)):
  r=val(2025,'All',n);a.bar(i,r.estimate,width=.6,color=c);a.errorbar(i,r.estimate,yerr=[[r.estimate-r.lo],[r.hi-r.estimate]],color='#222',capsize=4,lw=1.7);a.text(i,r.hi+1.7,f'{r.estimate:.1f}',ha='center',weight='bold',size=15)
 a.set(ylim=(0,85),ylabel='Reported annual contact, %',xticks=range(3),xticklabels=names);a.set_title('A  Eye care by survey definition',loc='left',pad=17);a.spines['top'].set_visible(False);a.set_axisbelow(True);a.yaxis.grid(True,color='#E4E7EA');a.text(.02,.94,'2025  |  17,806 complete pairs',transform=a.transAxes,color=DARK,size=11)
 a=axs[1];cats=['Both','Dilation only','Photography only','Neither'];colors=[NAVY,BLUE,TEAL,GRAY]
 for row,y in enumerate([2025,2022]):
  start=0
  for n,c in zip(cats,colors):
   r=val(y,'All',n);a.barh(row,r.estimate,left=start,height=.5,color=c,edgecolor='white',lw=1.5);a.text(start+r.estimate/2,row,f'{r.estimate:.1f}',ha='center',va='center',color='white' if n in ['Both','Photography only'] else '#26333E',size=12,weight='bold');start+=r.estimate
 a.set(yticks=[0,1],yticklabels=['2025\n20 states','2022\n12 states'],xlim=(0,100),ylim=(1.65,-.65),xlabel='Complete-pair respondents, %');a.set_title('B  Joint annual contact patterns',loc='left',pad=17);a.spines['left'].set_visible(False);a.tick_params(axis='y',length=0);a.text(.01,.02,'Photography only adds 4.7 percentage points in 2025',transform=a.transAxes,size=11,color=DARK)
 handles=[plt.Line2D([0],[0],marker='s',color=c,ms=10,lw=0) for c in colors];fig.legend(handles,cats,loc='lower center',bbox_to_anchor=(.74,.02),ncol=2,frameon=False,fontsize=11)
 save(fig,'Figure1_paired_eye_care_classification')
def f2():
 groups=[('All','All adults'),('nonmetro: 0.0','Metropolitan'),('nonmetro: 1.0','Nonmetropolitan'),('age_group: 18–44','Age 18–44'),('age_group: 45–64','Age 45–64'),('age_group: ≥65','Age ≥65'),('race: White NH','White, non-Hispanic'),('race: Black NH','Black, non-Hispanic'),('race: Hispanic','Hispanic'),('cost: 0.0','No cost barrier'),('cost: 1.0','Cost barrier'),('vision: 0.0','No serious vision difficulty'),('vision: 1.0','Serious vision difficulty')]
 fig,axs=plt.subplots(1,2,figsize=(15,9));fig.subplots_adjust(wspace=.48,left=.23,right=.87,bottom=.12,top=.9)
 for ax,measure,title,color,maxx in zip(axs,['Photography only','Photo among no dilation'],['A  Increment among all complete pairs','B  Photography without annual dilation'],[NAVY,TEAL],[11.5,29]):
  for i,(g,label) in enumerate(groups):
   r=val(2025,g,measure);ci(ax,r.estimate,i,r.lo,r.hi,color);ax.text(1.04,i,f'{r.estimate:.1f} ({r.lo:.1f}–{r.hi:.1f})',transform=ax.get_yaxis_transform(),va='center',size=10.5)
  ax.axvline(val(2025,'All',measure).estimate,color=color,lw=1,ls=(0,(4,4)),alpha=.55);ax.set(ylim=(len(groups)-.5,-.7),xlim=(0,maxx),yticks=range(len(groups)),yticklabels=[l for _,l in groups] if ax is axs[0] else ['']*len(groups),xlabel='Survey-weighted percentage (95% CI)');clean(ax);ax.set_title(title,loc='left',pad=18,size=13);ax.text(1.04,1.025,'Estimate (95% CI)',transform=ax.transAxes,size=10.5,weight='bold')
 fig.text(.23,.037,'2025 BRFSS  |  See Table S2a for subgroup denominators. Panel B is conditional on no reported annual dilation.',size=11,color=DARK);save(fig,'Figure2_photography_only_by_group')
def f3():
 fig,axs=plt.subplots(1,2,figsize=(14,6.6));fig.subplots_adjust(left=.18,right=.84,wspace=.7,top=.82,bottom=.30)
 for ax,exposure,title in zip(axs,['nonmetro','cost'],['A  Nonmetropolitan versus metropolitan','B  Cost barrier versus no cost barrier']):
  for i,(y,o,label,c) in enumerate([(2025,'dil','2025  Dilation',NAVY),(2025,'either','2025  Either modality',TEAL),(2022,'dil','2022  Dilation',BLUE),(2022,'either','2022  Either modality',PALE)]):
   r=adj[(adj.year==y)&adj.model.eq('M1')&adj.outcome.eq(o)&adj.exposure.eq(exposure)].iloc[0];ci(ax,r.difference_pp,i,r.lo,r.hi,c);ax.text(1.04,i,f'{r.difference_pp:.1f} ({r.lo:.1f} to {r.hi:.1f})',transform=ax.get_yaxis_transform(),va='center',size=10.5)
  ax.set(yticks=range(4),yticklabels=['2025  Dilation','2025  Either modality','2022  Dilation','2022  Either modality'] if ax is axs[0] else ['']*4,ylim=(3.5,-.5),xlim=(-25,25),xticks=[-20,-10,0,10,20],xlabel='Adjusted difference, percentage points');clean(ax);ax.axvline(0,color='#888',ls=(0,(4,4)),lw=1);ax.set_title(title,loc='left',pad=22,size=13)
  r=adj[(adj.year==2025)&adj.model.eq('M1')&adj.outcome.eq('Change in disparity')&adj.exposure.eq(exposure)].iloc[0];ax.text(0,-.27,f'2025 gap change  {r.difference_pp:+.2f} pp\n95% CI {r.lo:.2f} to {r.hi:.2f}',transform=ax.transAxes,color=DARK,size=11)
 fig.text(.18,.025,'Standardized differences from one four-category model per year.\nThe either-modality gap minus the dilation gap equals the photography-only difference.',size=10.5,color=DARK);save(fig,'Figure3_adjusted_contact_disparities')
def fs1():
 flow=json.loads((P/'results/flow.json').read_text());fig,axs=plt.subplots(1,2,figsize=(12,8));fig.subplots_adjust(wspace=.18)
 for ax,r in zip(axs,flow['samples']):
  ax.axis('off');ax.set_title(str(r['year'])+' BRFSS',loc='center',pad=15)
  boxes=[f'Public survey file\nn = {r["full_records"]:,}',f'{r["states"]} full-module US states\nAll respondents n = {r["state_sample"]:,}',f'Diagnosed nongestational diabetes\nn = {r["diabetes"]:,}',f'Complete eye-care recency pair\nn = {r["complete_pair"]:,}']
  for j,txt in enumerate(boxes):
   y=.84-j*.22;ax.add_patch(FancyBboxPatch((.07,y-.09),.86,.14,boxstyle='round,pad=.015,rounding_size=.008',facecolor='#EDF4F9' if j<3 else '#DDECE9',edgecolor=NAVY,lw=1.1,transform=ax.transAxes));ax.text(.5,y-.02,txt,ha='center',va='center',size=13,transform=ax.transAxes)
   if j<3:ax.annotate('',xy=(.5,y-.165),xytext=(.5,y-.105),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':NAVY,'lw':1.3})
  ax.text(.5,.01,f'Incomplete pairs: {r["diabetes"]-r["complete_pair"]:,}\nRetained as a separate missingness domain',ha='center',va='center',size=11,color=DARK,transform=ax.transAxes)
 save(fig,'FigureS1_participant_flow')
def fs2():
 data=pd.read_csv(P/'results/state_profiles.csv');data=data[(data.year==2025)&data.measure.eq('Photography only')].sort_values('estimate',ascending=False);fig,ax=plt.subplots(figsize=(10.5,9));fig.subplots_adjust(left=.25,right=.75,top=.92,bottom=.11)
 for i,r in enumerate(data.itertuples()):
  if r.flag=='adequate':ci(ax,r.estimate,i,r.lo,r.hi,NAVY)
  else:ax.plot(r.estimate,i,'o',mfc='white',mec=NAVY,ms=7)
  ax.text(1.04,i,f'{r.events:,} / {r.n:,}',transform=ax.get_yaxis_transform(),va='center',size=11)
 ax.set(yticks=range(len(data)),yticklabels=data.group,ylim=(len(data)-.5,-.6),xlim=(0,max(data.hi)*1.05),xlabel='Photography only among complete pairs, % (95% CI)');clean(ax);ax.axvline(val(2025,'All','Photography only').estimate,ls=(0,(4,4)),color=TEAL,lw=1.3);ax.set_title('State-specific incremental contact in 2025',loc='left',pad=22);ax.text(1.04,1.015,'Unweighted n / N',transform=ax.transAxes,weight='bold',size=11);fig.text(.25,.025,'Dashed line: overall 20-state estimate. Hollow markers: fewer than 10 events;\nintervals omitted from the figure. States are descriptive, not ranked on statistical significance.',size=10.5,color=DARK);save(fig,'FigureS2_state_specific_increment')
def fs3():
 fig,axs=plt.subplots(1,2,figsize=(14,6.5));fig.subplots_adjust(left=.09,right=.95,wspace=.78,bottom=.23,top=.85)
 a=axs[0];m=pd.read_csv(P/'results/module_missingness.csv')
 keys=[['Complete eye pair'],['Explicit eye unknown or refusal'],
       ['All module fields blank'],['Both eyes blank with other module data','One eye blank with module data']]
 colors=[NAVY,BLUE,GRAY,TEAL]
 for i,y in enumerate([2025,2022]):
  left=0
  for categories,color in zip(keys,colors):
   estimate=m.loc[m.year.eq(y)&m.category.isin(categories),'weighted_pct'].sum()
   a.barh(i,estimate,left=left,color=color,height=.42,edgecolor='white')
   if estimate>=4:a.text(left+estimate/2,i,f'{estimate:.1f}',ha='center',va='center',
                         color='white' if color==NAVY else '#213441',size=10)
   left+=estimate
 a.set(yticks=[0,1],yticklabels=['2025','2022'],xlim=(0,100),ylim=(1.65,-.65),xlabel='All diagnosed-diabetes respondents, %')
 a.set_title('A  Eye and module response patterns',loc='left',pad=20);a.spines['left'].set_visible(False);a.tick_params(axis='y',length=0)
 handles=[plt.Line2D([0],[0],marker='s',color=c,ms=8,lw=0) for c in colors]
 a.legend(handles,['Complete eye pair','Eye unknown or refused','All module fields blank','Other eye blanks'],
          loc='lower center',bbox_to_anchor=(.5,-.36),frameon=False,ncol=2,fontsize=9.5)
 a=axs[1];results=pd.read_csv(P/'results/nonresponse_sensitivity.csv')
 results=results[results.measure.eq('Photography only')]
 groups=['Complete pairs','IPW','IPW cap 2','IPW + insulin','IPW + insulin cap 2','MI core MAR','MI MAR']
 labels=['Complete pairs','IPW core only','IPW core only, cap 2','IPW with insulin','IPW with insulin, cap 2','MI core only','MI with insulin']
 for i,g in enumerate(groups):
  for y,off,c in [(2025,-.11,NAVY),(2022,.11,BLUE)]:
   r=results[(results.year==y)&results.method.eq(g)].iloc[0];ci(a,r.estimate,i+off,r.lo,r.hi,c)
 a.axhline(.5,color='#CED5DB',lw=1,ls='--')
 a.set(yticks=range(len(groups)),yticklabels=labels,ylim=(len(groups)-.45,-.55),
       xlim=(0,results.loc[results.method.isin(groups),'hi'].max()*1.1),xlabel='Photography-only contact, % (95% CI)')
 clean(a);a.set_title('B  Missing data sensitivity',loc='left',pad=20);a.legend([plt.Line2D([0],[0],color=c,marker='o') for c in [NAVY,BLUE]],['2025','2022'],frameon=False,loc='lower center',bbox_to_anchor=(.5,-.24),ncol=2);save(fig,'FigureS3_missingness_and_sensitivity')
if __name__=='__main__':
 for f in [f1,f2,f3,fs1,fs2,fs3]:f();print(f.__name__,'complete',flush=True)
