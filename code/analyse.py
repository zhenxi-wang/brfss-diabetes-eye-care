"""Design-based paired ratios and marginal standardized logistic associations."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import t
from scipy.special import expit,logit
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import patsy
P=Path(__file__).resolve().parents[1]
FIPS={2:'Alaska',10:'Delaware',13:'Georgia',17:'Illinois',18:'Indiana',19:'Iowa',23:'Maine',26:'Michigan',28:'Mississippi',29:'Missouri',30:'Montana',32:'Nevada',33:'New Hampshire',34:'New Jersey',35:'New Mexico',36:'New York',37:'North Carolina',38:'North Dakota',39:'Ohio',45:'South Carolina',46:'South Dakota',47:'Tennessee',48:'Texas',50:'Vermont',53:'Washington',55:'Wisconsin',56:'Wyoming'}
class Design:
 def __init__(self,d):
  self.n=len(d);self.w=d._LLCPWT.to_numpy();self.cluster,self.ck=pd.factorize(pd.MultiIndex.from_frame(d[['_STATE','_STSTR','_PSU']]))
  ch=pd.DataFrame({'c':self.cluster,'state':d._STATE,'stratum':d._STSTR}).drop_duplicates('c').sort_values('c')
  self.h,self.hk=pd.factorize(pd.MultiIndex.from_frame(ch[['state','stratum']]))
  self.m=np.bincount(self.h);self.df=len(ch)-len(self.m);self.singleton=int((self.m==1).sum());self.crit=t.ppf(.975,self.df)
 def cov(self,ifv,indices=None):
  a=np.asarray(ifv);a=a[:,None] if a.ndim==1 else a
  cluster=self.cluster if indices is None else self.cluster[indices]
  s=np.zeros((len(self.ck),a.shape[1]));np.add.at(s,cluster,a);hs=np.zeros((len(self.m),a.shape[1]));np.add.at(hs,self.h,s)
  multi=self.m>1;fac=np.where(multi,self.m/np.maximum(self.m-1,1),0.)
  v=s.T@(fac[self.h,None]*s)-hs.T@((fac/self.m)[:,None]*hs)
  if self.singleton:
   u=s[self.m[self.h]==1]-s.mean(axis=0);v+=u.T@u
  return v
 def ratio(self,num,den,scale=100):
  num=np.asarray(num,float);den=np.asarray(den,float);b=self.w@den
  if b<=0:return {'estimate':None,'lo':None,'hi':None,'se':None,'n':0,'events':0,'weighted_n':0,'flag':'empty'}
  est=self.w@num/b;iv=self.w*(num-est*den)/b;se=np.sqrt(max(self.cov(iv)[0,0],0))
  # Logit limits for binary proportions; no degenerate interval claimed at 0/1.
  if 0<est<1 and np.all((num==0)|(num==den)):
   ls=se/(est*(1-est));lo,hi=expit(logit(est)+np.array([-1,1])*self.crit*ls)
  else:lo,hi=est+np.array([-1,1])*self.crit*se
  n=int((den>0).sum());events=int((num>0).sum());flag='limited' if n<50 or events<10 or n-events<10 else 'adequate'
  return {'estimate':est*scale,'lo':lo*scale,'hi':hi*scale,'se':se*scale,'n':n,'events':events,'weighted_n':b,'weighted_events':float(self.w@num),'flag':flag}

def load(year):
 d=pd.read_parquet(P/f'data/derived/brfss{year}_panel.parquet');eligible=pd.read_csv(P/f'sources/eligible_states_{year}.csv').state.tolist();d['state_name']=d._STATE.map(FIPS)
 d=d[d.state_name.isin(eligible)&(d._LLCPWT>0)].copy().reset_index(drop=True)
 d['dm']=d.DIABETE4.eq(1);d['pair']=d.dm & d.EYEEXAM1.isin([1,2,3,4,8]) & d.DIABEYE1.isin([1,2,3,4,8])
 for limit,suffix in [(2,''),(3,'_2y')]:
  d['dil'+suffix]=d.EYEEXAM1.le(limit);d['photo'+suffix]=d.DIABEYE1.le(limit);d['either'+suffix]=d['dil'+suffix]|d['photo'+suffix];d['photo_only'+suffix]=~d['dil'+suffix]&d['photo'+suffix]
 d['age_group']=pd.cut(d._AGE80,[17,44,64,100],labels=['18–44','45–64','≥65']).astype('string')
 d['sex']=d._SEX.map({1:'Male',2:'Female'});race='_RACEGR3' if year==2025 else '_RACEGR4';d['race']=d[race].map({1:'White NH',2:'Black NH',3:'Other/multiracial NH',4:'Other/multiracial NH',5:'Hispanic'})
 d['education']=d.EDUCA.map({1:'<High school',2:'<High school',3:'<High school',4:'High school',5:'Some college',6:'College graduate'})
 d['income']=d.INCOME3.map({**{i:'<$35,000' for i in range(1,6)},6:'$35,000–74,999',7:'$35,000–74,999',**{i:'≥$75,000' for i in range(8,12)}}).fillna('Unknown')
 d['nonmetro']=d._METSTAT.map({1:0.,2:1.});d['cost']=d.MEDCOST1.map({1:1.,2:0.});d['insulin']=d.INSULIN1.map({1:1.,2:0.});d['vision']=d.BLIND.map({1:1.,2:0.});d['type']=d.DIABTYPE.map({1:'Type 1',2:'Type 2'})
 return d

def profile(d,design,mask,label,year,window=''):
 base=d.pair.to_numpy()&np.asarray(mask,bool);a=d['dil'+window].to_numpy();b=d['photo'+window].to_numpy();out=[]
 for name,num,den in [('Dilation',a,base),('Photography',b,base),('Either',a|b,base),('Both',a&b,base),('Dilation only',a&~b,base),('Photography only',~a&b,base),('Neither',~a&~b,base),('Missed among either',~a&b,base&(a|b)),('Photo among no dilation',b,base&~a)]:
  r=design.ratio(np.asarray(den)&num,den);r.update(year=year,group=label,measure=name,window=window or '1y');out.append(r)
 return out

def models(d,design,year,extended=False):
 cols=['age_group','sex','race','education','insulin','nonmetro','cost','state_name']+(['income','vision'] if extended else [])
 mask=d.pair & d[cols].notna().all(axis=1);s=d.loc[mask].copy();s['age_group']=s.age_group.astype(object);ix=s.index.to_numpy();w=s._LLCPWT.to_numpy();wn=w/w.mean();sw=wn.sum()
 formula='1 + nonmetro + cost + C(age_group) + C(sex) + C(race) + C(education) + insulin + C(state_name)'+(' + C(income) + vision' if extended else '')
 x=patsy.dmatrix(formula,s,return_type='dataframe');X=x.to_numpy();info=x.design_info;rows=[];coef=[];ifs={}
 for outcome in ['dil','either','photo_only']:
  y=s[outcome].to_numpy(float);fit=sm.GLM(y,X,family=sm.families.Binomial(),freq_weights=wn).fit(maxiter=100);beta=fit.params;p=fit.fittedvalues;bread=np.linalg.inv(X.T@((wn*p*(1-p))[:,None]*X));score=(wn*(y-p))[:,None]*X
  IFb=np.zeros((len(d),len(beta)));IFb[ix]=score@bread;cov=design.cov(IFb)
  for j,n in enumerate(x.columns):coef.append(dict(year=year,extended=extended,outcome=outcome,term=n,beta=beta[j],se=np.sqrt(max(cov[j,j],0)),n=len(s),converged=fit.converged))
  for exposure in ['nonmetro','cost']:
   z0=s.copy();z1=s.copy();z0[exposure]=0.;z1[exposure]=1.
   X0=np.asarray(patsy.build_design_matrices([info],z0)[0]);X1=np.asarray(patsy.build_design_matrices([info],z1)[0]);p0=expit(X0@beta);p1=expit(X1@beta);g=p1-p0;est=wn@g/sw
   grad=(wn[:,None]*(p1*(1-p1))[:,None]*X1-wn[:,None]*(p0*(1-p0))[:,None]*X0).sum(axis=0)/sw
   iv=np.zeros(len(d));iv[ix]=score@bread@grad+wn*(g-est)/sw;ifs[(outcome,exposure)]=iv
   se=np.sqrt(max(design.cov(iv)[0,0],0));lo,hi=(est+np.array([-1,1])*design.crit*se)*100
   rows.append(dict(year=year,extended=extended,outcome=outcome,exposure=exposure,n=len(s),difference_pp=100*est,lo=lo,hi=hi,se=100*se,p=2*t.sf(abs(est/se),design.df),adjusted_0=100*wn@p0/sw,adjusted_1=100*wn@p1/sw))
  np.savez_compressed(P/f'results/{year}_{"extended" if extended else "base"}_{outcome}_model.npz',beta=beta,cov=cov,columns=x.columns.to_numpy(str),index=ix)
 primary=[r for r in rows if r['outcome'] in ['dil','either']];q=multipletests([r['p'] for r in primary],method='holm')[1]
 for r,v in zip(primary,q):r['p_holm']=v
 for exposure in ['nonmetro','cost']:
  iv=ifs[('either',exposure)]-ifs[('dil',exposure)];a=next(r for r in rows if r['outcome']=='either' and r['exposure']==exposure);b=next(r for r in rows if r['outcome']=='dil' and r['exposure']==exposure);est=a['difference_pp']-b['difference_pp'];se=100*np.sqrt(max(design.cov(iv)[0,0],0));lo,hi=est+np.array([-1,1])*design.crit*se
  rows.append(dict(year=year,extended=extended,outcome='Change in disparity',exposure=exposure,n=len(s),difference_pp=est,lo=lo,hi=hi,se=se,p=2*t.sf(abs(est/se),design.df)))
 pd.DataFrame(coef).to_csv(P/f'results/{year}_{"extended" if extended else "base"}_coefficients.csv',index=False)
 return rows

def main():
 allprof=[];allmodels=[];flows=[];sens=[];miss=[];states=[];datasets={y:load(y) for y in [2025,2022]};common=set(datasets[2025].state_name)&set(datasets[2022].state_name)
 for year,d in datasets.items():
  ds=Design(d);f=dict(year=year,full_records=int(json.loads((P/f'sources/{year}_extraction.json').read_text())['records']),states=d.state_name.nunique(),state_sample=len(d),diabetes=int(d.dm.sum()),complete_pair=int(d.pair.sum()),weighted_complete=ds.w@d.pair.astype(float),psus=len(ds.ck),strata=len(ds.m),design_df=ds.df,singleton_strata=ds.singleton);flows.append(f);print(f,flush=True)
  allprof+=profile(d,ds,np.ones(len(d),bool),'All',year)
  for c in ['nonmetro','age_group','sex','race','education','income','cost','vision','insulin','type']:
   for v in sorted(d.loc[d.dm,c].dropna().unique(),key=str):allprof+=profile(d,ds,d[c].eq(v),f'{c}: {v}',year)
  sens+=profile(d,ds,np.ones(len(d),bool),'Two-year window',year,'_2y');sens+=profile(d,ds,d.type.eq('Type 2'),'Type 2 only',year);sens+=profile(d,ds,d.state_name.isin(common),'Common states',year)
  for st in sorted(d.state_name.unique()):
   states+=profile(d,ds,d.state_name.eq(st),st,year);sens+=profile(d,ds,~d.state_name.eq(st),'Exclude '+st,year)
  for name,c in [('Dilation','EYEEXAM1'),('Photography','DIABEYE1')]:
   for status,m in [('Valid',d[c].isin([1,2,3,4,8])),('Unknown',d[c].eq(7)),('Refused',d[c].eq(9)),('Blank',d[c].isna())]:
    r=ds.ratio(d.dm&m,d.dm);r.update(year=year,item=name,status=status);miss.append(r)
  # Nonparametric bounds on the increment over all diabetes respondents.
  a=d.EYEEXAM1;b=d.DIABEYE1;av=a.isin([1,2,3,4,8]);bv=b.isin([1,2,3,4,8]);positive_a=a.isin([1,2]);positive_b=b.isin([1,2]);lower=d.dm&av&~positive_a&positive_b;upper=d.dm&~positive_a&(~bv|positive_b)
  for name,num in [('Missing-pair lower bound',lower),('Missing-pair upper bound',upper),('Any known positive',d.dm&(positive_a|positive_b))]:
   r=ds.ratio(num,d.dm);r.update(year=year,group=name,measure='Sensitivity',window='1y');sens.append(r)
  known_positive=d.dm&(positive_a|positive_b);evaluable=d.pair|known_positive
  r=ds.ratio(known_positive,evaluable);r.update(year=year,group='Composite evaluable with one positive',measure='Sensitivity',window='1y');sens.append(r)
  allmodels+=models(d,ds,year);allmodels+=models(d,ds,year,True)
  d.to_parquet(P/f'data/derived/{year}_analysis.parquet',index=False)
 pd.DataFrame(allprof).to_csv(P/'results/profiles.csv',index=False);pd.DataFrame(states).to_csv(P/'results/state_profiles.csv',index=False);pd.DataFrame(sens).to_csv(P/'results/sensitivities.csv',index=False);pd.DataFrame(miss).to_csv(P/'results/missingness.csv',index=False);pd.DataFrame(allmodels).to_csv(P/'results/binary_adjusted_contrasts.csv',index=False);(P/'results/flow.json').write_text(json.dumps({'samples':flows,'common_states':sorted(common)},indent=2))
 print(pd.DataFrame(allprof).query('group == "All"')[['year','measure','estimate','lo','hi','n','events']].to_string(index=False),flush=True)
 print(pd.DataFrame(allmodels).query('extended == False')[['year','outcome','exposure','difference_pp','lo','hi','p']].to_string(index=False),flush=True)
if __name__=='__main__':
 main()
 from revision_analysis import main as revision_main
 revision_main()
