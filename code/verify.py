"""Independent checks against raw records, direct weighted counts and stratum sums."""
from pathlib import Path
import json,sys
import pandas as pd,numpy as np
from scipy.special import expit, softmax, logsumexp
from scipy.optimize import root, minimize
import patsy
from analyse import load,Design
P=Path(__file__).resolve().parents[1]

def direct_cov(d,z,lonely=True):
 z=np.atleast_2d(z).T if np.ndim(z)==1 else z;v=np.zeros((z.shape[1],z.shape[1]));mean=z.mean(axis=0)
 assert not d.duplicated(['_STATE','_STSTR','_PSU']).any()
 for _,ii in d.groupby(['_STATE','_STSTR']).indices.items():
  a=z[ii]
  if len(a)>1:b=a-a.mean(axis=0);v+=len(a)/(len(a)-1)*(b.T@b)
  elif lonely:b=a-mean;v+=b.T@b
 return v


def finite_jacobian(function, beta, step=1e-5):
    columns = []
    for direction in np.eye(len(beta)) * step:
        columns.append((function(beta + direction) - function(beta - direction)) / (2*step))
    return np.column_stack(columns)


def verify_joint(d, year):
    """Refit with independent score equations, numeric derivatives and explicit strata."""
    saved = np.load(P / f'results/{year}_joint_M1.npz')
    sample = d[d.pair].dropna(subset=['nonmetro','cost','age_group','sex','race','education','insulin'])
    x = patsy.dmatrix('1 + nonmetro + cost + C(age_group) + C(sex) + C(race) + C(education) + insulin + C(state_name)',
                     sample, return_type='dataframe')
    assert list(x.columns) == list(saved['columns'])
    matrix = x.to_numpy()
    w = sample._LLCPWT.to_numpy()
    w = w / w.mean()
    dil, photo = sample.dil.to_numpy(), sample.photo.to_numpy()
    y = np.column_stack([dil & photo, dil & ~photo, ~dil & photo])
    categories = np.array([[0, 0, 0], [1, 1, 0], [1, 1, 0], [0, 1, 1]])

    def predict(m, beta):
        return softmax(np.column_stack([np.zeros(len(m)), m @ beta.reshape(3, -1).T]), axis=1)

    def contributions(beta):
        residual = y - predict(matrix, beta)[:, 1:]
        return (w[:, None, None] * residual[:, :, None] * matrix[:, None, :]).reshape(len(w), -1)

    def score(beta):
        residual = y - predict(matrix, beta)[:, 1:]
        return ((w[:, None] * residual).T @ matrix).ravel()

    def loss(beta):
        nonreference = matrix @ beta.reshape(3, -1).T
        return np.average(logsumexp(np.column_stack([np.zeros(len(w)), nonreference]), axis=1)
                          - np.sum(y * nonreference, axis=1), weights=w)

    initial = minimize(loss, np.zeros(saved['beta'].size), jac=lambda b: -score(b)/sum(w),
                       method='BFGS', options={'gtol': 1e-8, 'maxiter': 1800})
    answer = root(score, initial.x, tol=1e-9)
    assert np.max(abs(score(answer.x))) < 1e-5, answer.message
    np.testing.assert_allclose(answer.x, saved['beta'].ravel(), atol=1e-6, rtol=0)
    bread = np.linalg.inv(-finite_jacobian(score, answer.x))
    target = pd.read_csv(P / 'results/adjusted_contrasts.csv').query('year == @year and model == "M1"')
    for exposure in ['nonmetro', 'cost']:
        off, on = matrix.copy(), matrix.copy()
        off[:, x.columns.get_loc(exposure)] = 0
        on[:, x.columns.get_loc(exposure)] = 1

        def differences(beta):
            return (predict(on, beta) - predict(off, beta)) @ categories

        def marginal(beta):
            return np.average(differences(beta), axis=0, weights=w)

        est = marginal(answer.x)
        df = len(d) - d.groupby(['_STATE', '_STSTR']).ngroups
        assert df > 0
        iv = np.zeros((len(d), len(est)))
        iv[sample.index] = (contributions(answer.x) @ bread @ finite_jacobian(marginal, answer.x).T
                            + (w / w.sum())[:, None] * (differences(answer.x) - est))
        standard_errors = 100 * np.sqrt(np.diag(direct_cov(d, iv)))
        for i, outcome in enumerate(['dil', 'either', 'Change in disparity']):
            r = target[(target.exposure == exposure) & (target.outcome == outcome)].iloc[0]
            np.testing.assert_allclose(100*est[i], r.difference_pp, atol=1e-7, rtol=0)
            np.testing.assert_allclose(standard_errors[i], r.se, rtol=1e-7)
        np.testing.assert_allclose(iv[:, 1] - iv[:, 0], iv[:, 2], atol=1e-9, rtol=0)


def verify_nonresponse(d, year, include_insulin=False):
    import statsmodels.api as sm
    from revision_analysis import diabetes_sample, auxiliary_design
    sample = diabetes_sample(d)
    complete = sample.pair.to_numpy()
    x = auxiliary_design(sample, include_insulin=include_insulin)
    w = sample._LLCPWT.to_numpy()
    fit = sm.GLM(complete.astype(float), x, freq_weights=w/w.mean(), family=sm.families.Binomial()).fit()
    assert fit.converged
    rows = pd.read_csv(P / 'results/nonresponse_sensitivity.csv').query('year == @year and measure == "Photography only"')
    y = sample.photo_only.to_numpy()[complete]
    prefix = 'IPW + insulin' if include_insulin else 'IPW'
    for cap, method in [(np.inf, prefix), (2, prefix+' cap 2'), (1.5, prefix+' cap 1.5')]:
        def estimate(beta):
            weights = w[complete] * np.minimum(1 / expit(x[complete] @ beta), cap)
            return np.average(y, weights=weights)

        value = estimate(fit.params)
        inverse = np.minimum(1 / fit.predict()[complete], cap)
        weights = w[complete] * inverse
        nuisance = (fit.model.score_obs(fit.params) @ np.linalg.inv(-fit.model.hessian(fit.params))
                    @ finite_jacobian(estimate, fit.params).ravel())
        nuisance[complete] += (y-value)*weights/weights.sum()
        influence = np.zeros(len(d))
        influence[sample.index] = nuisance
        r = rows[rows.method == method].iloc[0]
        np.testing.assert_allclose(100*value, r.estimate, atol=1e-7, rtol=0)
        np.testing.assert_allclose(100*np.sqrt(direct_cov(d, influence)[0, 0]), r.se, atol=1e-6, rtol=0)

def main():
 checks=[];diag=[];rng=np.random.default_rng(20261113)
 for y in [2025,2022]:
  panel=pd.read_parquet(P/f'data/derived/brfss{y}_panel.parquet');selected=set(rng.choice(len(panel),100,replace=False));start=0;n=0
  for a in pd.read_sas(P/f'data/raw/LLCP{y}.xpt',format='xport',chunksize=50000,encoding='utf8'):
   ids=[i for i in selected if start<=i<start+len(a)]
   for i in ids:
    for c in panel.columns:
     aa=a.iloc[i-start][c];bb=panel.iloc[i][c]
     assert (pd.isna(aa) and pd.isna(bb)) or aa==bb,(y,i,c)
    n+=1
   start+=len(a)
  assert n==100;checks.append(f'{y}: 100 random original XPT records matched every extracted field and survey identifier')
  d=load(y);ds=Design(d);pairs=d[d.pair];w=pairs._LLCPWT;profile=pd.read_csv(P/'results/profiles.csv');q=profile[(profile.year==y)&profile.group.eq('All')]
  direct={name:float(np.average(v,weights=w)*100) for name,v in {'Dilation':pairs.dil,'Photography':pairs.photo,'Either':pairs.either,'Both':pairs.dil&pairs.photo,'Dilation only':pairs.dil&~pairs.photo,'Photography only':pairs.photo_only,'Neither':~pairs.either}.items()}
  for name,v in direct.items():assert abs(v-q[q.measure.eq(name)].iloc[0].estimate)<1e-10
  assert abs(direct['Either']-direct['Dilation']-direct['Photography only'])<1e-10
  assert abs(sum(direct[k] for k in ['Both','Dilation only','Photography only','Neither'])-100)<1e-10
  checks.append(f'{y}: all seven paired classification percentages independently matched direct raw weighted counts and partition identities')
  for name,c in [('Dilation','dil'),('Either','either'),('Photography only','photo_only')]:
   mask=d.pair.to_numpy();num=mask&d[c].to_numpy();den=mask;est=ds.w@num/(ds.w@den);iv=ds.w*(num-est*den)/(ds.w@den)
   v=direct_cov(d,iv)[0,0];assert np.isclose(v,ds.cov(iv)[0,0],rtol=1e-10,atol=1e-15);r=q[q.measure.eq(name)].iloc[0];assert np.isclose(100*np.sqrt(v),r.se,rtol=1e-9)
   diag.append(dict(year=y,measure=name,se_center=100*np.sqrt(v),se_omit_singletons=100*np.sqrt(direct_cov(d,iv,False)[0,0]),singleton_strata=ds.singleton))
  checks.append(f'{y}: Taylor variances independently matched explicit within-stratum centered PSU sums; singleton omission sensitivity recorded')
  # Validate reported marginal contrasts using saved coefficients and numerical differentiation.
  cols=['age_group','sex','race','education','insulin','nonmetro','cost','state_name'];s=d[d.pair&d[cols].notna().all(axis=1)].copy();s['age_group']=s.age_group.astype(object)
  X=patsy.dmatrix('1 + nonmetro + cost + C(age_group) + C(sex) + C(race) + C(education) + insulin + C(state_name)',s,return_type='dataframe');W=s._LLCPWT.to_numpy();target=pd.read_csv(P/'results/binary_adjusted_contrasts.csv');target=target[(target.year==y)&~target.extended]
  for outcome in ['dil','either']:
   saved=np.load(P/f'results/{y}_base_{outcome}_model.npz');beta=saved['beta'];assert list(saved['columns'])==list(X.columns);xx=X.to_numpy();p=expit(xx@beta);wn=W/W.mean();score=(wn*(s[outcome].to_numpy()-p))[:,None]*xx;assert np.max(np.abs(score.sum(axis=0)))<1e-5
   exposure='cost';X0=xx.copy();X1=xx.copy();j=list(X.columns).index(exposure);X0[:,j]=0;X1[:,j]=1
   def contrast(b):return np.average(expit(X1@b)-expit(X0@b),weights=W)
   est=contrast(beta);r=target[(target.outcome==outcome)&(target.exposure==exposure)].iloc[0];assert np.isclose(100*est,r.difference_pp)
   h=1e-5;grad=np.array([(contrast(beta+np.eye(len(beta))[j]*h)-contrast(beta-np.eye(len(beta))[j]*h))/(2*h) for j in range(len(beta))]);bread=np.linalg.inv(xx.T@((wn*p*(1-p))[:,None]*xx));g=expit(X1@beta)-expit(X0@beta);iv=np.zeros(len(d));iv[s.index]=score@bread@grad+wn*(g-est)/wn.sum();se=100*np.sqrt(direct_cov(d,iv)[0,0]);assert np.isclose(se,r.se,rtol=2e-7)
  checks.append(f'{y}: binary sensitivity cost-barrier differences and SEs matched finite-difference gradients and independent stratum sums')
  verify_joint(d, y)
  checks.append(f'{y}: independently refitted joint model, both adjusted disparities and gap-change SEs matched using numerical derivatives and explicit strata')
  for include_insulin in [False, True]:
   verify_nonresponse(d, y, include_insulin)
  x=pd.read_csv(P/'results/ipw_ranked_influence.csv').query('year == @y')
  np.testing.assert_allclose(x['without_top_pct'],
                             (x['estimate_pct']-x['direct_contribution_pp'])/(1-x['share_total_weight_pct']/100),
                             atol=1e-9, rtol=0)
  checks.append(f'{y}: independent response models and finite-difference variances matched six core-only and insulin-augmented IPW estimates and SEs; influence ratios reconciled with their numerator and denominator')
 pd.DataFrame(diag).to_csv(P/'results/singleton_sensitivity.csv',index=False)
 (P/'qa/statistical_verification.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2));print('\n'.join(checks))
if __name__=='__main__':main()
