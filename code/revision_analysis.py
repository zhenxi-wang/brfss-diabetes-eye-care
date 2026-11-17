"""Nonresponse, temporal-comparability and model sensitivity analyses for the BRFSS revision."""
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
from scipy.special import expit, logit
from scipy.stats import chi2, t

from analyse import Design, load
from joint_models import (CategoricalFit, MEASURES, PATTERNS, StandardizedModel, add_covariates,
                          adjusted_rows, interval, onehot, outcome_pattern)
from response_diagnostics import MODULE_FIELDS, module_diagnostics, ranked_weight_influence

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
VALID = [1, 2, 3, 4, 8]
IMPUTATIONS = 100
SCENARIOS = {'MI MAR': 1., 'MNAR 0.5x': .5, 'MNAR 1.5x': 1.5, 'MNAR 2x': 2.}
MEASURE_NAMES = ['Dilation', 'Photography', 'Either modality', 'Photography only']
VALUES = np.column_stack(list(MEASURES.values()))
BASE_COMPLETE = ['nonmetro', 'cost', 'age_group', 'sex', 'race', 'education', 'insulin']
CORE_AUXILIARIES = ['state_name', 'nonmetro', 'cost', 'age_group', 'sex', 'race', 'education',
                    'income', 'vision', 'checkup', 'provider', 'health']
SELECTION = {
    'age_group': 'Age', 'sex': 'Sex', 'race': 'Race and ethnicity', 'education': 'Education',
    'income': 'Income', 'nonmetro': 'Nonmetropolitan', 'cost': 'Cost barrier',
    'insulin': 'Insulin use', 'vision': 'Serious vision difficulty', 'checkup': 'Annual routine checkup',
    'provider': 'Personal health care provider',
}


def report(estimate, se, df, measure, method, year, bounded=False, **extra):
    if bounded:
        critical = t.ppf(.975, df)
        low, high = expit(logit(estimate) + np.array([-1, 1]) * critical * se / (estimate * (1 - estimate)))
    else:
        low, high = interval(estimate, se, df)
    return dict(year=year, method=method, measure=measure, estimate=100*estimate, lo=100*low, hi=100*high,
                se=100*se, df=df, **extra)


def auxiliary_design(d, for_imputation=False, include_insulin=False):
    auxiliary = d.copy()
    predictors = CORE_AUXILIARIES.copy()
    if include_insulin:
        # Retain the original column ordering to reproduce the previous fit and draws.
        predictors.insert(predictors.index('income'), 'insulin')
    for column in predictors:
        auxiliary[column] = auxiliary[column].astype('string').fillna('Unknown')
    terms = [f'C({column})' for column in predictors]
    if for_imputation:
        auxiliary['log_weight'] = np.log(d._LLCPWT / d._LLCPWT.median())
        # Pool fine strata via their observed eligible-state diabetes sample size;
        # state, survey weight and this stratum characteristic enter predictions.
        log_size = np.log(d.groupby(['_STATE', '_STSTR'])['_PSU'].transform('size'))
        auxiliary['log_stratum_n'] = log_size - log_size.mean()
        terms.extend(['log_weight', 'log_stratum_n'])
    return np.asarray(patsy.dmatrix('1 + ' + ' + '.join(terms), auxiliary))


def diabetes_sample(d):
    sample = d.loc[d.dm].copy()
    sample['checkup'] = sample.CHECKUP1.map({1: 'Yes', 2: 'No', 3: 'No', 4: 'No', 8: 'No'})
    sample['provider'] = sample.PERSDOC3.map({1: 'Yes', 2: 'Yes', 3: 'No'})
    sample['health'] = sample.GENHLTH.map({1: 'Excellent/very good', 2: 'Excellent/very good',
                                         3: 'Good', 4: 'Fair/poor', 5: 'Fair/poor'})
    return sample


def response_adjustment(sample, design, year, include_insulin=False):
    index = sample.index.to_numpy()
    observed = sample.pair.to_numpy()
    paired = sample.loc[observed]
    weights = sample._LLCPWT.to_numpy()
    x = auxiliary_design(sample, include_insulin=include_insulin)
    fit = CategoricalFit(x, onehot(observed.astype(int), size=2), weights)
    predicted = fit.probabilities(fit.beta)[:, 1]
    factor = 1 / predicted[observed]
    specification = 'Core + insulin' if include_insulin else 'Core only'
    prefix = 'IPW + insulin' if include_insulin else 'IPW'
    dlogw = -(1 - predicted[observed, None]) * x[observed]
    response_influence = fit.score @ fit.bread
    outcomes = outcome_pattern(paired)
    all_rows = []
    diagnostics = []
    for cap in [np.inf, 2., 1.5]:
        method = prefix if np.isinf(cap) else f'{prefix} cap {cap:g}'
        new_weights = weights[observed] * np.minimum(factor, cap)
        sensitivity = dlogw * (factor < cap)[:, None]
        values = VALUES[outcomes]
        mean = np.average(values, axis=0, weights=new_weights)
        ratio_influence = (values - mean) * (new_weights / new_weights.sum())[:, None]
        influence = response_influence @ (sensitivity.T @ ratio_influence)
        influence[observed] += ratio_influence
        se = np.sqrt(np.diag(design.cov(influence, index)))
        for name, est, error in zip(MEASURE_NAMES, mean, se):
            all_rows.append(report(est, error, design.df, name, method, year, bounded=True))
        keep = paired[BASE_COMPLETE].notna().all(axis=1).to_numpy()
        model = StandardizedModel(paired.loc[keep], 'M1')
        model_fit, estimates = model.fit(outcomes[keep], weights=new_weights[keep])
        assert np.array_equal(model.index, paired.index.to_numpy()[keep])
        for exposure, (est, coefficient, target, _, _) in estimates.items():
            direct = coefficient + target
            # Differentiating the weighted estimating equation and the weighted
            # target mean propagates the fitted response model and its covariance.
            nuisance = response_influence @ (sensitivity[keep].T @ direct)
            nuisance[np.flatnonzero(observed)[keep]] += direct
            se = np.sqrt(np.diag(design.cov(nuisance, index)))
            for outcome, i in [('dilation gap', 0), ('either-modality gap', 2), ('gap change', 3)]:
                all_rows.append(report(est[i], se[i], design.df, f'{exposure} {outcome}', method, year))
        diagnostics.append(dict(year=year, method=method, specification=specification, model_n=len(model.index),
                                weighted_response=np.average(observed, weights=weights),
                                probability_min=predicted.min(), probability_p01=np.quantile(predicted, .01),
                                probability_p50=np.median(predicted), probability_max=predicted.max(),
                                multiplier_p99=np.quantile(np.minimum(factor, cap), .99),
                                multiplier_max=np.minimum(factor, cap).max(),
                                capped_n=int((factor > cap).sum()),
                                kish_original=weights[observed].sum()**2/np.sum(weights[observed]**2),
                                kish_adjusted=new_weights.sum()**2/np.sum(new_weights**2),
                                response_information_condition=fit.condition,
                                max_abs_score=max(fit.max_score, model_fit.max_score)))
    return all_rows, diagnostics, ranked_weight_influence(paired, factor, year, specification)


def rubin(estimates, variances, design_df):
    estimates, variances = np.asarray(estimates), np.asarray(variances)
    m = len(estimates)
    mean = estimates.mean(axis=0)
    within = variances.mean(axis=0)
    between = estimates.var(axis=0, ddof=1)
    total = within + (1 + 1/m) * between
    fraction = (1 + 1/m) * between / total
    observed_df = (design_df + 1) * design_df / (design_df + 3) * (1 - fraction)
    df = 1 / (fraction**2 / (m - 1) + 1 / observed_df)
    return mean, np.sqrt(total), df, np.sqrt(between/m), fraction


def imputation(sample, design, year, include_insulin=True):
    rng = np.random.default_rng(20261116 + year)
    weights = sample._LLCPWT.to_numpy()
    index = sample.index.to_numpy()
    responses = sample[['EYEEXAM1', 'DIABEYE1']]
    valid = responses.isin(VALID).to_numpy()
    positive = responses.isin([1, 2]).to_numpy()
    compatible = np.all(~valid[:, None, :] | (positive[:, None, :] == PATTERNS[None, :, :]), axis=2)
    x = auxiliary_design(sample, for_imputation=True, include_insulin=include_insulin)
    # Sparse unknown-education strata quasi-separate one 2022 response pattern.
    # A weak Gaussian penalty stabilizes the imputation fit in both years.
    fit = CategoricalFit(x, compatible, weights, ridge=1 / 2.5**2)
    covariance = design.cov(fit.score @ fit.bread, index) + fit.bread @ np.diag(fit.penalty) @ fit.bread
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if eigenvalues.min() < -1e-9:
        raise RuntimeError('Imputation model covariance is not positive semidefinite')
    root = eigenvectors * np.sqrt(eigenvalues.clip(min=0))
    model = StandardizedModel(sample.dropna(subset=BASE_COMPLETE), 'M1')
    model_rows = sample.index.get_indexer(model.index)
    labels = [*MEASURE_NAMES, 'nonmetro dilation gap', 'nonmetro either-modality gap', 'nonmetro gap change',
              'cost dilation gap', 'cost either-modality gap', 'cost gap change']
    scenarios = SCENARIOS if include_insulin else {'MI core MAR': 1.}
    combined = {method: ([], []) for method in scenarios}
    start = dict.fromkeys(scenarios)
    module_blank = sample[MODULE_FIELDS].isna().all(axis=1).to_numpy()
    domain_rates = []
    maximum_score = fit.max_score
    maximum_condition = fit.condition
    for replicate in range(IMPUTATIONS):
        beta = fit.beta + root @ rng.standard_normal(len(fit.beta))
        _, probability = fit.probabilities(beta, observed=True)
        uniform = rng.random(len(sample))
        for method, multiplier in scenarios.items():
            adjusted = probability.copy()
            # Tilt only records with missing photography. Respect observed dilation
            # and preserve the relative probabilities of other compatible patterns.
            eligible = ~valid[:, 1] & compatible[:, 3]
            mass = probability[eligible, 3]
            target = np.minimum(multiplier * mass, 1)
            adjusted[eligible] *= ((1 - target)/(1 - mass))[:, None]
            adjusted[eligible, 3] = target
            np.testing.assert_allclose(adjusted.sum(axis=1), 1, atol=1e-10)
            outcome = np.sum(uniform[:, None] > adjusted.cumsum(axis=1)[:, :-1], axis=1)
            assert np.all(compatible[np.arange(len(sample)), outcome])
            values = VALUES[outcome]
            if multiplier == 1:
                domain_rates.append(np.average(values[module_blank, 3], weights=weights[module_blank]))
            mean = np.average(values, axis=0, weights=weights)
            iv = (values - mean) * (weights/weights.sum())[:, None]
            variance = list(np.diag(design.cov(iv, index)))
            means = list(mean)
            new_fit, estimates = model.fit(outcome[model_rows], start=start[method])
            start[method] = new_fit.beta
            maximum_score = max(maximum_score, new_fit.max_score)
            maximum_condition = max(maximum_condition, new_fit.condition)
            for est, coefficient, target_term, _, _ in estimates.values():
                columns = [0, 2, 3]
                means.extend(est[columns])
                variance.extend(np.diag(design.cov(coefficient + target_term, model.index))[columns])
            combined[method][0].append(means)
            combined[method][1].append(variance)
        if (replicate + 1) % 25 == 0:
            print(f'{year}: {replicate + 1}/{IMPUTATIONS} imputations, {len(scenarios)} scenario(s)', flush=True)
    rows = []
    for method, (estimates, variances) in combined.items():
        mean, se, df, mc, fraction = rubin(estimates, variances, design.df)
        half = IMPUTATIONS // 2
        half_difference = abs(np.mean(estimates[:half], axis=0) - np.mean(estimates[half:], axis=0))
        for i, name in enumerate(labels):
            rows.append(report(mean[i], se[i], df[i], name, method, year, bounded=i < len(MEASURE_NAMES),
                               imputations=IMPUTATIONS, mcse=100*mc[i], missing_information=fraction[i],
                               half_difference=100*half_difference[i], model_n=len(model.index)))
    return rows, dict(year=year, specification='Core + insulin' if include_insulin else 'Core only',
                      model_n=len(model.index), imputation_width=len(fit.beta),
                      max_abs_score=maximum_score, max_information_condition=maximum_condition,
                      all_module_blank_weighted_pct=100*weights[module_blank].sum()/weights.sum(),
                      all_module_blank_imputed_pct=100*np.mean(domain_rates),
                      all_module_blank_imputed_contribution_pp=100*np.mean(domain_rates)*weights[module_blank].sum()/weights.sum())


def observed_nonresponse(sample, design, year):
    valid_dil = sample.EYEEXAM1.isin(VALID)
    valid_photo = sample.DIABEYE1.isin(VALID)
    domain = np.ones(len(sample), bool)
    weights = sample._LLCPWT.to_numpy()
    patterns = {'Both valid': valid_dil & valid_photo, 'Dilation valid only': valid_dil & ~valid_photo,
                'Photography valid only': ~valid_dil & valid_photo, 'Neither valid': ~valid_dil & ~valid_photo}
    counts = []
    for name, mask in patterns.items():
        counts.append(dict(year=year, pattern=name, n=int(mask.sum()), N=len(sample),
                           weighted_pct=100*np.average(mask, weights=weights)))
    characteristics = []
    paired = sample.pair.to_numpy()
    for variable, title in SELECTION.items():
        series = sample[variable].replace({0.: 'No', 1.: 'Yes'}).astype('string').fillna('Unknown')
        for label in sorted(series.unique()):
            mask = series.eq(label).to_numpy()
            first = np.average(mask[paired], weights=weights[paired])
            second = np.average(mask[~paired], weights=weights[~paired])
            sd = np.sqrt((first*(1-first) + second*(1-second))/2)
            characteristics.append(dict(year=year, characteristic=title, level=label, complete_pct=100*first,
                                        incomplete_pct=100*second, smd=(first-second)/sd if sd else 0.))
    state_rows = []
    for state, group in sample.groupby('state_name'):
        sw = group._LLCPWT.to_numpy()
        state_rows.append(dict(year=year, state=state, diabetes_n=len(group), complete_n=int(group.pair.sum()),
                               complete_pct=100*np.average(group.pair, weights=sw),
                               photo_unknown_pct=100*np.average(group.DIABEYE1.eq(7), weights=sw),
                               dil_unknown_pct=100*np.average(group.EYEEXAM1.eq(7), weights=sw)))
    model_sample = sample.loc[sample.pair].dropna(subset=BASE_COMPLETE)
    overlap = []
    x = patsy.dmatrix('nonmetro + C(state_name) + C(age_group) + C(sex) + C(race) + C(education) + insulin',
                      model_sample, return_type='dataframe').to_numpy()
    weights = model_sample._LLCPWT.to_numpy()
    fit = CategoricalFit(x, onehot(model_sample.cost.to_numpy(int), size=2), weights)
    model_sample = model_sample.assign(cost_propensity=fit.probabilities(fit.beta)[:, 1])
    for col in ['state_name', 'age_group', 'sex', 'race', 'education', 'nonmetro']:
        for label, g in model_sample.groupby(col):
            for cost, group in g.groupby('cost'):
                q = np.quantile(group.cost_propensity, [.01, .5, .99])
                overlap.append(dict(year=year, characteristic=col, level=label, cost=int(cost), n=len(group),
                                    weighted_pct=100*group._LLCPWT.sum()/g._LLCPWT.sum(),
                                    propensity_p01=q[0], propensity_p50=q[1], propensity_p99=q[2]))
    return counts, characteristics, state_rows, overlap


def temporal_comparison(datasets):
    common = sorted(set.intersection(*(set(d.state_name) for d in datasets.values())))
    years = [2022, 2025]
    designs = {}
    parts = []
    crude = []
    for year in years:
        d = datasets[year].loc[lambda x: x.state_name.isin(common)].reset_index(drop=True)
        designs[year] = Design(d)
        r = designs[year].ratio(d.pair & d.photo_only, d.pair)
        crude.append(r)
        parts.append(d.loc[d.pair].assign(year=year, design_index=lambda x: x.index))
    s = pd.concat(parts, ignore_index=True)
    x = patsy.dmatrix('C(year) + C(state_name) + C(age_group) + C(sex) + C(race) + C(education)',
                      s, return_type='dataframe')
    s = s.loc[x.index]
    weights = s._LLCPWT.to_numpy()
    y = s.year.to_numpy()
    share = np.zeros(len(s))
    for year in years:
        share[y == year] = .5 * weights[y == year]/weights[y == year].sum()
    fit = CategoricalFit(x.to_numpy(), onehot(s.photo_only.to_numpy(int), size=2), weights)
    values, gradients = [], []
    for year in years:
        target = np.asarray(patsy.build_design_matrices([x.design_info], s.assign(year=year))[0])
        value, gradient = fit.marginal(target, share, np.array([[0], [1]]))
        values.append(value.ravel())
        gradients.append(gradient.ravel())
    values = np.column_stack([*values, values[1] - values[0]])
    means = share @ values
    gradients = np.column_stack([*gradients, gradients[1] - gradients[0]])
    effect = fit.score @ (fit.bread @ gradients)
    covariances = []
    for year in years:
        mask = y == year
        center = np.average(values[mask], axis=0, weights=weights[mask])
        influence = effect[mask] + share[mask, None] * (values[mask] - center)
        covariances.append(designs[year].cov(influence, s.loc[mask, 'design_index'].to_numpy()))
    v = np.diag(covariances[0] + covariances[1])
    # Satterthwaite df over independent year-specific variance contributions.
    df = v**2 / sum(np.diag(c)**2/designs[y].df for c, y in zip(covariances, years))
    out = []
    for i, year in enumerate(years):
        r = crude[i]
        out.append(dict(year=year, method='Crude', measure='Prevalence', estimate=r['estimate'], lo=r['lo'],
                        hi=r['hi'], se=r['se'], df=designs[year].df, n=r['n']))
        out.append(report(means[i], np.sqrt(v[i]), df[i], 'Prevalence', 'Composition standardized', year,
                          bounded=True, n=int((y == year).sum())))
    v0, v1 = [r['se']**2/10000 for r in crude]
    df_crude = (v0 + v1)**2 / (v0**2/designs[2022].df + v1**2/designs[2025].df)
    out.extend([
        report((crude[1]['estimate']-crude[0]['estimate'])/100, np.sqrt(v0+v1), df_crude, '2025 minus 2022',
               'Crude', 'Both', n=sum(r['n'] for r in crude)),
        report(means[2], np.sqrt(v[2]), df[2], '2025 minus 2022', 'Composition standardized', 'Both', n=len(s)),
    ])
    return out


def state_heterogeneity():
    """Compare observed state variation with within-state sampling error."""
    estimates = pd.read_csv(RESULTS / 'state_profiles.csv').query('measure == "Photography only"')
    rows = []
    for year, data in estimates.groupby('year', sort=False):
        for domain, states in [('All states', data), ('Adequate counts', data[data.flag == 'adequate'])]:
            raw = states.estimate.to_numpy() / 100
            se = states.se.to_numpy() / (100 * raw * (1 - raw))
            weights = 1 / se**2
            log_odds = logit(raw)
            center = np.average(log_odds, weights=weights)
            statistic = np.dot(weights, (log_odds - center)**2)
            df = len(states) - 1
            rows.append(dict(year=year, domain=domain, states=len(states), minimum_pct=100*raw.min(),
                             maximum_pct=100*raw.max(), Q=statistic, df=df, p=chi2.sf(statistic, df),
                             I2_pct=100*max(0, (statistic - df)/statistic)))
    pd.DataFrame(rows).to_csv(RESULTS / 'state_heterogeneity.csv', index=False)


def main():
    datasets = {year: load(year) for year in (2025, 2022)}
    all_tables = {stem: [] for stem in ['adjusted_contrasts', 'joint_diagnostics', 'nonresponse_sensitivity',
                                       'ipw_diagnostics', 'mi_diagnostics', 'response_patterns',
                                       'selection_characteristics', 'state_response', 'cost_overlap',
                                       'module_missingness', 'module_route_patterns', 'module_response_support',
                                       'ipw_ranked_influence']}
    for year, d in datasets.items():
        design = Design(d)
        sample = diabetes_sample(d)
        rows, diagnostics = adjusted_rows(d, design, year, RESULTS)
        all_tables['adjusted_contrasts'].extend(rows)
        all_tables['joint_diagnostics'].extend(diagnostics)
        for name, values in zip(['response_patterns', 'selection_characteristics', 'state_response', 'cost_overlap'],
                                observed_nonresponse(sample, design, year)):
            all_tables[name].extend(values)
        for column, name in zip(VALUES.T, MEASURE_NAMES):
            ix = sample.index.to_numpy()
            outcome = outcome_pattern(sample)
            pair = sample.pair.to_numpy()
            denominator = d.pair.to_numpy()
            numerator = np.zeros(len(d), bool)
            numerator[ix] = pair & column[outcome].astype(bool)
            r = design.ratio(numerator, denominator)
            all_tables['nonresponse_sensitivity'].append(dict(
                year=year, method='Complete pairs', measure=name, estimate=r['estimate'],
                lo=r['lo'], hi=r['hi'], se=r['se'], df=design.df))
        for row in rows:
            if row['model'] != 'M1' or row['exposure'] == 'hispanic' or row['outcome'] == 'photo_only':
                continue
            label = {'dil': 'dilation gap', 'either': 'either-modality gap', 'Change in disparity': 'gap change'}
            all_tables['nonresponse_sensitivity'].append(dict(
                year=year, method='Complete pairs', measure=f"{row['exposure']} {label[row['outcome']]}",
                estimate=row['difference_pp'], lo=row['lo'], hi=row['hi'], se=row['se'], df=design.df))
        for name, values in zip(['module_missingness', 'module_route_patterns', 'module_response_support'],
                                module_diagnostics(sample, year)):
            all_tables[name].extend(values)
        for include_insulin in [False, True]:
            ipw, ipw_diagnostics, ranked = response_adjustment(sample, design, year, include_insulin)
            all_tables['nonresponse_sensitivity'].extend(ipw)
            all_tables['ipw_diagnostics'].extend(ipw_diagnostics)
            all_tables['ipw_ranked_influence'].extend(ranked)
            mi, mi_diagnostics = imputation(sample, design, year, include_insulin)
            all_tables['nonresponse_sensitivity'].extend(mi)
            all_tables['mi_diagnostics'].append(mi_diagnostics)
    all_tables['common_state_temporal'] = temporal_comparison(datasets)
    for filename, rows in all_tables.items():
        pd.DataFrame(rows).to_csv(RESULTS / f'{filename}.csv', index=False)
    state_heterogeneity()


if __name__ == '__main__':
    main()

