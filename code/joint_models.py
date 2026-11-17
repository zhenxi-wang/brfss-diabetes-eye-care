"""Survey-weighted categorical regressions and marginal-standardization influence functions."""
import numpy as np
import pandas as pd
import patsy
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

PATTERNS = np.array([[0, 0], [1, 1], [1, 0], [0, 1]])
MEASURES = {
    'dil': PATTERNS[:, 0],
    'photo': PATTERNS[:, 1],
    'either': np.any(PATTERNS, axis=1).astype(float),
    'photo_only': (PATTERNS[:, 0] == 0) & (PATTERNS[:, 1] == 1),
}
BASE = 'nonmetro + cost + C(age_group) + C(sex) + C(race) + C(education) + insulin + C(state_name)'
FORMULAS = {
    'M0': 'nonmetro + cost + C(state_name)',
    'M1': BASE,
    'M2': BASE + ' + C(income) + vision',
    'M3': BASE.replace(' + insulin', ''),
    'M4': BASE + ' + C(duration) + C(type_proxy)',
}
EXPOSURES = {
    'nonmetro': ('nonmetro', 0., 1.),
    'cost': ('cost', 0., 1.),
    'hispanic': ('race', 'White NH', 'Hispanic'),
}


def add_covariates(d):
    """Harmonized auxiliaries; diagnosis duration is approximate (age minus reported age at diagnosis)."""
    d = d.copy()
    duration = (d._AGE80 - d.DIABAGE4).where(
        d.DIABAGE4.between(1, 97) & d._AGE80.lt(80) & d._AGE80.ge(d.DIABAGE4))
    d['duration'] = pd.cut(duration, [-1, 4, 9, 100], labels=['0–4', '5–9', '10+']).astype('string').fillna('Unknown')
    d['type_proxy'] = d['type'].fillna('Unknown')
    return d


class CategoricalFit:
    """Reference-category logit fit, optionally using partially observed outcome patterns."""
    def __init__(self, x, compatible, weights, start=None, ridge=0):
        self.x = x
        self.compatible = compatible
        self.weights = weights / np.mean(weights)
        self.k = compatible.shape[1] - 1
        self.width = x.shape[1]
        self.penalty = np.full(self.k * self.width, ridge, dtype=float)
        self.penalty[::self.width] = 0
        if start is None:
            start = np.zeros(self.k * self.width)
        solved = minimize(self.objective, start, jac=True, hess=self.hessian,
                          method='trust-exact', options={'gtol': 1e-7, 'maxiter': 100})
        beta = solved.x
        _, gradient = self.objective(beta)
        if 1e-7 < np.max(np.abs(gradient)) < 1e-3:
            # Near the optimum the objective can round to an unchanged value.
            # A Newton step can still reduce the score at this precision.
            beta = beta - np.linalg.solve(self.hessian(beta), gradient)
            _, gradient = self.objective(beta)
        if np.max(np.abs(gradient)) > 2e-5:
            raise RuntimeError(f'Categorical model failed to converge: {solved.message}, score={np.max(np.abs(gradient)):.3g}')
        self.beta = beta
        self.max_score = np.max(np.abs(gradient))
        information = self.hessian(self.beta)
        self.condition = np.linalg.cond(information)
        if np.linalg.eigvalsh(information)[0] <= 0:
            raise RuntimeError('Categorical likelihood information is not positive definite')
        self.bread = np.linalg.inv(information)
        probability, observed = self.probabilities(self.beta, observed=True)
        self.score = ((observed[:, 1:] - probability[:, 1:]) * self.weights[:, None])[:, :, None] * x[:, None, :]
        self.score = self.score.reshape(len(x), -1)

    def probabilities(self, beta, x=None, observed=False):
        x = self.x if x is None else x
        logits = np.column_stack([np.zeros(len(x)), x @ beta.reshape(self.k, self.width).T])
        total = logsumexp(logits, axis=1, keepdims=True)
        probability = np.exp(logits - total)
        if observed:
            restricted = np.where(self.compatible, logits, -np.inf)
            observed_total = logsumexp(restricted, axis=1, keepdims=True)
            return probability, np.exp(restricted - observed_total)
        return probability

    def objective(self, beta):
        logits = np.column_stack([np.zeros(len(self.x)), self.x @ beta.reshape(self.k, self.width).T])
        total = logsumexp(logits, axis=1, keepdims=True)
        restricted = np.where(self.compatible, logits, -np.inf)
        observed_total = logsumexp(restricted, axis=1, keepdims=True)
        difference = np.exp(logits - total) - np.exp(restricted - observed_total)
        gradient = (difference[:, 1:] * self.weights[:, None]).T @ self.x
        loss = np.dot(self.weights, (total - observed_total).ravel())
        return loss + .5 * np.dot(self.penalty, beta**2), gradient.ravel() + self.penalty * beta

    def hessian(self, beta):
        probability, observed = self.probabilities(beta, observed=True)
        p, q = probability[:, 1:], observed[:, 1:]
        out = np.empty((self.k, self.width, self.k, self.width))
        for i in range(self.k):
            for j in range(i + 1):
                curvature = (p[:, i] - q[:, i]) * (i == j) - p[:, i] * p[:, j] + q[:, i] * q[:, j]
                block = self.x.T @ ((self.weights * curvature)[:, None] * self.x)
                out[i, :, j, :] = block
                out[j, :, i, :] = block
        information = out.reshape(self.k * self.width, -1)
        return information + np.diag(self.penalty)

    def marginal(self, x, target_weights, measures):
        probability = self.probabilities(self.beta, x)
        values = probability @ measures
        share = target_weights / target_weights.sum()
        gradient = np.concatenate([
            x.T @ ((share * probability[:, j])[:, None] * (measures[j] - values))
            for j in range(1, self.k + 1)
        ], axis=0)
        return values, gradient


def interval(estimate, se, df):
    critical = t.ppf(.975, df)
    return estimate - critical * se, estimate + critical * se


def onehot(outcome, size=4):
    return outcome[:, None] == np.arange(size)


def outcome_pattern(d):
    return np.select([d.dil & d.photo, d.dil, d.photo], [1, 2, 3], default=0)


class StandardizedModel:
    """Common fit/sample for all outcomes and contrasts in an adjustment specification."""
    def __init__(self, sample, specification, exposures=('nonmetro', 'cost')):
        x = patsy.dmatrix(FORMULAS[specification], sample, return_type='dataframe')
        self.index = x.index.to_numpy()
        self.x = x.to_numpy()
        self.columns = x.columns.to_numpy(str)
        self.original = sample.loc[self.index, '_LLCPWT'].to_numpy()
        self.prediction = {}
        for exposure in exposures:
            variable, first, second = EXPOSURES[exposure]
            self.prediction[exposure] = []
            for value in (first, second):
                target = sample.loc[self.index].assign(**{variable: value})
                self.prediction[exposure].append(np.asarray(patsy.build_design_matrices([x.design_info], target)[0]))

    def fit(self, outcome, weights=None, target_weights=None, start=None):
        weights = self.original if weights is None else weights
        target_weights = weights if target_weights is None else target_weights
        fit = CategoricalFit(self.x, onehot(outcome), weights, start)
        measures = np.column_stack(list(MEASURES.values()))
        result = {}
        for exposure, (x0, x1) in self.prediction.items():
            values0, gradient0 = fit.marginal(x0, target_weights, measures)
            values1, gradient1 = fit.marginal(x1, target_weights, measures)
            contrast = values1 - values0
            share = target_weights / target_weights.sum()
            adjusted0 = share @ values0
            adjusted1 = share @ values1
            estimate = share @ contrast
            np.testing.assert_allclose(estimate[2] - estimate[0], estimate[3], atol=1e-8, rtol=0)
            coefficient = fit.score @ (fit.bread @ (gradient1 - gradient0))
            target = share[:, None] * (contrast - estimate)
            result[exposure] = (estimate, coefficient, target, adjusted0, adjusted1)
        return fit, result


def adjusted_rows(d, design, year, directory):
    rows = []
    diagnostics = []
    paired = add_covariates(d.loc[d.pair])
    categories = pd.Series(outcome_pattern(paired), index=paired.index)
    for specification in FORMULAS:
        exposures = ('nonmetro', 'cost', 'hispanic') if specification == 'M1' else ('nonmetro', 'cost')
        # Keep the M0/M3/M4 comparisons on M1's covariate-complete sample.
        # M0 still adjusts only for the focal exposure and state.
        sample = paired.dropna(subset=['nonmetro', 'cost', 'age_group', 'sex', 'race', 'education', 'insulin'])
        models = {e: StandardizedModel(sample, specification, (e,)) for e in exposures[:2]} if specification == 'M0' else {}
        if specification == 'M0':
            for exposure, model in models.items():
                other = 'nonmetro' if exposure == 'cost' else 'cost'
                drop = list(model.columns).index(other)
                model.x = np.delete(model.x, drop, axis=1)
                model.columns = np.delete(model.columns, drop)
                model.prediction[exposure] = [np.delete(x, drop, axis=1) for x in model.prediction[exposure]]
        else:
            models[specification] = StandardizedModel(sample, specification, exposures)
        for label, model in models.items():
            fit, results = model.fit(categories.loc[model.index].to_numpy())
            covariance = design.cov(fit.score @ fit.bread, model.index)
            np.savez_compressed(directory / f'{year}_joint_{label if specification != "M0" else "M0_" + label}.npz',
                                beta=fit.beta.reshape(3, -1), cov=covariance, columns=model.columns,
                                categories=['Both', 'Dilation only', 'Photography only'])
            diagnostics.append(dict(year=year, model=specification, fit=label, n=len(model.index),
                                    max_abs_score=fit.max_score, information_condition=fit.condition))
            for exposure, (estimate, coefficient, target, adjusted0, adjusted1) in results.items():
                influence = coefficient + target
                np.testing.assert_allclose(influence[:, 2] - influence[:, 0], influence[:, 3], atol=1e-10, rtol=0)
                covariance = design.cov(influence, model.index)
                for i, outcome in enumerate(MEASURES):
                    if outcome == 'photo' or (exposure == 'hispanic' and outcome != 'photo_only'):
                        continue
                    est = estimate[i] * 100
                    se = np.sqrt(covariance[i, i]) * 100
                    lo, hi = interval(est, se, design.df)
                    row = dict(year=year, model=specification, extended=specification == 'M2',
                               outcome=outcome, exposure=exposure, n=len(model.index),
                               difference_pp=est, lo=lo, hi=hi, se=se, p=2*t.sf(abs(est/se), design.df),
                               adjusted_0=adjusted0[i]*100, adjusted_1=adjusted1[i]*100)
                    rows.append(row)
                    if outcome == 'photo_only' and exposure != 'hispanic':
                        rows.append(dict(row, outcome='Change in disparity'))
    for specification in ('M1', 'M2'):
        primary = [r for r in rows if r['model'] == specification and r['outcome'] in ('dil', 'either')]
        for row, p in zip(primary, multipletests([r['p'] for r in primary], method='holm')[1]):
            row['p_holm'] = p
    return rows, diagnostics

