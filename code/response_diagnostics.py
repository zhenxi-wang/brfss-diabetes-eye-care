"""Diagnostics for module records, response support and concentration of response weights."""
import numpy as np
import pandas as pd

MODULE_FIELDS = ['DIABTYPE', 'INSULIN1', 'CHKHEMO3', 'EYEEXAM1', 'DIABEYE1', 'DIABEDU1', 'FEETSORE']
EYE_FIELDS = ['EYEEXAM1', 'DIABEYE1']
PATTERN_LABELS = ['Complete eye pair', 'All module fields blank', 'Explicit eye unknown or refusal',
                  'Both eyes blank with other module data', 'One eye blank with module data']
DISPOSITIONS = {1100: 'Completed', 1200: 'Partial'}


def module_diagnostics(sample, year):
    module_blank = sample[MODULE_FIELDS].isna().all(axis=1)
    eyes = sample[EYE_FIELDS]
    explicit = eyes.isin([7, 9]).any(axis=1)
    two_blanks = eyes.isna().all(axis=1)
    group = np.select([sample.pair, module_blank, explicit, two_blanks],
                      PATTERN_LABELS[:4], default=PATTERN_LABELS[4])
    assert eyes.isna().any(axis=1).to_numpy()[group == PATTERN_LABELS[4]].all()
    assert not (sample.pair & module_blank).any()
    w = sample._LLCPWT.to_numpy()
    categories = []
    for label in PATTERN_LABELS:
        mask = group == label
        partial = sample.DISPCODE.eq(1200).to_numpy() & mask
        completed = sample.DISPCODE.eq(1100).to_numpy() & mask
        categories.append(dict(year=year, category=label, n=int(mask.sum()),
                               weighted_pct=100*w[mask].sum()/w.sum(),
                               completed_n=int(completed.sum()), partial_n=int(partial.sum()),
                               other_disposition_n=int((mask & ~partial & ~completed).sum()),
                               partial_weighted_pct=100*w[partial].sum()/w[mask].sum() if mask.any() else np.nan))

    records = sample.assign(pattern=group, disposition=sample.DISPCODE.map(DISPOSITIONS).fillna('Other'),
                            questionnaire=sample.QSTVER)
    routes = []
    for (state, version, disposition), region in records.groupby(
        ['state_name', 'questionnaire', 'disposition'], dropna=False
    ):
        for label, subgroup in region.groupby('pattern'):
            routes.append(dict(year=year, state=state, questionnaire=version, disposition=disposition,
                               category=label, n=len(subgroup), domain_n=len(region),
                               weighted_pct=100*subgroup._LLCPWT.sum()/region._LLCPWT.sum()))

    support = []
    for field, valid_codes, unknown_codes in [
        ('DIABTYPE', [1, 2], [7, 9]),
        ('INSULIN1', [1, 2], [7, 9]),
        ('CHKHEMO3', [*range(1, 77), 88, 98], [77, 99]),
    ]:
        code = sample[field]
        for label, mask in [('Valid', code.isin(valid_codes)), ('Explicit unknown or refusal', code.isin(unknown_codes)),
                            ('Blank', code.isna())]:
            pairs = sample.pair & mask
            support.append(dict(year=year, field=field, level=label, n=int(mask.sum()),
                                complete_n=int(pairs.sum()), weighted_n=w[mask].sum(),
                                weighted_response_pct=100*w[pairs].sum()/w[mask].sum() if mask.any() else np.nan))
    return categories, routes, support


def ranked_weight_influence(paired, factor, year, specification):
    """Remove selected terms from a Hájek ratio with all fitted propensities held fixed."""
    w = paired._LLCPWT.to_numpy()
    event = paired.photo_only.to_numpy()
    ranking = np.argsort(-factor, kind='stable')
    rows = []
    for cap in [np.inf, 2., 1.5]:
        adjusted = w * np.minimum(factor, cap)
        denominator = adjusted.sum()
        numerator = adjusted @ event
        estimate = numerator/denominator
        for n in [1, 5, 10, 20, 50]:
            top = ranking[:n]
            wt = adjusted[top].sum()
            events = adjusted[top] @ event[top]
            omitted = (numerator-events)/(denominator-wt)
            rows.append(dict(year=year, specification=specification, cap='Uncapped' if np.isinf(cap) else str(cap),
                             top_n=n, top_events=int(event[top].sum()),
                             top_insulin_unknown=int(paired.insulin.iloc[top].isna().sum()),
                             smallest_ranked_multiplier=factor[top[-1]], max_multiplier=factor[top[0]],
                             share_total_weight_pct=100*wt/denominator,
                             share_event_weight_pct=100*events/numerator,
                             direct_contribution_pp=100*events/denominator, estimate_pct=100*estimate,
                             without_top_pct=100*omitted, removal_change_pp=100*(omitted-estimate)))
    return rows

