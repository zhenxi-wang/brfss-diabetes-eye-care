"""Generate publication tables for analyses added after review."""
import pandas as pd


def revision_tables(results, table, ci):
    nonresponse = pd.read_csv(results / 'nonresponse_sensitivity.csv')
    rows = []
    measures = ['Dilation', 'Photography', 'Either modality', 'Photography only',
                'nonmetro gap change', 'cost gap change']
    labels = [*measures[:4], 'Geography gap change', 'Cost gap change']
    for year in (2025, 2022):
        for method in ['Complete pairs', 'IPW', 'IPW cap 2', 'IPW cap 1.5', 'IPW + insulin',
                       'IPW + insulin cap 2', 'IPW + insulin cap 1.5', 'MI core MAR', 'MI MAR',
                       'MNAR 0.5x', 'MNAR 1.5x', 'MNAR 2x']:
            selected = nonresponse[(nonresponse.year == year) & (nonresponse.method == method)]
            show = [ci(selected[selected.measure == name].iloc[0], 2) for name in measures]
            rows.append([year, method, *show])
    tables = [table('Supplementary Table S11 Eye care estimates under alternative missing data assumptions',
                    ['Year', 'Method', *labels], rows, [.43, .85, .83, .83, .83, .88, .925, .925],
                    'Percentages and 95% CIs for contact. Gap changes are in percentage points. Complete-pair percentages use respondents with both outcomes observed. Other prevalence estimates target all eligible respondents with diagnosed diabetes. Adjusted gap analyses require observed base model covariates. IPW uses core predictors unless marked + insulin. Caps apply to inverse response multipliers. MI MAR and MNAR use the imputation model including insulin. Full contrasts and Monte Carlo diagnostics are in nonresponse_sensitivity.csv.')]

    temporal = pd.read_csv(results / 'common_state_temporal.csv')
    rows = [[r.method, '2025 minus 2022' if r.year == 'Both' else r.year, f'{r.n:,}', ci(r, 2)]
            for r in temporal.itertuples()]
    tables.append(table('Supplementary Table S12 Common state temporal comparison',
                        ['Method', 'Year or contrast', 'n', 'Estimate (95% CI)'], rows, [1.95, 1.6, .8, 2.15],
                        'Photography-only percentage for each year and difference in percentage points. The design and weights remain year specific. The standardized estimates use an equally weighted mixture of the two years’ weighted covariate distributions across Delaware, Georgia, Indiana, Maine, Missouri, North Dakota and Wisconsin. Covariates are state, age, sex, race and ethnicity, and education. All analyses use complete eye-care pairs.'))

    selection = pd.read_csv(results / 'selection_characteristics.csv')
    left = selection[selection.year == 2025].set_index(['characteristic', 'level'])
    right = selection[selection.year == 2022].set_index(['characteristic', 'level'])
    rows = []
    for key in left.index.union(right.index, sort=False):
        values = []
        for column in (left, right):
            if key not in column.index:
                values.extend(['—'] * 3)
                continue
            r = column.loc[key]
            values.extend([f'{r.complete_pct:.1f}', f'{r.incomplete_pct:.1f}', f'{r.smd:+.2f}'])
        rows.append([f'{key[0]} | {key[1]}', *values])
    tables.append(table('Supplementary Table S13 Characteristics by complete pair status',
                        ['Characteristic', '2025 complete %', '2025 incomplete %', '2025 SMD',
                         '2022 complete %', '2022 incomplete %', '2022 SMD'],
                        rows, [2.09, .76, .76, .685, .76, .76, .685],
                        'Weighted distributions within complete and incomplete groups. SMD is the signed standardized mean difference for each indicator (complete minus incomplete, divided by the square root of the average of the two Bernoulli variances). Unknown values are shown when present.'))

    states = pd.read_csv(results / 'state_response.csv')
    heterogeneity = pd.read_csv(results / 'state_heterogeneity.csv').query('domain == "All states"')
    summaries = ' '.join(f'In {r.year}, state estimates ranged from {r.minimum_pct:.2f}% to {r.maximum_pct:.2f}% '
                         f'(Q={r.Q:.1f}, df={r.df}, I²={r.I2_pct:.0f}%).' for r in heterogeneity.itertuples())
    rows = [[r.year, r.state, f'{r.complete_n:,}/{r.diabetes_n:,}', f'{r.complete_pct:.1f}',
             f'{r.photo_unknown_pct:.1f}', f'{r.dil_unknown_pct:.1f}'] for r in states.itertuples()]
    tables.append(table('Supplementary Table S14 State specific eye item response',
                        ['Year', 'State', 'Pairs / diabetes', 'Pairs %', 'Photo unknown %', 'Dilation unknown %'],
                        rows, [.55, 1.5, 1.25, .75, 1.2, 1.25],
                        'Response percentages use all diagnosed-diabetes respondents within each state. Unknown denotes an explicit do-not-know response. Within-year heterogeneity in complete-pair photography-only estimates used inverse-variance Q on the logit scale and I² as descriptive measures. '
                        + summaries + ' Results excluding sparse state estimates are available in state_heterogeneity.csv.'))

    diag = pd.read_csv(results / 'ipw_diagnostics.csv')
    rows = [[r.year, r.method, f'{r.probability_min:.3f} / {r.probability_p01:.3f}', f'{r.probability_max:.3f}',
             f'{r.multiplier_p99:.2f} / {r.multiplier_max:.2f}', f'{r.capped_n:,}',
             f'{r.kish_original:,.0f} / {r.kish_adjusted:,.0f}'] for r in diag.itertuples()]
    tables.append(table('Supplementary Table S15 Response probabilities and weight diagnostics',
                        ['Year', 'Method', 'Probability min / P1', 'Probability max', 'Multiplier P99 / max',
                         'Capped n', 'Kish ESS original / adjusted'],
                        rows, [.5, .87, 1.1, .92, 1.06, .65, 1.4],
                        'Predicted probabilities are evaluated among all eligible diagnosed-diabetes respondents. Multipliers and Kish effective sample sizes (ESS) describe complete pairs. IPW uses core predictors unless marked + insulin. Pooling explicit unknown answers and blanks in the insulin sensitivity creates apparent modelled overlap. The separate blank levels of each earlier module item had zero complete pairs (Table S19).'))

    adj = pd.read_csv(results / 'adjusted_contrasts.csv')
    binary = pd.read_csv(results / 'binary_adjusted_contrasts.csv')
    baseline = adj[adj.model.isin(['M0', 'M3', 'M4'])]
    rows = []
    for year in (2025, 2022):
        for exposure, label in [('nonmetro', 'Geography'), ('cost', 'Cost')]:
            for model in ('M0', 'M3', 'M4', 'Binary base', 'Binary extended'):
                if model.startswith('Binary'):
                    group = binary[(binary.year == year) & (binary.exposure == exposure)
                                   & (binary.extended == (model == 'Binary extended'))]
                else:
                    group = baseline[(baseline.year == year) & (baseline.exposure == exposure)
                                     & (baseline.model == model)]
                values = [ci(group[group.outcome == outcome].iloc[0], 2, 'difference_pp')
                          for outcome in ['dil', 'either', 'Change in disparity']]
                rows.append([year, label, model, *values])
    tables.append(table('Supplementary Table S16 Alternative adjustment models',
                        ['Year', 'Exposure', 'Model', 'Dilation gap', 'Either gap', 'Gap change'],
                        rows, [.48, .72, .98, 1.44, 1.44, 1.44],
                        'Marginal differences in percentage points (95% CI). M0 includes state and the focal exposure. M3 excludes insulin from M1. M4 adds approximate diabetes duration and reported type to M1, with unknown categories. M0, M3 and M4 retain the M1 sample. Binary models reproduce the previous separate fits and use a paired influence-function difference for gap change.'))

    overlap = pd.read_csv(results / 'cost_overlap.csv')
    rows = []
    common_labels = {'state_name': 'State', 'age_group': 'Age', 'sex': 'Sex', 'race': 'Race/ethnicity',
                     'education': 'Education', 'nonmetro': 'Metropolitan status'}
    for (year, characteristic), group in overlap.groupby(['year', 'characteristic'], sort=False):
        positive = group[group.cost == 1]
        negative = group[group.cost == 0]
        rows.append([year, common_labels[characteristic], int(positive.n.min()), int(negative.n.min()),
                     f'{positive.weighted_pct.min():.1f}–{positive.weighted_pct.max():.1f}',
                     f'{min(positive.propensity_p01):.3f}–{max(positive.propensity_p99):.3f}',
                     f'{min(negative.propensity_p01):.3f}–{max(negative.propensity_p99):.3f}'])
    tables.append(table('Supplementary Table S17 Cost barrier overlap',
                        ['Year', 'Grouping', 'Minimum n with barrier', 'Minimum n without barrier',
                         'Weighted with barrier % range', 'Propensity range with barrier', 'Propensity range without barrier'],
                        rows, [.47, 1.05, .89, .92, 1.02, 1.075, 1.075],
                        'Ranges are across levels within each grouping (minimum first percentile to maximum 99th percentile for fitted propensities). Logistic propensity diagnostics condition on state, nonmetropolitan status, age, sex, race and ethnicity, education, and insulin. See cost_overlap.csv for every state and subgroup. Overlap in these summaries does not guarantee overlap in every covariate combination.'))

    sensitivities = pd.read_csv(results / 'sensitivities.csv')
    two_year = sensitivities[sensitivities.group.eq('Two-year window')]
    rows = []
    for measure in ['Dilation', 'Photography', 'Either', 'Both', 'Dilation only', 'Photography only', 'Neither']:
        r25, r22 = [two_year[(two_year.year == y) & (two_year.measure == measure)].iloc[0] for y in (2025, 2022)]
        rows.append([measure, f'{r25.events:,}/{r25.n:,}', ci(r25, 2), f'{r22.events:,}/{r22.n:,}', ci(r22, 2)])
    tables.append(table('Supplementary Table S18 Two year eye care classification',
                        ['Measure', '2025 n/N', '2025 % (95% CI)', '2022 n/N', '2022 % (95% CI)'],
                        rows, [1.4, 1.1, 1.45, 1.1, 1.45],
                        'All measures use the same complete-pair denominator as the annual analysis, with the additional one-to-two-year category counted as contact.'))

    missing = pd.read_csv(results / 'module_missingness.csv')
    support = pd.read_csv(results / 'module_response_support.csv')
    rows = [[r.year, r.category, f'{r.n:,}', f'{r.weighted_pct:.2f}', f'{r.completed_n:,}', f'{r.partial_n:,}']
            for r in missing.itertuples()]
    insulin = support[support.field.eq('INSULIN1') & support.level.ne('Valid')]
    insulin_text = ' '.join(
        f'In {year}, {g.loc["Blank", "n"]:,} insulin blanks yielded no complete pairs. '
        f'Of {g.loc["Explicit unknown or refusal", "n"]:,} explicit unknown or refused insulin responses, '
        f'{g.loc["Explicit unknown or refusal", "complete_n"]:,} had complete eye pairs.'
        for year, group in insulin.groupby('year') for g in [group.set_index('level')])
    tables.append(table('Supplementary Table S19 Diabetes module response and interview disposition',
                        ['Year', 'Recorded pattern', 'n', 'Weighted %', 'Completed n', 'Partial n'],
                        rows, [.52, 2.61, .61, .88, .94, .94],
                        'Each respondent with diagnosed diabetes in eligible states is assigned one pattern. '
                        'All module fields blank means that all seven public-use module fields were blank. '
                        'An explicit eye response of do not know or refused takes precedence over a blank on the other item. '
                        'CDC disposition codes 1100 and 1200 describe the interview overall. '
                        'A record with no module answers does not reveal where the interview stopped or why the items are blank. '
                        + insulin_text + ' Counts by state, questionnaire version and disposition and support by earlier module item accompany the table.'))

    impact = pd.read_csv(results / 'ipw_ranked_influence.csv')
    rows = []
    for (year, specification, size), block in impact.groupby(['year', 'specification', 'top_n'], sort=False):
        uncapped, capped = [block.loc[block.cap.eq(cap)].iloc[0] for cap in ('Uncapped', '2.0')]
        rows.append([year, specification, size, uncapped.top_events,
                     f'{uncapped.share_total_weight_pct:.2f}', f'{uncapped.share_event_weight_pct:.2f}',
                     f'{uncapped.direct_contribution_pp:.2f} / {capped.direct_contribution_pp:.2f}',
                     f'{uncapped.without_top_pct:.2f} / {capped.without_top_pct:.2f}'])
    tables.append(table('Supplementary Table S20 Influence of the largest response multipliers',
                        ['Year', 'Response model', 'Top n', 'Photo only n', 'Weight share %',
                         'Event weight share %', 'Contribution pp uncapped / cap 2', 'Estimate % without top n uncapped / cap 2'],
                        rows, [.42, .94, .43, .56, .67, .86, 1.22, 1.4],
                        'Rows rank complete pairs by their uncapped inverse response multiplier within model and year. '
                        'Weight share uses the denominator of the uncapped ratio. Event-weight share uses its photography-only numerator. '
                        'Contribution divides the selected observations’ weighted photography-only numerator by the full adjusted denominator of the corresponding model. '
                        'The estimate without top n removes their numerator and denominator terms with propensities held fixed. '
                        'The capped columns use the same ranked observations. These are influence diagnostics. '
                        'Full results also include the 1.5 cap.'))
    return tables

