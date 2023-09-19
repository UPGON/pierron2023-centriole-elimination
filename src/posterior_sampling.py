import stan
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import arviz as az

import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
matplotlib.rcParams.update({'font.size': 7,
                            'font.family': 'Helvetica'})

model = 'poisson_glm'
intensity = 'intensity_scaled_ep'
with open(f'../models/{model}.stan', 'r') as f:
    model_str = f.read()

comps = [
    ('WT', 'Tubulin'),
    ('sas-1ts', 'Tubulin'),
    ('WT', 'SAS-4'),
    ('sas-1ts', 'SAS-4'),
]
phase_comps = [
    ((0, 1), 'EP', 'LP'),
    ((1, 2), 'LP', 'Diplo'),
]


def format_ax(ax):
    ax.set_xlabel('Phase')
    ax.set_ylabel('Intensity')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    return ax


def main():
    data = pd.read_csv('../data/data.csv')

    fig, axs = plt.subplots(1, 4, figsize=(4 * 1.7, 1.7))
    summary = []

    for i, (genotype, protein) in enumerate(comps):
        sub = data.loc[(data['protein'] == protein) &
                       (data['genotype'] == genotype)].copy()

        mean_ep = sub.loc[sub['phase'] == 'EP', 'intensity'].mean()
        # Multiply by 100 so that the Poisson counts are not 0 or 1 because of the floating point
        sub['intensity_scaled_ep'] = 100 * sub['intensity'] / mean_ep

        ax = axs[i]
        ax = format_ax(ax)

        if intensity == 'intensity_scaled_ep':
            ax.set_ylim(0, 200)

        sns.stripplot(data=sub, x='phase', y=intensity, ax=ax, legend=False, color='black', size=2, jitter=True)

        for comp in phase_comps:
            x_pos, phase1, phase2 = comp
            ssub = sub.loc[(data['phase'].isin(comp))]
            identifier = f'{genotype}_{protein}_{phase1}_{phase2}_{intensity}'
            data_dict = {
                "N": len(ssub),
                "J": len(pd.unique(ssub['phase'])),
                "x": pd.factorize(ssub['phase'])[0],
                "y": np.array(ssub[intensity].astype(int)),
            }

            posterior = stan.build(model_str, data=data_dict, random_seed=1993)
            fit = posterior.sample(num_chains=4, num_samples=8000)

            fit.to_frame().round(3).to_csv(f'../samples/posterior_{identifier}.tsv')
            az.hdi(fit, hdi_prob=.95).to_dataframe().round(3).to_csv(f'../samples/hdi_{identifier}_95.tsv')

            post_preds = pd.DataFrame(np.repeat([0, 1], 1000), columns=['group'])
            post_preds['b0'] = fit.to_frame()['b0']
            post_preds['b1'] = fit.to_frame()['b1']
            b0_mean = post_preds['b0'].mean().round(2)
            b0_std = post_preds['b0'].std().round(2)
            b1_mean = post_preds['b1'].mean().round(2)
            b1_std = post_preds['b1'].std().round(2)

            stats = {'genotype': genotype,
                     'protein': protein,
                     'phase_ref': phase1,
                     'phase_next': phase2,
                     'b0_mean': b0_mean,
                     'b0_std': b0_std,
                     'b1_mean': b1_mean,
                     'b1_std': b1_std,
                     }
            summary.append(stats)

            ax.plot(x_pos, [np.exp(post_preds['b0'].mean()),
                            np.exp(post_preds['b0'].mean() + post_preds['b1'].mean())],
                    color='black', alpha=1, lw=1, ls='--')

            ax.set_title(
                f"{genotype}, {protein}")

    fig.tight_layout()
    fig.savefig(f'../plots/png/figure7ef_{intensity}.png', dpi=300)
    fig.savefig(f'../plots/eps/figure7ef_{intensity}.eps', dpi=300, format='eps')
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f'../samples/summary_{intensity}.csv')


if __name__ == '__main__':
    main()
