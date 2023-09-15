import stan
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
matplotlib.rcParams.update({'font.size': 7,
                            'font.family': 'Helvetica'})


def main():
    model = 'poisson_glm'
    with open(f'../models/{model}.stan', 'r') as f:
        model_str = f.read()

    comps = [
        ('WT', 'Tubulin'),
        ('WT', 'SAS-4'),
        ('sas-1ts', 'Tubulin'),
        ('sas-1ts', 'SAS-4'),
    ]
    phase_comps = [
        ((0, 1), 'EP', 'LP'),
        ((1, 2), 'LP', 'Diplo'),
    ]

    data = pd.read_csv('../data/data.csv')
    genotype_type = pd.CategoricalDtype(categories=['WT', 'sas-1ts'], ordered=True)
    data['genotype'] = data['genotype'].astype(genotype_type)
    summary = []
    fig, axs = plt.subplots(1, 4, figsize=(4 * 1.7, 1.7))
    for i, (genotype, protein) in enumerate(comps):
        ax = axs[i]

        ax.set_xlabel('Phase')
        ax.set_ylabel('Intensity')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        sub = data.loc[(data['protein'] == protein) &
                       (data['genotype'] == genotype)]

        sns.stripplot(data=sub, x='phase', y='intensity', ax=ax, legend=False, color='black')

        for comp in phase_comps:
            x_pos, phase1, phase2 = comp
            ssub = sub.loc[(data['phase'].isin(comp))]

            data_dict = {
                "N": len(ssub),
                "J": len(pd.unique(ssub['phase'])),
                "x": pd.factorize(ssub['phase'])[0],
                "y": np.array(ssub['intensity'].astype('int')),
            }

            posterior = stan.build(model_str, data=data_dict, random_seed=1993)
            fit = posterior.sample(num_chains=4, num_samples=8000)

            post_preds = pd.DataFrame(np.repeat([0, 1], 1000), columns=['group'])
            post_preds['b0'] = fit.to_frame()['b0']
            post_preds['b1'] = fit.to_frame()['b1']
            b0_mean = post_preds['b0'].mean().round(2)
            b0_std = post_preds['b0'].std().round(2)
            b1_mean = post_preds['b1'].mean().round(2)
            b1_std = post_preds['b1'].std().round(2)
            post_preds = post_preds.sample(100)

            summary.append(f"{phase1} vs {phase2}, a={b0_mean:.2f}±{b0_std}, b={b1_mean}±{b1_std}\n")

            ax.plot(x_pos, [np.exp(post_preds['b0'].mean()),
                            np.exp(post_preds['b0'].mean() + post_preds['b1'].mean())],
                    color='black', alpha=1, lw=1, ls='--')

            ax.set_title(
                f"{genotype}, {protein}")

    fig.tight_layout()
    fig.savefig(f'../plots/png/figure7ef.png', dpi=300)
    fig.savefig(f'../plots/eps/figure7ef.eps', dpi=300, format='eps')
    with open('../samples/summary.txt', 'w') as f:
        f.writelines(summary)


if __name__ == '__main__':
    main()
