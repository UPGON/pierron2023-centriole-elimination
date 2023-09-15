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
    display_posterior = False

    model = 'poisson_glm'
    with open(f'../models/{model}.stan', 'r') as f:
        model_str = f.read()

    comps = [
        ('WT', 'Tubulin', ('EP', 'LP')),
        ('WT', 'Tubulin', ('LP', 'Diplo')),
        ('WT', 'SAS-4', ('EP', 'LP')),
        ('WT', 'SAS-4', ('LP', 'Diplo')),
        ('sas-1ts', 'Tubulin', ('EP', 'LP')),
        ('sas-1ts', 'Tubulin', ('LP', 'Diplo')),
        ('sas-1ts', 'SAS-4', ('EP', 'LP')),
        ('sas-1ts', 'SAS-4', ('LP', 'Diplo')),
    ]

    data = pd.read_csv('../data/data.csv')
    genotype_type = pd.CategoricalDtype(categories=['WT', 'sas-1ts'], ordered=True)
    data['genotype'] = data['genotype'].astype(genotype_type)

    for genotype, protein, (phase1, phase2) in comps:
        print(genotype, protein, phase1, phase2)
        sub = data.loc[(data['protein'] == protein) &
                       (data['genotype'] == genotype) &
                       (data['phase'].isin((phase1, phase2)))]

        data_dict = {
            "N": len(sub),
            "J": len(pd.unique(sub['phase'])),
            "x": pd.factorize(sub['phase'])[0],
            "y": np.array(sub['intensity'].astype('int')),
        }

        posterior = stan.build(model_str, data=data_dict, random_seed=1993)
        fit = posterior.sample(num_chains=4, num_samples=8000)

        postpred = pd.DataFrame(np.repeat([0, 1], 1000), columns=['group'])
        postpred['b0'] = fit.to_frame()['b0']
        postpred['b1'] = fit.to_frame()['b1']
        b0_mean = postpred['b0'].mean().round(2)
        b0_std = postpred['b0'].std().round(2)
        b1_mean = postpred['b1'].mean().round(2)
        b1_std = postpred['b1'].std().round(2)
        postpred = postpred.sample(100)

        fig, ax = plt.subplots(1, 1, figsize=(1.7, 1.7))
        sns.stripplot(data=sub, x='phase', y='intensity', hue='genotype', ax=ax, legend=False, color='black')

        if display_posterior:
            for i, rec in postpred.iterrows():
                ax.plot([0, 1], [np.exp(rec['b0']),
                                 np.exp(rec['b0'] + rec['b1'])],
                        color='black', alpha=.2, lw=.5)

        ax.plot([0, 1], [np.exp(postpred['b0'].mean()),
                         np.exp(postpred['b0'].mean() + postpred['b1'].mean())],
                color='black', alpha=1, lw=1)

        ax.set_xlabel('Phase')
        ax.set_ylabel('Intensity')
        ax.set_ylim(0)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        alpha = r'$\alpha$'
        beta = r'$\beta$'
        ax.set_title(
            f"{protein}, {genotype}, {phase1} vs {phase2}\n{alpha}={b0_mean:.2f}±{b0_std}\n{beta}={b1_mean}±{b1_std}")

        fig.tight_layout()
        fig.savefig(f'../plots/png/figure_{genotype}_{protein}_{phase1}_{phase2}.png', dpi=300)
        fig.savefig(f'../plots/eps/figure_{genotype}_{protein}_{phase1}_{phase2}.eps', dpi=300, format='eps')


if __name__ == '__main__':
    main()
