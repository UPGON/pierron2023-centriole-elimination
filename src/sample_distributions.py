import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def main():
    data = pd.read_csv('../data/data.csv')
    genotypes = pd.unique(data['genotype'])
    proteins = pd.unique(data['protein'])
    phases = pd.unique(data['phase'])
    colors = {'WT': 'k', 'sas-1ts': 'b'}
    print(proteins)
    print(genotypes)
    fig, axs = plt.subplots(nrows=len(phases), ncols=len(proteins), figsize=(4, 6), dpi=300)
    for h, pha in enumerate(phases):
        for p, prot in enumerate(proteins):
            ax = axs[h, p]
            for g, gen in enumerate(genotypes):
                sub = data.loc[(data['protein'] == prot) & (data['genotype'] == gen) & (data['phase'] == pha)]
                print(sub)
                ax.hist(sub['intensity'], color=colors[gen], bins=np.arange(0, 35+2, step=2),
                        density=True, histtype='bar', alpha=.5, label=gen)
                ax.set_xlim(0, 30)
                ax.set_ylim(0, 1)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_title(f"{prot}, {gen}, {pha}")
                if p == 0 and h == 0:
                    ax.legend()
                    ax.set_xlabel('Intensity')
                    ax.set_ylabel('Counts')
    fig.tight_layout()
    fig.savefig('../plots/plot_distributions.png')
    print(data)


if __name__ == '__main__':
    main()
