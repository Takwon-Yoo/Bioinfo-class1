import sys
import directory as DIR
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


ORF_palette = {
    'ORF1a': '#fbde8e',
    'ORF1b': '#cce8ff',  # '#99ccff',
    'S': '#d62728',
    'ORF3a': '#8c564b',
    'E': '#1b1464',  # '#17becf',
    'M': '#1f77b4',
    'ORF6': '#9467bd',
    'ORF7a': '#ff7f0e',
    'ORF7b': '#067f4c',
    'ORF8': '#e377c2',
    'N': '#bcbd22',
    'ORF10': '#3bcf49',  # '#6f9379',
}

font = {'size': 16}
matplotlib.rc('font', **font)


def main_read_counts(align):
    if align == 'star':
        plotdir = f'{DIR.STAR_DIR}read_counts.png'
        datadir = f'{DIR.STAR_DIR}'
    else:
        plotdir = f'{DIR.BWA_DIR}read_counts.png'
        datadir = f'{DIR.BWA_DIR}'
    xticks = DIR.HPI
    ways = ['TR', 'RF']

    fig, axs = plt.subplots(nrows=1, ncols=2, sharex='all', figsize=(8, 6))

    for idx, way in enumerate(ways):
        ax = axs.flat[idx]
        plot_read_counts(ax, way, datadir)

    plt.xticks(np.arange(0, 8), labels=xticks)
    plt.tight_layout()
    plt.savefig(plotdir)
    plt.close()
    return None


def plot_read_counts(ax, way, datadir):
    df = pd.read_csv(f'{datadir}{way}_avg_read_counts.txt', sep='\t', index_col=0)
    df.drop(labels=[col for col in df.columns if 'hpi' not in col], axis=1, inplace=True)
    print(df)
    df = np.log10(df + 1)

    for gene in df.index.tolist():
        try:
            ax.plot(df.loc[gene, :], linewidth=4, color=ORF_palette[gene], label=gene)
        except KeyError:
            pass
    ax.set_xlabel('hpi')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if way == 'TR':
        ax.text(0.5, 5, f'mRNA-seq')
        ax.set_ylabel(r'''$\mathregular{log_{10}}$(ReadCounts+1)''')
    else:
        ax.text(0.5, 3.8, f'RPF-seq')
    plt.legend(bbox_to_anchor=(1.05, 1.2))

    return None


def main_TE(align):
    xticks = DIR.HPI

    if align == 'star':
        plotdir = f'{DIR.STAR_DIR}TE.png'
        datadir = f'{DIR.STAR_DIR}'
    else:
        plotdir = f'{DIR.BWA_DIR}TE.png'
        datadir = f'{DIR.BWA_DIR}'

    df_mRNA = pd.read_csv(f'{datadir}TR_avg_read_counts.txt', sep='\t', index_col=0)
    orf1a = df_mRNA.loc['ORF1a', :]
    orf1b = orf1a.rename('ORF1b')
    df_mRNA = df_mRNA.append(orf1b)
    df_mRNA.loc['ORF7b', :] = df_mRNA.loc['ORF7a', :]
    df_mRNA.drop(labels=[col for col in df_mRNA.columns if 'hpi' not in col], axis=1, inplace=True)

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex='all', figsize=(6, 6))

    plot_TE(ax, datadir, df_mRNA)

    title = 'Translation efficiency (TE) = ' + r'''${\frac{RPF\ Readcount+1}{mRNA\ Readcount+1}}$'''

    plt.suptitle(title, fontsize=17)
    plt.xticks(np.arange(0, 8), labels=xticks)
    plt.tight_layout()
    plt.savefig(plotdir)
    plt.close()
    return None


def plot_TE(ax, datadir, df_mRNA):
    df = pd.read_csv(f'{datadir}RF_avg_read_counts.txt', sep='\t', index_col=0)
    df.drop(labels=[col for col in df.columns if 'hpi' not in col], axis=1, inplace=True)
    print(df, df_mRNA)
    for gene in df.index.tolist():
        try:
            ax.plot(np.log10((df.loc[gene, :] + 1) / (df_mRNA.loc[gene, :] + 1)), linewidth=3, color=ORF_palette[gene])
        except KeyError:
            pass
    ax.set_xlabel('hpi')
    # ax.set_yticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_ylabel(r'''$\mathregular{log_{10}}$(TE)''')

    ax.text(0.5, 1.2, f'TE(RPF/mRNA)')

    return None


MODE = sys.argv[1]
PARAM = sys.argv[2:]
if __name__ == '__main__':
    if MODE in locals().keys():
        locals()[MODE](*PARAM)
    else:
        sys.exit('error: cmd=%s' % MODE)