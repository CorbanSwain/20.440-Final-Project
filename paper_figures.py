#!/bin/python3
# paper_figures.py
# Corban Swain 2018

from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
from pyvenn import venn
from mpl_toolkits import axes_grid1
import seaborn as sns
import pandas as pd
import textwrap


saturation_val = 1
pal = sns.color_palette('tab20b', 20, saturation_val)
pal = pal[0:4] + pal[5:8] + [pal[10], pal[13], pal[18]]
default_palette = pal
ptsz = 2 ** 2


def plot_hv_line(direction, position, lims=[-100, 100], axis=None, **kwargs):
    z = np.ones(2) * position;
    a = np.array(lims)
    if axis is None:
        axis = plt.gca()
    if direction is 'h':
        axis.plot(a, z, **kwargs)
    elif direction is 'v':
        axis.plot(z, a, **kwargs)
    else:
        raise ValueError('direction must be \'h\' or \'v\'')


def make_corr_cluster(cds, filt=None, title='untitled'):
    df = cds.profile_panda(filt)
    cancers = cds.dss_names

    # categorical palette
    lut = dict(zip(cancers, default_palette))

    column_names = df.columns.get_level_values('cancer_type')
    cancer_colors = pd.Series(column_names, index=df.columns).map(lut)

    c = sns.clustermap(df.corr(), center=0, cmap='vlag',
                       row_colors=cancer_colors, col_colors=cancer_colors,
                       figsize=(8, 8), xticklabels=False, row_cluster=False,
                       col_cluster=False, yticklabels=False)
    plt.savefig(os.path.join('figures', title), dpi=300)


def make_corr_cluster_2(cds, filt=None, title='untitled'):
    df = cds.median_profile_panda(filt)
    cancers = cds.dss_names

    # categorical palette
    lut = dict(zip(cancers, default_palette))

    column_names = df.columns
    cancer_colors = pd.Series(column_names, index=df.columns).map(lut)

    c = sns.clustermap(df.corr(), center=0, cmap='vlag',
                       row_colors=cancer_colors, col_colors=cancer_colors,
                       figsize=(5, 5), linewidth=0.75)
    plt.savefig(os.path.join('figures', title), dpi=300)


def make_barchart(datasets, title=None):
    # Significance count histogram
    n_datasets = len(datasets)
    panda_dict = {
        'expressed': np.zeros(n_datasets, dtype=int),
        'significant': np.zeros(n_datasets, dtype=int),
        'cancer': np.zeros(n_datasets, dtype=np.object)
    }
    for i, ds in enumerate(datasets):
        num_expressed = np.count_nonzero(ds.results['is_expressed'])
        _, _, is_signif = ds.results['t_test']
        num_signif = np.count_nonzero(is_signif)

        panda_dict['significant'][i] = num_signif
        panda_dict['expressed'][i] = num_expressed
        panda_dict['cancer'][i] = ds.cancer_type

    df = pd.DataFrame(panda_dict)
    # df = df.sort_values('significant', ascending=False)

    plt.figure(1, (9, 4))
    ax = plt.subplot(111)
    sns.barplot(y='expressed', x='cancer',
                data=df, label='Total Expressed', color=(0.8, 0.8, 0.8))
    sns.barplot(y='significant', x='cancer',
                data=df, label='Differentially Expressed',
                palette=default_palette)
    ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.1))
    plt.tight_layout()
    ax.set(ylabel='Number of lncRNAs', xlabel='')
    sns.despine(ax=ax, right=True, top=True, left=True)
    plt.savefig(os.path.join('figures', title), dpi=300)


def make_type_volcanos(cds, title=None):
    fig = plt.figure(figsize=(8, 5))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(2, 5), axes_pad=0.1,
                           label_mode='L', share_all=True)
    df = cds.long_panda
    x = df['fc']
    y = -np.log10(df['q_values'])
    s = df['t_stat'] ** 2
    s = s / max(s) * 15
    maxy = max(y[df['is_expressed']])
    maxx = max(np.abs(x[df['is_expressed']]))
    for i in range(cds.n_dss):
        base_sel = np.logical_and(df['cancer_num_id'] == i,
                                  df['is_expressed'])
        not_sig = np.logical_and(base_sel, np.logical_not(df['is_signif']))
        sig = np.logical_and(base_sel, df['is_signif'])
        ax = grid[i]
        sns.regplot(x[not_sig], y[not_sig], ax=ax, color=(0.3, 0.3, 0.3),
                    fit_reg=False,
                    scatter_kws=dict(s=ptsz, alpha=0.6))
        sns.regplot(x[sig], y[sig], ax=ax, color=default_palette[8],
                    fit_reg=False,
                    scatter_kws=dict(s=ptsz, alpha=0.6))
        kwa = dict(axis=ax, zorder=0, c='k', alpha=0.1, ls='--', lw=0.75)
        plot_hv_line('h', 3, **kwa)
        plot_hv_line('v', -1, [0, maxy], **kwa)
        plot_hv_line('v', 1, [0, maxy], **kwa)
        if i == 5:
            ax.set_xlabel('log2 fold change')
            ax.set_ylabel('-log10 q-value')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
        ax.set_xlim(-maxx, maxx)
        ax.set_ylim(0, maxy + 2.5)
        ax.text(-maxx + 2, maxy + 2, cds.dss_names[i], va='top', ha='left',
                fontweight='bold', fontsize=7.5,
                bbox={'boxstyle': 'round', 'facecolor': 'w', 'pad': 0.3})
        sns.despine(ax=ax, right=True, top=True)

    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title), dpi=300)


def make_num_volcanos(cds, title=None):
    fig = plt.figure(figsize=(7, 3))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(1, 7), axes_pad=0.2,
                           label_mode='L', share_all=True)
    df = cds.long_panda
    x = df['fc']
    y = -np.log10(df['q_values'])
    maxy = max(y[df['is_expressed']])
    maxx = max(np.abs(x[df['is_expressed']]))

    num_sig = range(1, 8)
    for i, ns in enumerate(num_sig):
        base_sel = df['is_expressed']
        not_sig = np.logical_and(base_sel,
                                 np.logical_not(df['num_signif'] >= ns))
        sig = np.logical_and(base_sel,
                             df['num_signif'] >= ns)
        ax = grid[i]
        sns.set_style('whitegrid')
        sns.regplot(x[not_sig], y[not_sig], ax=ax, color=(0.3, 0.3, 0.3),
                    fit_reg=False, scatter_kws={'s': ptsz})
        sns.regplot(x[sig], y[sig], ax=ax, color=default_palette[8],
                    fit_reg=False, scatter_kws={'s': ptsz})
        kwa = dict(axis=ax, zorder=0, c='k', alpha=0.1, ls='--', lw=0.75)
        plot_hv_line('h', 3, **kwa)
        plot_hv_line('v', -1, [0, maxy], **kwa)
        plot_hv_line('v', 1, [0, maxy], **kwa)
        ax.set_xlabel('log2 fold change')
        ax.set_ylabel('-log10 q-value')
        ax.set_xlim(-maxx, maxx)
        ax.set_ylim(0, maxy + 2.5)
        ax.text(-maxx + 0.5, maxy + 2, '>= %d cancers' % ns,
                va='top', ha='left', fontsize=10)
        sns.despine(ax=ax, right=True, top=True)

    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title), dpi=300)


def make_violin(cds, title=None):
    fig = plt.figure(figsize=(6, 4))

    df = cds.long_panda
    n_genes = np.count_nonzero(cds.any_expressed)
    allex_long = np.tile(cds.any_expressed, cds.n_dss)
    df_dict = {}

    temp = df['fc']
    for i, n in enumerate(cds.dss_names):
        df_dict[n] = temp[
            np.logical_and(np.logical_and(df['cancer_num_id'] == i,
                                          df['is_expressed']),
                           df['num_signif'] >= 1)
        ]

    v_df = pd.DataFrame(df_dict)
    v_df = v_df[cds.dss_names]
    ax = fig.subplots(1, 1)
    sns.violinplot(data=v_df, ax=ax, orient='v', showfliers=False, bw=0.15,
                   linewidth=1, inner=None, cut=0, width=1.25,
                   palette=default_palette)
    ax.set_xlim(-0.5, 9.8)
    sns.despine(ax=ax, right=True, top=True)
    ax.set_ylabel('log2 fold change')
    ax.set_xlabel('Cancer Type')

    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title), dpi=300)

    # PLot 2
    fig = plt.figure(figsize=(8, 5))

    for i, n in enumerate(cds.dss_names):
        df_dict[n] = temp[
            np.logical_and(np.logical_and(df['cancer_num_id'] == i,
                                          df['is_expressed']),
                           df['num_signif'] == 0)
        ]


    v_df = pd.DataFrame(df_dict)
    v_df = v_df[cds.dss_names]
    ax = fig.subplots(1, 1)
    sns.violinplot(data=v_df, ax=ax, orient='v', showfliers=False, bw=0.15,
                   linewidth=2, inner=None, cut=0, width=1,
                   palette=default_palette)

    sns.despine(ax=ax, right=True, top=True)
    ax.set_ylabel('log2 fold change')
    ax.set_xlabel('Cancer Type')

    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title + '_not_signif'), dpi=300)


def make_cancer_table(cds):
    for ds in cds.dss:
        num_signif = np.count_nonzero(ds.results['t_test'][2])
        num_expr = np.count_nonzero(ds.results['is_expressed'])
        print('%4s, %37s, %3d, %3d, %3d, %4d, %4d, %2d%%'
              % (ds.cancer_type,
                 TanricDataset.sampleid2name(ds.cancer_type),
                 ds.n_normal_samples,
                 ds.n_tumor_samples,
                 ds.n_pairs,
                 num_signif,
                 num_expr,
                 round(num_signif / num_expr * 100)))

def make_pie_chart(cds, title=None):

    values = [np.count_nonzero(cds.num_signif == i) for i in range(8)]
    values[0] = np.count_nonzero(cds.any_expressed) - sum(values[1:])
    valuesnil = np.count_nonzero(np.logical_not(cds.any_expressed))
    tot = sum(values) + valuesnil
    groups1 = ['Not Expressed (%d%%)' % round(valuesnil / tot * 100),
               'Expressed, \nNot D.E. (%d%%)' % round(values[0] / tot * 100),
               '1 Cancer (%d%%)' % round(values[1] / tot * 100),
               '2+ Cancers (%d%%)' % round(sum(values[2:]) / tot * 100)]
    groups2 = ['%d Cancers (%d%%)'
               % (i, round(values[i] / sum(values[2:]) * 100))
               for i in range(2, 8)]

    plt.figure(figsize=(8, 2.8))
    ax = plt.subplot(121)
    wed, _ = plt.pie([valuesnil, ] + values[:2] + [sum(values[2:]), ],
                     labels=groups1,
                     colors=sns.color_palette('tab10', 10, 1)[:4],
                     explode=[0, 0, 0, 0.1])
    for w in wed:
        w.set_linewidth(1)
        w.set_edgecolor('k')
    ax.add_artist(plt.Circle((0, 0), 0.6, facecolor='w', edgecolor='k',
                             linewidth=1))
    ax.text(0, 0, '{:,}'.format(cds.n_genes), fontsize=12, fontweight='bold',
            va='center', ha='center')
    ax.axis('equal')

    ax = plt.subplot(122)
    wed, _ = plt.pie(values[2:], labels=groups2,
                     colors=sns.light_palette(sns.color_palette('tab10',
                                                                10, 1)[3],
                                              len(groups2)))
    for w in wed:
        w.set_linewidth(1)
        w.set_edgecolor('k')
    ax.add_artist(plt.Circle((0, 0), 0.6, facecolor='w', edgecolor='k',
                             linewidth=1))
    ax.axis('equal')
    ax.text(0, 0, '{:,}'.format(sum(values[2:])), fontsize=12,
            fontweight='bold',
            va='center', ha='center')
    plt.tight_layout(pad=1.6)
    plt.savefig(os.path.join('figures', title), dpi=300)

