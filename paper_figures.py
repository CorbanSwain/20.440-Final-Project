#!/bin/python3
# paper_figures.py
# Corban Swain 2018

from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2, venn2_circles
from pyvenn import venn
from mpl_toolkits import axes_grid1
import seaborn as sns
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance
from sklearn.metrics import silhouette_score
import textwrap


saturation_val = 1
pal = sns.color_palette('tab20b', 20, saturation_val)
pal = pal[0:4] + pal[5:8] + [pal[10], pal[13], pal[18]]
default_palette = pal
ptsz = 2 ** 2
lw = 1
savefig_args = dict(transparent=True, dpi=300)


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


clust_args = dict(vmin=-1, vmax=1, figsize=(4, 4), center=0, cmap='vlag',
                  cbar_kws=dict(ticks=[-1, 0, 1], format='%d', drawedges=False),
                  dgline_kws=dict(linewidth=1))


def make_corr_cluster(cds, filt=None, title='untitled'):
    df = cds.profile_panda(filt)
    cancers = cds.dss_names

    # categorical palette
    lut = dict(zip(cancers, default_palette))
    column_names = df.columns.get_level_values('cancer_type')
    cancer_colors = pd.Series(column_names, index=df.columns).map(lut)

    c = sns.clustermap(df.corr(),
                       row_colors=cancer_colors, col_colors=cancer_colors,
                       xticklabels=False, row_cluster=False,
                       col_cluster=False, yticklabels=False, **clust_args)
    plt.savefig(os.path.join('figures', title), **savefig_args)


def make_corr_cluster_2(cds, filt=None, title='untitled'):
    df = cds.median_profile_panda(filt)
    cancers = cds.dss_names

    corr_df = df.corr()
    dists = distance.pdist(np.array(corr_df))
    linkage = hierarchy.linkage(dists, method='average')

    # categorical palette
    lut = dict(zip(cancers, default_palette))

    column_names = df.columns
    cancer_colors = pd.Series(column_names, index=df.columns).map(lut)
    lw = None
    c = sns.clustermap(corr_df, row_linkage=linkage, col_linkage=linkage,
                       row_colors=cancer_colors, col_colors=cancer_colors,
                       linewidth=lw, **clust_args)
    plt.savefig(os.path.join('figures', title), **savefig_args)

    plt.figure(figsize=(3, 3.2))
    dist_mat = distance.squareform(dists + dists.T)
    x = np.arange(1, len(df.columns)) + 1
    y = np.zeros(x.shape)
    for i, maxc in enumerate(x):
        nodes = hierarchy.fcluster(linkage, t=maxc, criterion='maxclust')
        y[i] = silhouette_score(dist_mat + dist_mat.T, nodes,
                                metric='euclidean')
    ax = plt.subplot(111)
    ax.plot(x, y, linewidth=1.5, marker='o', color=default_palette[0],
            ms=7)
    ax.plot(x[4], y[4], linewidth=1.5, marker='o',
            markerfacecolor=default_palette[7],
            markeredgecolor=default_palette[0],
            markeredgewidth=1,
            ms=7)
    sns.despine(ax=ax)
    ax.set_ylabel('Silhouette Score')
    ax.set_xlabel('Number of Clusters')
    ax.tick_params(direction='in', length=2.5, pad=2.5)
    plt.tight_layout(pad=2.5)
    plt.savefig(os.path.join('figures', title + '_sillohoute'), **savefig_args)


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

    plt.figure(1, (5, 2.8))
    ax = plt.subplot(111)
    bar_args = dict(x='cancer', linewidth=lw, edgecolor='k', data=df)
    sns.barplot(y='expressed', facecolor=(1, 1, 1, 0),
                label='Total Expressed', **bar_args)
    sns.barplot(y='significant', label='Differentially Expressed',
                palette=default_palette, **bar_args)
    ax.legend(ncol=2, loc='upper center', framealpha=0,
              bbox_to_anchor=(0.5, -0.15))
    plt.tight_layout()
    ax.set(ylabel='Number of lncRNAs', xlabel='')
    sns.despine(ax=ax, right=True, top=True)
    plt.savefig(os.path.join('figures', title), **savefig_args)


def make_type_volcanos(cds, title=None):
    fig = plt.figure(figsize=(6.2, 3.1))
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
            ax.set_xlabel('')
            ax.set_ylabel('')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
        ax.set_xlim(-maxx, maxx)
        ax.set_ylim(0, maxy + 2.5)
        ax.text(-maxx + 1.5, maxy + 2.2, cds.dss_names[i], va='top', ha='left',
                fontweight='900', fontsize=7.5)
        ax.tick_params(direction='in', length=2.5, pad=2.5)
        sns.despine(ax=ax)

    ax = fig.add_subplot(111)
    ax.set_xlabel('log$_2$ fold change', labelpad=12)
    ax.set_ylabel('log$_{10}$ q-value', labelpad=12)
    sns.despine(ax=ax, left=True, bottom=True)
    ax.tick_params(labelcolor=(), top=False, bottom=False, left=False,
                   right=False, labeltop=False, labelbottom=False,
                   labelleft=False, labelright=False)
    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title), **savefig_args)


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
    figsz = (4.5, 2.8)
    fig = plt.figure(figsize=figsz)

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
    ax.tick_params(direction='in', length=2.5, pad=2.5)
    ax.set_xlim(-0.5, 9.8)
    sns.despine(ax=ax, right=True, top=True)
    ax.set_ylabel('log$_{2}$ fold change')
    ax.set_xlabel('Cancer Type')

    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title), **savefig_args)

    # PLot 2
    fig = plt.figure(figsize=figsz)

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
                   linewidth=1, inner=None, cut=0, width=1,
                   palette=default_palette)
    ax.tick_params(direction='in', length=2.5, pad=2.5)
    sns.despine(ax=ax, right=True, top=True)
    ax.set_ylabel('log$_{2}$ fold change')
    ax.set_xlabel('Cancer Type')

    plt.tight_layout(pad=1.15)
    plt.savefig(os.path.join('figures', title + '_not_signif'), **savefig_args)


def make_cancer_table(cds):
    for ds in cds.dss:
        num_signif = np.count_nonzero(ds.results['t_test'][2])
        num_expr = np.count_nonzero(ds.results['is_expressed'])
        fmt1 = '%4s, %37s, %3d, %3d, %3d, %4d, %4d, %2d%%'
        fmt2 = '%s,%s,%d,%d,%d,%d,%d,%d%%'
        print(fmt2
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

    plt.figure(figsize=(6, 2.2))
    ax = plt.subplot(121)
    wed, txts = plt.pie([valuesnil, ] + values[:2] + [sum(values[2:]), ],
                     labels=groups1,
                     colors=sns.color_palette('tab10', 10, 1)[:4],
                     explode=[0, 0, 0, 0.1])
    t = txts[1]
    x, y = t.get_position()
    t.set_position((x, y - 0.05))
    for w in wed:
        w.set_linewidth(1)
        w.set_edgecolor('k')
    ax.add_artist(plt.Circle((0, 0), 0.6, facecolor='w', edgecolor='k',
                             linewidth=1))
    ax.text(0, 0, '{:,}'.format(cds.n_genes), fontsize=12, fontweight='700',
            va='center', ha='center')
    ax.axis('equal')

    ax = plt.subplot(122)
    wed, txts = plt.pie(values[2:], labels=groups2,
                        colors=sns.light_palette(
                            sns.color_palette('tab10', 10, 1)[3], len(groups2)))
    t = txts[5]
    x, y = t.get_position()
    t.set_position((x, y + 0.1))
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
    plt.savefig(os.path.join('figures', title), **savefig_args)


def report_expression(cds, names):
    df = cds.long_panda
    for n in names:
        row_sel = df['gene_name'] == n
        info = df[['gene_name', 'expression', 'cancer_type', 'fc',
                   'is_expressed']][row_sel]
        print('\n\n' + n)
        print(info)



test_ls = ['PART1','MIR22HG','TCL6','SMIM10L2B','SNHG12','DUXAP8','EWSAT1',
           'MEG3','MIR99AHG','CYTOR','TINCR','LINP1','DANCR','BMS1P17',
           'MIR34AHG','DNM3OS','MIR600HG','PVT1','MIR200CHG','PWAR6','PTCSC3',
           'LUCAT1','DGCR9']

def make_table_heatmaps(cds, title=None):
    df = cds.long_panda
    grps = [['KIRC', 'KIRP'], ['BRCA', 'HNSC', 'LUAD', 'STAD']]
    for gp in grps:
        overlaps = np.ones(cds.n_genes, dtype=bool)
        for cnc in gp:
            cancer_sel = df['cancer_type'] == cnc
            overlaps = np.logical_and(overlaps, np.array(df['is_signif'][
                cancer_sel]))

        if gp is grps[0]:
            overlaps = np.zeros(cds.n_genes, dtype=bool)
            for i, gn in enumerate(cds.gene_names):
                if gn in test_ls:
                    overlaps[i] = True
        pd_dict = {}
        for cnc in gp:
            cancer_sel = df['cancer_type'] == cnc
            fc = np.array(df['fc'][cancer_sel])
            pd_dict[cnc] = fc[overlaps]

        idxs = [cds.gene_names[i] for i in range(cds.n_genes) if overlaps[i]]
        df2 = pd.DataFrame(pd_dict, index=idxs)

        if gp is grps[0]:
            fig = plt.figure(figsize=(3, 4.8))
        else:
            fig = plt.figure(figsize=(3, 4))

        sns.heatmap(df2, center=0, cmap='vlag', linewidths=1, annot=True,
                    fmt='.1f', cbar=False)
        plt.tight_layout(pad=2.5)
        plt.savefig(os.path.join('figures', title + '_'.join(gp)),
                    **savefig_args)


    max_genes_sel = df['num_signif'] == max(df['num_signif'])
    max_genes = np.unique(np.array(df['gene_name'][max_genes_sel]))
    vals = []
    for gn in max_genes:
        gn_selec = df['gene_name'] == gn
        selec = np.logical_and(max_genes_sel, gn_selec)
        vals.append(np.mean(df['fc'][selec]))

    df3 = pd.DataFrame(vals, index=max_genes, columns=['Mean Fold Change', ])
    fig = plt.figure(figsize=(2.8, 3))
    sns.heatmap(df3, center=0, cmap='vlag', linewidths=1, annot=True,
                fmt='.1f', cbar=False)
    plt.tight_layout(pad=2.5)
    plt.savefig(os.path.join('figures', title + '_max_expr'),
                **savefig_args)


def make_venn_diagrams(cds, title=None):
    grps = [['KIRC', 'KIRP'], ['BRCA', 'HNSC', 'LUAD', 'STAD']]
    df = cds.long_panda
    gp = grps[0]
    sets = []
    for cnc in gp:
        cancer_sel = df['cancer_type'] == cnc
        temp = np.array(df['is_signif'][cancer_sel])
        sets.append(set(np.where(temp)[0]))

    plt.figure(figsize=(3, 2))
    v = venn2(sets, tuple(gp))
    v.get_patch_by_id('10').set_color(default_palette[7])
    v.get_patch_by_id('01').set_color(default_palette[6])
    v.get_patch_by_id('11').set_color(default_palette[2])
    c = venn2_circles(sets, linestyle='solid')
    [cl.set_lw(1) for cl in c]
    plt.savefig(os.path.join('figures', title + '_gp_1'), **savefig_args)

    gp = grps[1]
    sets = []
    for cnc in gp:
        cancer_sel = df['cancer_type'] == cnc
        temp = np.array(df['is_signif'][cancer_sel])
        sets.append(set(np.where(temp)[0]))

    labels = venn.get_labels(sets)
    fig, ax = venn.venn4(labels, figsize=(6, 6), names=gp, colors=
                         [(default_palette[i][0],
                           default_palette[i][1],
                           default_palette[i][2],
                           0.3)
                          for i in [3, 5, 7, 8]])
    plt.savefig(os.path.join('figures', title + '_gp_2'), **savefig_args)