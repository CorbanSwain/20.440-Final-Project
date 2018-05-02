#!/bin/python3
# paper_figures.py
# Corban Swain 2018

from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
from mpl_toolkits import axes_grid1
import seaborn as sns
import pandas as pd
import textwrap


def make_correl_cluster(datasets):
    n_datasets = len(datasets)
    n_genes = TanricDataset.n_genes
    signif = np.zeros((n_genes, n_datasets), dtype=bool)
    fc_means = np.zeros((n_genes, n_datasets))
    n_total_samples = sum(ds.n_tumor_samples for ds in datasets)
    signif_all = np.zeros((n_genes, n_total_samples))
    fc_all = np.zeros((n_genes, n_total_samples))

    for i_ds, ds in enumerate(datasets):
        _, _, iss = ds.results['t_test']
        signif[:, i_ds] = iss



def make_simple_charts(datasets):
    # Significance count histogram
    n_datasets = len(datasets)
    panda_dict = {
        'expressed': np.zeros(n_datasets, dtype=int),
        'significant': np.zeros(n_datasets, dtype=int),
        'cancer': np.zeros(n_datasets, dtype=np.object)
    }
    row_names = {}
    for i, ds in enumerate(datasets):
        num_valid = np.count_nonzero(ds.results['is_nonzero'])
        num_expressed = np.count_nonzero(ds.results['is_expressed'])
        _, _, is_signif = ds.results['t_test']
        num_signif = np.count_nonzero(is_signif)

        panda_dict['significant'][i] = num_signif
        panda_dict['expressed'][i] = num_expressed
        panda_dict['cancer'][i] = ds.cancer_type

    df = pd.DataFrame(panda_dict)
    df = df.sort_values('significant', ascending=False)

    # sns.set(style='whitegrid')
    fig = plt.figure(1, (9, 4))
    ax = plt.subplot(111)
    sns.set_color_codes('muted')
    sns.barplot(y='expressed', x='cancer',
                data=df, label='Number Expressed', color='b')
    sns.set_color_codes('dark')
    sns.barplot(y='significant', x='cancer',
                data=df, label='Number Significant', color='b')
    sns.set_color_codes()
    ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.1))
    plt.tight_layout()
    ax.set(ylabel='Number of lncRNAs', xlabel='')
    sns.despine(ax=ax, right=True, top=True, left=True)
    plt.show()

    s = np.zeros((TanricDataset.n_genes, n_datasets), dtype=bool)
    cancer_list = []
    set_1_names = ('KIRC', 'KICH', 'KIRP')
    set_1_idxs = []
    set_2_names = ('HNSC', 'STAD', 'BLCA')
    set_2_idxs = []
    for i, ds in enumerate(datasets):
        _, _, s[:, i] = ds.results['t_test']
        ct = ds.cancer_type
        cancer_list.append(ct)
        if ct in set_1_names:
            set_1_idxs.append(i)
        if ct in set_2_names:
            set_2_idxs.append(i)

    # FIXME - scale to correct size for each subset
    fig, axs = plt.subplots(1, 2)
    set_names = [set_1_names, set_2_names]
    # for i, group in enumerate([set_1_idxs, set_2_idxs]):
        # venn_sets = [set(np.where(s[:, j])[0]) for j in group]
        # venn3(venn_sets, set_names[i], ax=axs[i])
    # plt.show()

    set_names = [
        tuple(),
        tuple(),
        ('HNSC', 'LUAD', 'STAD'),
        ('BRCA', 'PRAD'),
        ('LIHC', 'PRAD', 'STAD'),
        ('HNSC', 'LUAD', 'LUSC'),
        ('HNSC', 'LIHC', 'STAD'),
        tuple(),
        ('KICH', 'KIRP', 'THCA')
    ]
    set_idxs = [[] for _ in range(len(set_names))]
    for i, ds in enumerate(datasets):
        ct = ds.cancer_type
        cancer_list.append(ct)
        # for j, names in enumerate(set_names):
        #      if ct in names:
        #         set_idxs[j].append(i)

    # FIXME - scale to correct size for each subset
    # fig, axs = plt.subplots(3, 3)
    # for i, group in enumerate(set_idxs):
    #     if set_names[i]:
    #         venn_sets = [set(np.where(s[:, j])[0]) for j in group]
    #         if len(set_names[i]) == 3:
    #             venn3(venn_sets, set_names[i], ax=axs[i // 3][i % 3])
    #         if len(set_names[i]) == 2:
    #             venn2(venn_sets, set_names[i], ax=axs[i // 3][i % 3])
    # plt.show()

    num_signif = np.sum(s, 1)
    max_signif_idxs = np.where(num_signif >= 9)[0]
    max_signif_idxs = max_signif_idxs[np.argsort(num_signif[max_signif_idxs])]
    max_signif_idxs = np.flip(max_signif_idxs, 0)
    panda_dict = {}
    print(cancer_list)
    for idx in max_signif_idxs:
        check_str = ' | '.join(['+' if v else ' ' for v in s[idx, :]])
        g_name = bytes.decode(TanricDataset.gene_info['code'][idx], 'utf-8')
        panda_dict[g_name] = s[idx, :]
        print('%20s | %2d | %s'
              % (g_name,
                 np.count_nonzero(s[idx, :]),
                 check_str))

    # name_dict = {}
    # for i, n in enumerate(cancer_list):
    #    name_dict[i] = n

    # df2 = pd.DataFrame(panda_dict)
    # df2 = df2.rename(name_dict, axis='index')

    #sns.clustermap(df2)
    # plt.show()
    # Pie CHart of Overlay
    x = []
    nums = np.sort(np.unique(num_signif))
    labels = []
    for num in nums:
        if num > 1:
            count = np.count_nonzero(num_signif == num)
            x.append(count)
            labels.append('%d Cancers (%d - %.0f%%)'
                          % (num, count,
                             count / np.count_nonzero(num_signif > 1) * 100))

    fig, axs = plt.subplots(1, 2, figsize=(14, 7))
    tot = TanricDataset.n_genes
    x2 = [np.count_nonzero(num_signif == 0),
          np.count_nonzero(num_signif == 1),
          np.count_nonzero(num_signif > 1)]
    axs[0].pie(x2,
               labels=['Not Significant\n(%d - %.0f%%)'
                       % (x2[0], x2[0]/tot*100),
                       'Significant in\n1 Cancer\n(%d - %.0f%%)'
                       % (x2[1], x2[1]/tot*100),
                       'Significant in\n2 or more Cancers\n(%d - %.0f%%)'
                       % (x2[2], x2[2]/tot*100)])
    _, texts = axs[1].pie(x, labels=labels)
    for i, t in enumerate(texts):
        x, y = t.get_position()
        if (i + 2) == 11:
            t.set_position((x, y + 0.1))
        if (i + 2) == 10:
            t.set_position((x, y + 0.05))

    plt.show()

    print('\n\nTable')
    for ds in datasets:
        code = ds.cancer_type
        full_name = TanricDataset.sampleid2name(code)
        nn = ds.n_normal_samples
        nt = ds.n_tumor_samples
        ne = np.count_nonzero(ds.results['is_expressed'])
        _, _, is_signif = ds.results['t_test']
        ns = np.count_nonzero(is_signif)
        fmt = '%s | %38s | %3d | %3d | %4d / %4d = %.0f%%'
        print(fmt % (code, full_name, nn, nt, ns, ne, 100 * ns / ne))
    # tabls of top overlapping genes


def make_random_plots(data):
    values = data[0]
    fignum = 0
    fignum += 1
    plt.figure(fignum)
    plt.hist(values.reshape(-1), bins=200)
    plt.show()


def make_ma_plots(datasets, fcf):
    fignum = 0
    n_datasets = len(datasets)
    n_points = n_datasets * TanricDataset.n_genes
    fc = np.zeros(n_points)
    mean = np.zeros(n_points)
    numsignif = np.zeros(TanricDataset.n_genes)
    signif = np.zeros(n_points)
    signif2 = np.zeros((TanricDataset.n_genes, n_datasets), dtype=int)
    valid = np.zeros(n_points, dtype=bool)

    fignum += 1
    fig = plt.figure(fignum, (15, 10))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(4, 3),
                           axes_pad=0.25, label_mode='L',)
    for i, ds in enumerate(datasets):
        _, _, is_signif = ds.results['t_test']
        is_valid = ds.results['is_nonzero']
        norm_mean = np.mean(ds.normal_samples, 1) + fcf
        tumor_mean = np.mean(ds.tumor_samples, 1) + fcf

        selec = np.arange(TanricDataset.n_genes, dtype=int) + \
            TanricDataset.n_genes * i
        fc[selec] = np.log2(tumor_mean / norm_mean)
        valid[selec] = is_valid
        mean[selec] = tumor_mean
        signif[selec] = is_signif
        signif2[:, i] = is_signif
        numsignif += is_signif.astype(int)

        sz = 4
        ax = grid[i]
        x = tumor_mean
        y = fc[selec]
        ns = np.logical_and(np.logical_not(is_signif), is_valid)
        ax.scatter(x[ns], y[ns], s=sz**2, c='k', alpha=0.5)
        ax.scatter(x[is_signif], y[is_signif], s=sz**2, c='r', alpha=0.8)
        ax.text(0.99, 0.980,
                '%s\n%.1f%% (%d/%d) significant'
                % (TanricDataset.sampleid2name(ds.cancer_type),
                   np.count_nonzero(is_signif) / np.count_nonzero(is_valid) *
                   100, np.count_nonzero(is_signif),
                   np.count_nonzero(is_valid)),
                transform=ax.transAxes,
                family='Lao Sangam MN',
                fontsize=10,
                va='top',
                ha='center',
                bbox=(dict(boxstyle='square',
                           facecolor='white',
                           ec='k',
                           alpha=1,
                           lw=0.75)))
        ax.plot(np.array([-1000, 1000]), np.array([0, 0]), 'k',
                linewidth='0.75',
                zorder=0, alpha=0.5)
        ax.set_xlabel('Mean')
        ax.set_ylabel('Log Fold Change')
        ax.set_xlim(-0.1, 15)
        ax.set_ylim(-15, 15)

    plt.tight_layout()
    plt.show()

    gene_nums = np.tile(np.arange(TanricDataset.n_genes, dtype=int), n_datasets)
    numsignif = np.tile(numsignif, n_datasets)
    numsignif *= signif.astype(int)

    store = (valid, numsignif, mean, fc)

    numsignif = numsignif[valid]
    fc = fc[valid]
    mean = mean[valid]
    gene_nums = gene_nums[valid]

    fignum += 1
    fig = plt.figure(fignum, (15, 10))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(4, 3),
                           axes_pad=0.25, label_mode='L',)

    n_valid = np.count_nonzero(valid)
    for i in range(12):
        signif_local = (numsignif == (i + 1))
        notsignif = (numsignif == 0)
        n_signif = np.count_nonzero(signif_local)
        n_signif_genes = int(n_signif / (i+1))
        ax = grid[i]
        ax.scatter(mean[notsignif], fc[notsignif], s=sz**2, c='k', alpha=0.3)
        if 0 < n_signif_genes < 16:
            gn = gene_nums[signif_local]
            for j, idx in enumerate(np.unique(gn)):
                selec = np.where(np.logical_and(gene_nums == idx,
                                                signif_local))[0]
                n = bytes.decode(TanricDataset.gene_info['code'][idx],
                                 'utf-8')
                ax.scatter(mean[selec], fc[selec], s=sz**2, alpha=1, label=n)
            ax.legend(facecolor='white', framealpha=0.8, loc='lower right',
                      ncol=2)
        else:
            ax.scatter(mean[signif_local], fc[signif_local], s=sz**2, c='r',
                       alpha=1)
        ax.plot(np.array([-1000, 1000]), np.array([0, 0]), 'k',
                linewidth='0.75',
                zorder=0, alpha=0.5)
        ax.text(0.99, 0.980,
                'Present in %d Cancers\n %5.1f%% (%d/%d) significant'
                % (i + 1, n_signif / n_valid * 100, n_signif // (i + 1),
                   n_valid // n_datasets),
                transform=ax.transAxes,
                family='Lao Sangam MN',
                fontsize=10,
                va='top',
                ha='right',
                bbox=(dict(boxstyle='square',
                           facecolor='white',
                           ec='k',
                           alpha=1,
                           lw=0.75)))
        ax.set_xlabel('Mean')
        ax.set_ylabel('Log Fold Change')
        # ax.set_xlim(-0.1, max(fc[numsignif > 0]) + 0.5)
        ax.set_xlim(-0.1, 15)
        ax.set_ylim(-15, 15)

    plt.tight_layout()
    plt.show()

    valid, numsignif, mean, fc = store

    d = {}
    for i, ds in enumerate(datasets):
        selec = np.arange(TanricDataset.n_genes, dtype=int) + \
                TanricDataset.n_genes * i
        fc_signif = fc[selec]
        _, _, is_signif = ds.results['t_test']
        sig = np.logical_and(is_signif, valid[selec])
        fc_signif = fc_signif[sig]
        name = ds.cancer_type
        name = '\n'.join(textwrap.wrap(name, 15))
        d[name] = pd.Series(fc_signif)

    df = pd.DataFrame(d)
    means = df.mean()
    df = df[df.columns[means.argsort()]]
    fignum += 1
    fig = plt.figure(fignum, (16, 8))
    ax = plt.subplot(111)
    sns.violinplot(data=df, ax=ax, orient='v', showfliers=False, bw=0.15,
                   linewidth=1, inner=None, cut=0, width=1)
    sns.swarmplot(data=df, ax=ax, color='k', size=1.5, orient='v', alpha=0.2)
    a = np.array(ax.get_xlim())
    b = np.array([0, 0])
    plt.plot(a, b, 'k', lw=0.75, alpha=0.5, zorder=0)
    plt.ylim(-13, 13)
    plt.ylabel('Log Fold Change')
    plt.tight_layout()
    # fig_path = os.path.join('figures', 'cancer_diff_exp_distribution.png')
    # plt.savefig(fig_path, dpi=300)
    plt.show()


def make_volcano_plots(datasets, test, fcf):
    fignum = 10
    n_datasets = len(datasets)
    n_points = n_datasets * TanricDataset.n_genes
    fc = np.zeros(n_points)
    qls = np.zeros(n_points)
    pls = np.zeros(n_points)
    mean = np.zeros(n_points)
    numsignif = np.zeros(TanricDataset.n_genes)
    signif = np.zeros(n_points)
    signif2 = np.zeros((TanricDataset.n_genes, n_datasets), dtype=int)
    valid = np.zeros(n_points, dtype=bool)

    fignum += 1
    fig = plt.figure(fignum, (15, 10))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(2, 5),
                           axes_pad=0.25, label_mode='L',)
    for i, ds in enumerate(datasets):
        t_stat, p, is_signif = ds.results['t_test']
        q = ds.results['q_values']
        is_valid = ds.results['is_expressed']


        selec = np.arange(TanricDataset.n_genes, dtype=int) + \
            TanricDataset.n_genes * i


        if test is 'wsr':
            nidx, tidx = ds.sample_pairs
            nsps = ds.exprdata[:, nidx]
            tsps = ds.exprdata[:, tidx]
            tumor_mean = np.mean(tsps, 1)
            fc_temp = (tsps + fcf) / (nsps + fcf)
            # fc[selec] = np.log2(np.mean(fc_temp, 1))
            fc[selec] = np.log2((np.mean(tsps, 1) + fcf) /
                                (np.mean(nsps, 1) + fcf))
            # fc[selec] = t_stat * np.sign(np.mean(tsps, 1) - np.mean(nsps, 1))
        else:
            # FIXME, - this is an improper calculation
            tumor_mean = np.mean(ds.tumor_samples, 1) + fcf
            norm_mean = np.mean(ds.normal_samples, 1) + fcf
            fc[selec] = np.log2(tumor_mean / norm_mean)


        q[is_valid] = -np.log10(q[is_valid])
        qls[selec] = q
        p[is_valid] = -np.log10(p[is_valid])
        pls[selec] = p
        valid[selec] = is_valid
        mean[selec] = tumor_mean
        signif[selec] = is_signif
        signif2[:, i] = is_signif
        numsignif += is_signif.astype(int)

        sz = 3.5
        ax = grid[i]
        y = qls[selec]
        # y = pls[selec]
        x = fc[selec]
        ns = np.logical_and(np.logical_not(is_signif), is_valid)
        ax.scatter(x[ns], y[ns], s=sz**2, c='k', alpha=0.8)
        ax.scatter(x[is_signif], y[is_signif], s=sz**2, c='r', alpha=1)
        ax.text(0.99, 0.980,
                '%s\n%.1f%% (%d/%d) significant'
                % (TanricDataset.sampleid2name(ds.cancer_type),
                   np.count_nonzero(is_signif) / np.count_nonzero(is_valid) *
                   100, np.count_nonzero(is_signif),
                   np.count_nonzero(is_valid)),
                transform=ax.transAxes,
                family='Lao Sangam MN',
                fontsize=8,
                va='top',
                ha='right',
                bbox=(dict(boxstyle='square',
                           facecolor='white',
                           ec='k',
                           alpha=1,
                           lw=0.75)))
        ax.plot(np.array([0, 0]), np.array([0, 20]), 'k',
                linewidth='0.75',
                zorder=0, alpha=0.5)
        ax.set_ylabel('-Log10 q-Value')
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylim(-0.1, 18)
        #ax.set_xlim(-15, 15)

    plt.tight_layout()
    fig_path = os.path.join('figures', 'volcano_by_cancer_scramble.png')
    plt.savefig(fig_path, dpi=300)
    plt.show()

    gene_nums = np.tile(np.arange(TanricDataset.n_genes, dtype=int), n_datasets)
    numsignif = np.tile(numsignif, n_datasets)
    numsignif *= signif.astype(int)

    store = (valid, numsignif, mean, fc)

    numsignif = numsignif[valid]
    fc = fc[valid]
    qls = qls[valid]
    pls = pls[valid]
    mean = mean[valid]
    gene_nums = gene_nums[valid]

    fignum += 1
    fig = plt.figure(fignum, (15, 10))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(2, 5),
                           axes_pad=0.25, label_mode='L',)

    n_valid = np.count_nonzero(valid)
    for i in range(n_datasets):
        signif_local = (numsignif == (i + 1))
        notsignif = numsignif == 0
        n_signif = np.count_nonzero(signif_local)
        n_signif_genes = int(n_signif / (i+1))
        ax = grid[i]
        # FIXME - var names should be better
        aa = qls[notsignif]
        # aa = pls[notsignif]
        bb = fc[notsignif]
        ax.scatter(bb, aa, s=sz**2, c='k', alpha=0.8)
        if 0 < n_signif_genes < 10:
            gn = gene_nums[signif_local]
            for j, idx in enumerate(np.unique(gn)):
                selec = np.where(np.logical_and(gene_nums == idx,
                                                signif_local))[0]
                n = bytes.decode(TanricDataset.gene_info['code'][idx],
                                 'utf-8')
                a2 = qls[selec]
                # a2 = pls[selec]
                b2 = fc[selec]
                ax.scatter(b2, a2, s=sz**2, alpha=1, label=n)
            ax.legend(facecolor='white', framealpha=0.8, loc='lower right',
                      ncol=1)
        else:
            a2 = qls[signif_local]
            # a2 = pls[signif_local]
            b2 = fc[signif_local]
            ax.scatter(b2, a2, s=sz**2, c='r',
                       alpha=1)
        ax.plot(np.array([0, 0]), np.array([0, 20]), 'k',
                linewidth='0.75',
                zorder=0, alpha=0.5)
        ax.text(0.99, 0.980,
                'Present in %d Cancers\n %5.1f%% (%d/%d) significant'
                % (i + 1, n_signif / n_valid * 100, n_signif // (i + 1),
                   n_valid // n_datasets),
                transform=ax.transAxes,
                family='Lao Sangam MN',
                fontsize=8,
                va='top',
                ha='right',
                bbox=(dict(boxstyle='square',
                           facecolor='white',
                           ec='k',
                           alpha=1,
                           lw=0.75)))
        ax.set_ylabel('Log10 q-value')
        ax.set_xlabel('Log2 Fold Change')
        # ax.set_xlim(-0.1, max(fc[numsignif > 0]) + 0.5)
        ax.set_ylim(-0.1, 18)
        #ax.set_xlim(-15, 15)

    plt.tight_layout()
    fig_path = os.path.join('figures', 'volcano_by_number_scramble.png')
    plt.savefig(fig_path, dpi=300)
    plt.show()

    valid, numsignif, mean, fc = store

    d = {}
    for i, ds in enumerate(datasets):
        selec = np.arange(TanricDataset.n_genes, dtype=int) + \
                TanricDataset.n_genes * i
        fc_signif = fc[selec]
        _, _, is_signif = ds.results['t_test']
        sig = np.logical_and(is_signif, valid[selec])
        fc_signif = fc_signif[sig]
        name = ds.cancer_type
        name = '\n'.join(textwrap.wrap(name, 15))
        d[name] = pd.Series(fc_signif)

    df = pd.DataFrame(d)
    means = df.mean()
    df = df[df.columns[means.argsort()]]
    fignum += 1
    fig = plt.figure(fignum, (16, 8))
    ax = plt.subplot(111)
    sns.violinplot(data=df, ax=ax, orient='v', showfliers=False, bw=0.15,
                   linewidth=1, inner=None, cut=0, width=1)
    sns.swarmplot(data=df, ax=ax, color='k', size=1.5, orient='v', alpha=0.2)
    a = np.array(ax.get_xlim())
    b = np.array([0, 0])
    plt.plot(a, b, 'k', lw=0.75, alpha=0.5, zorder=0)
    # plt.ylim(-13, 13)
    plt.ylabel('Log Fold Change')
    plt.tight_layout()
    fig_path = os.path.join('figures',
                            'cancer_diff_exp_distribution_scramble.png')
    plt.savefig(fig_path, dpi=300)
    plt.show()


