#!/bin/python3
# tanric_analysis.py
# Corban Swain, 2018

import os
import numpy as np
import matplotlib.pyplot as plt
from sys import stdout
import datetime
import scipy.io as sio
from utils import *
from scipy.stats import ttest_ind
from mpl_toolkits import axes_grid1
import seaborn as sns
import pandas as pd
import textwrap


def import_all_data(min_normal_samples):
    print('\nPerforming Data Import ...')
    tanric_dir = os.path.join('data', 'tanric_data')
    names = []
    with open(os.path.join(tanric_dir, 'names.txt')) as names_file:
        for line in names_file.readlines():
            if not line.find('TCGA') is -1:
                names.append(line.strip(' \n'))

    datasets = []
    for name in names:  # TODO - this can be parallelized
        metafile = os.path.join(tanric_dir, name + '-rnaexpr-META.tsv')
        tds = TanricDataset(file2dict(metafile))

        # criteria for including datasets
        if tds.n_normal_samples < min_normal_samples:
            continue

        cachefile = os.path.join('data', 'np_cache', name + '.npy')
        try:
            data = np.load(cachefile)
        except FileNotFoundError:
            datafile = os.path.join(tanric_dir, name + '-rnaexpr.tsv')
            # this takes a while
            data = np.genfromtxt(datafile, delimiter='\t', dtype=None,
                                 encoding=None, names=True, deletechars='')
            np.save(cachefile, data)
        tds.parse_exprdata(data)
        stdout.write('\r\tLoaded in dataset for %s (%3d, %3d)'
                     % (tds.cancer_type, tds.n_normal_samples,
                     tds.n_tumor_samples))
        stdout.flush()
        datasets.append(tds)
    stdout.write('\r\tDone.\n')
    return datasets


def assess_validity(datasets, expr_cutoff):
    for ds in datasets:
        ds.results['is_nonzero'] = np.logical_and(
            np.any(ds.normal_samples > 1e-6, 1),
            np.any(ds.tumor_samples > 1e-6, 1))
        ds.results['is_expressed'] = np.mean(ds.exprdata, 1) > expr_cutoff


def perform_t_test(datasets, t_filter, procedure, **kwargs):
    print('\nPerforming t-tests ...')
    # FIXME - May need to do some memory management here
    # FIXME - Validity testing should be done in its own function
    find_signif = procedure.signif_func()
    n_counts = np.zeros((TanricDataset.n_genes,), dtype=int)
    all_valid = np.zeros((TanricDataset.n_genes,), dtype=bool)
    for ds in datasets:
        is_valid = ds.results[t_filter]
        all_valid = np.logical_or(all_valid, is_valid)
        norm_valid = ds.normal_samples[is_valid]
        tumor_valid = ds.tumor_samples[is_valid]

        t = np.zeros((ds.n_genes,))
        p = np.zeros((ds.n_genes,))
        t[is_valid], p[is_valid] = ttest_ind(tumor_valid, norm_valid, axis=1)

        is_signif = np.zeros((ds.n_genes,), dtype=bool)
        is_signif_valid = find_signif(p[is_valid], **kwargs)
        is_signif[is_valid] = is_signif_valid
        n_counts[is_signif] += 1
        stdout.write('\r\t%s: # implicated = %d'
                     % (ds.cancer_type, np.count_nonzero(is_signif)))
        stdout.flush()
        ds.results['t_test'] = (t, p, is_signif)
    stdout.write('\r\tDone.\n')

    print('\n\tCount Summary')
    print('\t\t%4d lncRNAs with significant expression in at least one '
          'cancer type.\n' % np.count_nonzero(all_valid))
    for i in range(1, max(n_counts) + 1):
        n = np.count_nonzero(n_counts == i)
        print('\t\t%4d lncRNAs implicated in %2d cancer types.' % (n, i))


def make_composite_dataset(datasets, filter_method, metric,
                           samples, fold_change_fudge):
    print('\tMaking Composite Dataset ...')

    # FIXME - There is a better way than all these if statements,
    # maybe implement a set of lambda functions and call them for all of the
    # operations.
    if metric in [Metric.FC_PAIR, Metric.MEAN2MEAN]:
        assert samples is Samples.TUMOR

    # Filtering and counting
    n_datasets = len(datasets)
    n_cancer_samples = 0
    n_normal_samples = 0
    n_pairs = 0
    all_valid = np.ones(TanricDataset.n_genes, dtype=bool)
    num_signif = np.zeros(TanricDataset.n_genes, dtype=int)
    if filter_method is Filter.T_TEST:
        any_signif = np.zeros(TanricDataset.n_genes, dtype=bool)

    num_annotations_v1 = np.zeros(TanricDataset.n_genes, dtype=int)
    num_annotations_v2 = np.zeros(TanricDataset.n_genes, dtype=int)
    for i_gene in range(TanricDataset.n_genes):
        num_annotations_v1[i_gene] += len(TanricDataset.transcripts[i_gene])
        # FIXME - need to try except here
        try:
            num_annotations_v2[i_gene] += max(
                ts.n_aliases for ts in TanricDataset.transcripts[i_gene])
        except ValueError:
            num_annotations_v2[i_gene] = 0

    for ds in datasets:
        n_cancer_samples += ds.n_tumor_samples
        n_normal_samples += ds.n_normal_samples
        n_pairs += ds.n_pairs
        _, _, is_signif = ds.results['t_test']
        num_signif += is_signif.astype(int)
        if filter_method is Filter.NONE:
            is_valid = ds.results['is_nonzero']
        elif filter_method in [Filter.THRESHOLD, Filter.T_TEST]:
            is_valid = ds.results['is_expressed']
        all_valid = np.logical_and(all_valid, is_valid)
        if filter_method is Filter.T_TEST:
            any_signif = np.logical_or(any_signif, is_signif)

    if filter_method is Filter.T_TEST:
        all_valid = np.logical_and(all_valid, any_signif)


    # setting up combined dataset
    genes_names = TanricDataset.gene_info['code'].copy(order='C')
    gene_descr = TanricDataset.gene_info['description'].copy(order='C')
    genes_names = genes_names[all_valid]
    gene_descr = gene_descr[all_valid]
    num_signif = num_signif[all_valid]
    num_annotations_v1 = num_annotations_v1[all_valid]
    num_annotations_v2 = num_annotations_v2[all_valid]
    n_genes = np.count_nonzero(all_valid)
    print('\t%d Genes considered significant' % n_genes)

    if metric is Metric.FC_PAIR:
        n_samples = n_pairs
    elif metric is Metric.MEAN2MEAN:
        n_samples = n_datasets
    else:
        if samples is Samples.NORMAL:
            n_samples = n_normal_samples
        elif samples is Samples.TUMOR:
            n_samples = n_cancer_samples
        elif samples is Samples.ALL:
            n_samples = n_cancer_samples + n_normal_samples

    combined_values = np.zeros((n_genes, n_samples), dtype=np.float64)
    group_labels = np.zeros(n_samples, dtype='|S20')
    label_list = []
    label_counter = 0
    group_numbers = np.zeros(n_samples, dtype=int)
    insert_idx = 0

    for ds in datasets:
        if metric == Metric.FC_PAIR:
            end_idx = insert_idx + ds.n_pairs
        elif metric == Metric.MEAN2MEAN:
            end_idx = insert_idx + 1
        else:
            if samples is Samples.NORMAL:
                end_idx = insert_idx + ds.n_normal_samples
            elif samples is Samples.TUMOR:
                end_idx = insert_idx + ds.n_tumor_samples
            elif samples is Samples.ALL:
                end_idx = insert_idx + ds.n_samples

        insert_range = np.arange(insert_idx, end_idx, dtype=int)
        insert_idx = end_idx

        if metric is Metric.RPKM:
            value_matrix = ds.exprdata[all_valid]
        elif metric in [Metric.MEAN2MEAN, Metric.FC_MEAN]:
            norm_expr_val = np.mean(ds.normal_samples[all_valid], 1)
            norm_expr_val += fold_change_fudge
            value_matrix = ((ds.exprdata[all_valid].T +
                            fold_change_fudge) / norm_expr_val).T
            value_matrix = np.log2(value_matrix)
        elif metric is Metric.FC_PAIR:
            norm_sel, tumor_sel = ds.sample_pairs
            # FIXME - this could be implemented better, doing excess operations
            normal = ds.exprdata[:, norm_sel] + fold_change_fudge
            tumor = ds.exprdata[:, tumor_sel] + fold_change_fudge
            value_matrix = tumor[all_valid, :] / normal[all_valid, :]
            value_matrix = np.log2(value_matrix)

        if samples is Samples.NORMAL:
            combined_values[:, insert_range] = value_matrix[:, ds.normal_sel]
        elif samples is Samples.TUMOR:
            if metric is Metric.MEAN2MEAN:
                combined_values[:, insert_range] \
                    = np.mean(value_matrix[:, ds.tumor_sel], 1, keepdims=True)
            elif metric in [Metric.FC_MEAN, Metric.RPKM]:
                combined_values[:, insert_range] = value_matrix[:, ds.tumor_sel]
            elif metric is Metric.FC_PAIR:
                combined_values[:, insert_range] = value_matrix
        elif samples is Samples.ALL:
            combined_values[:, insert_range] = value_matrix

        cancer_name = ds.cancer_type
        norm_name = 'Normal-' + cancer_name
        for i, idx in enumerate(insert_range):
            if samples is Samples.NORMAL:
                group_labels[idx] = norm_name
                group_numbers[idx] = label_counter
            elif samples is Samples.TUMOR:
                group_labels[idx] = cancer_name
                group_numbers[idx] = label_counter
            elif samples is Samples.ALL:
                if i < ds.n_normal_samples:
                    group_labels[idx] = norm_name
                    group_numbers[idx] = label_counter
                else:
                    group_labels[idx] = cancer_name
                    group_numbers[idx] = label_counter + 1

        if samples is Samples.ALL:
            label_counter = label_counter + 2
            label_list.append(norm_name)
            label_list.append(cancer_name)
        else:
            label_counter = label_counter + 1
            if samples is Samples.NORMAL:
                label_list.append(norm_name)
            elif samples is Samples.TUMOR:
                label_list.append(cancer_name)

    print('\tDone.')
    return (combined_values,
            group_labels,
            group_numbers,
            strlist2np(label_list),
            genes_names,
            gene_descr,
            num_signif,
            num_annotations_v1,
            num_annotations_v2,
            n_genes)


def make_random_plots(data):
    values = data[0]

    fignum = 0
    fignum += 1
    plt.figure(fignum)
    plt.hist(values.reshape(-1), bins=200)
    plt.show()


def make_ma_plots(datasets, t_filter, fcf):
    fignum = 0
    n_datasets = len(datasets)
    n_points = n_datasets * TanricDataset.n_genes
    fc = np.zeros(n_points)
    mean = np.zeros(n_points)
    numsignif = np.zeros(TanricDataset.n_genes)
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
        numsignif += is_signif.astype(int)

        sz = 4
        ax = grid[i]
        x = tumor_mean
        y = fc[selec]
        ns = np.logical_and(np.logical_not(is_signif), is_valid)
        ax.scatter(x[ns], y[ns], s=sz**2, c='k', alpha=0.5)
        ax.scatter(x[is_signif], y[is_signif], s=sz**2, c='r', alpha=0.8)
        ax.text(1, 1,
                '%s\n%.1f%% (%d) significant'
                % (TanricDataset.sampleid2name(ds.cancer_type),
                   np.count_nonzero(is_signif) / np.count_nonzero(is_valid) *
                   100, np.count_nonzero(is_signif)),
                transform=ax.transAxes,
                family='helvetica',
                fontsize=10,
                va='top',
                ha='right',
                bbox=(dict(boxstyle='square',
                           facecolor='white',
                           ec='k',
                           alpha=1,
                           lw=0.75)))
        ax.plot(np.array([-1000, 1000]), np.array([0, 0]), 'k',
                linewidth='0.75',
                zorder=0, alpha=0.5)
        ax.set_xlim(-0.1, 10)
        ax.set_ylim(-10, 10)

    plt.tight_layout()
    plt.show()

    numsignif = np.repeat(numsignif, n_datasets)

    store = (valid, numsignif, mean, fc)

    numsignif = numsignif[valid]
    fc = fc[valid]
    mean = mean[valid]

    fignum += 1
    fig = plt.figure(fignum, (15, 10))
    grid = axes_grid1.Grid(fig, rect=111, nrows_ncols=(4, 3),
                           axes_pad=0.25, label_mode='L',)

    for i in range(12):
        signif = numsignif == (i + 1)
        notsignif = np.logical_not(signif)

        ax = grid[i]
        ax.scatter(mean[notsignif], fc[notsignif], s=sz**2, c='k', alpha=0.3)
        ax.scatter(mean[signif], fc[signif], s=sz**2, c='r', alpha=0.8)
        ax.plot(np.array([-1000, 1000]), np.array([0, 0]), 'k',
                linewidth='0.75',
                zorder=0, alpha=0.5)
        ax.text(1, 1,
                'Present in %d Cancers' % (i + 1),
                transform=ax.transAxes,
                family='helvetica',
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
        ax.set_xlim(-0.1, 10)
        ax.set_ylim(-10, 10)


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
        name = TanricDataset.sampleid2name(ds.cancer_type)
        name = '\n'.join(textwrap.wrap(name, 15))
        d[name] = pd.Series(fc_signif)

    df = pd.DataFrame(d)
    means = df.median()
    df = df[df.columns[means.argsort()]]
    fignum += 1
    plt.figure(fignum, (10, 14))
    ax = plt.subplot(111)
    sns.boxplot(data=df, ax=ax, orient='h', showfliers=False)
    sns.swarmplot(data=df, ax=ax, color='k', size=2, orient='h')
    a = np.array(ax.get_ylim())
    b = np.array([0, 0])
    plt.plot(b, a, 'k', lw=0.75, alpha=0.5, zorder=0)
    plt.xlim(-6, 6)
    plt.xlabel('Log Fold Change')
    plt.tight_layout()
    plt.show()


def save_for_matlab(datasets, comp_ds, settings):
    print('\nSaving results for MATLAB ...')
    spinner.start()
    n_sets = len(datasets)
    matlab_struct = {}
    for ds in datasets:
        t, p, fc, is_valid, is_signif = ds.results['t_test']
        logfc = np.zeros(fc.shape)
        logfc[is_valid] = np.log2(fc[is_valid])
        sample_names = matlab_cell_arr(ds.sample_names)
        temp_dict = {'tStat': col_vec(t),
                     'pValue': col_vec(p),
                     'fc': col_vec(fc),
                     'logfc': col_vec(logfc),
                     'isValid': col_vec(is_valid),
                     'isSignif': col_vec(is_signif),
                     'isNormal': row_vec(ds.normal_sel),
                     'isTumor': row_vec(ds.tumor_sel),
                     'expression': ds.exprdata,
                     'nSamples': ds.n_samples,
                     'sampleNames': row_vec(sample_names)}
        matlab_struct[ds.cancer_type] = temp_dict

    n_genes = TanricDataset.n_genes
    gene_ids = matlab_cell_arr(TanricDataset.gene_ids)

    # Using ... .copy(order='C') to ensure the slices are contiguous; see
    # https://stackoverflow.com/questions/26778079/valueerror-ndarray-is-not
    # -c-contiguous-in-cython
    gene_codes = TanricDataset.gene_info['code'].copy(order='C')
    gene_codes = matlab_cell_arr(gene_codes)
    gene_descriptions = TanricDataset.gene_info['description'].copy(order='C')
    gene_descriptions = matlab_cell_arr(gene_descriptions)

    comp_set, labels, nums, comb_genes = comp_ds
    labels = matlab_cell_arr(labels)
    comb_genes = matlab_cell_arr(comb_genes)

    matpath = os.path.join('data', 'matlab_io', '%s_analysis_v%s.mat'
                           % ('%s', settings['version']))
    spinner.stop()
    print('\t' + matpath % 't_test')
    spinner.start()
    sio.savemat(matpath % 't_test',
                {'S': matlab_struct,
                 'nGenes': n_genes,
                 'geneIDs': col_vec(gene_ids),
                 'geneCodes': col_vec(gene_codes),
                 'geneDescriptions': col_vec(gene_descriptions),
                 'analysisMetadata': settings})
    spinner.stop()
    print('\t' + matpath % 'combined')
    spinner.start()
    sio.savemat(matpath % 'combined',
                {'analysisMetadata': settings,
                 'combGenes': col_vec(comb_genes),
                 'combData': comp_set,
                 'combLabels': labels,
                 'combNumIDs': nums})
    spinner.stop()
    print('\tDone.')


def save_for_matlab_2(filename, combined_data, settings):
    print('\tSaving results for MATLAB ...')
    spinner.start()

    ml_settings = {}
    for k, v in settings.items():
        if isinstance(v, Enum):
            ml_settings[k] = v.value
        else:
            ml_settings[k] = v

    values, g_labels, g_numbers, label_list, gene_names, gene_descr, \
        num_signif, n_anat1, n_anat2, n_genes = combined_data
    n_samples = len(g_labels)
    g_labels = matlab_cell_arr(g_labels)
    label_list = matlab_cell_arr(label_list)
    gene_names = matlab_cell_arr(gene_names)
    gene_descr = matlab_cell_arr(gene_descr)

    matpath = os.path.join('data', 'matlab_io',
                           'multi_analysis', filename)

    sio.savemat(matpath,
                {'analysisMetadata': ml_settings,
                 'values': values,
                 'numSamples': n_samples,
                 'sampleNames': row_vec(g_labels),
                 'sampleGroupNumbers': row_vec(g_numbers),
                 'sampleList': label_list,
                 'numGenes': n_genes,
                 'geneNames': col_vec(gene_names),
                 'geneDescrips': col_vec(gene_descr),
                 'geneNumSignif': col_vec(num_signif),
                 'geneNumAnat1': col_vec(n_anat1),
                 'geneNumAnat2': col_vec(n_anat2)})
    spinner.stop()
    print('\tDone.')


def make_multi_analysis(datasets, settings):
    file_list = []

    for filt in Filter:
        for metric in Metric:
            for samples in Samples:
                print('\nANALYSIS - %s - %s - %s' % (filt, metric, samples))

                settings['filter_method'] = filt
                settings['metric'] = metric
                settings['samples'] = samples

                try:
                    data = make_composite_dataset(
                        datasets,
                        filt,
                        metric,
                        samples,
                        settings['fold_change_fudge'])
                except AssertionError:
                    print('\tEXCEPTION: Couldn\'t Process')
                    continue

                if settings['do_save']:
                    filename = '%s-%s-%s_v%s.mat' % (filt,
                                                     metric,
                                                     samples,
                                                     settings['version'])

                    save_for_matlab_2(filename, data, settings)
                    file_list.append(filename)

    if settings['do_save']:
        path = os.path.join('data', 'matlab_io',
                            'multi_analysis', 'file_list.mat')
        sio.savemat(path, {'fileNames': matlab_cell_arr(file_list)})


if __name__ == "__main__":
    settings = {
        'min_norm_samples': 5,
        'version': '5.0',  # TODO - some type of automatic versioning
        'expression_cutoff': 0.3,  # 0.3 used in TANRIC paper
        'filter_method': None,
        't_filter': 'is_expressed', # is_nonzero, is_expressed
        'multi_hyp_procedure': MultiHypProc.BONFERONI,
        'alpha_crit': 0.005,
        'metric': None,
        'samples': None,
        'fold_change_fudge': 5e-3,
        'do_save': True,
        'do_plot': False,
        'analysis_date': str(datetime.datetime.now())
    }

    if settings['multi_hyp_procedure'] is MultiHypProc.BEN_HOCH:
        add_args = {'q': settings['alpha_crit']}
    else:
        add_args = {'a': settings['alpha_crit']}

    datasets = import_all_data(settings['min_norm_samples'])

    TanricDataset.get_transcript_info()

    TanricDataset.get_gene_info()

    assess_validity(datasets, settings['expression_cutoff'])

    perform_t_test(datasets,
                   settings['t_filter'],
                   settings['multi_hyp_procedure'],
                   **add_args)

    make_multi_analysis(datasets, settings)

    plt.style.use('seaborn-notebook')
    make_ma_plots(datasets,
                  settings['t_filter'],
                  settings['fold_change_fudge'])



