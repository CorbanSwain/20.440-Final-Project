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
        ds.results['is_nonzero'] = np.mean(ds.exprdata, 1) > 1e-6
        ds.results['is_expressed'] = np.mean(ds.exprdata, 1) > expr_cutoff


def perform_t_test(datasets, procedure, **kwargs):
    print('\nPerforming t-tests ...')
    # FIXME - May need to do some memory management here
    # FIXME - Validity testing should be done in its own function
    find_signif = mh_tests[procedure]
    n_counts = np.zeros((TanricDataset.n_genes,), dtype=int)
    all_valid = np.zeros((TanricDataset.n_genes,), dtype=bool)
    for ds in datasets:
        is_valid = ds.results['is_expressed']
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
    print('\nMaking Composite Dataset ...')
    spinner.start()
    assert metric == 'fc_pair' and samples == 'tumor'
    assert metric == 'mean2mean' and samples == 'tumor'
    t_test_filter = filter_method == 't_test'

    # Filtering and counting
    n_datasets = len(datasets)
    n_cancer_samples = 0
    n_normal_samples = 0
    n_pairs = 0
    all_valid = np.ones(TanricDataset.n_genes, dtype=bool)
    if t_test_filter:
        any_signif = np.zeros(TanricDataset.n_genes, dtype=bool)
    for ds in datasets:
        n_cancer_samples += ds.n_tumor_samples
        n_normal_samples += ds.n_normal_samples
        n_pairs += ds.n_pairs
        if filter_method == 'none':
            is_valid = ds.results['is_nonzero']
        else:
            is_valid = ds.results['is_expressed']
        all_valid = np.logical_and(all_valid, is_valid)
        if t_test_filter:
            _, _, _, is_signif = ds.results['t_test']
            any_signif = np.logical_or(any_signif, is_signif)
    if t_test_filter:
        all_valid = np.logical_and(all_valid, any_signif)

    # setting up combined dataset
    comb_genes = TanricDataset.gene_info['code'].copy(order='C')
    comb_genes = comb_genes[all_valid]
    n_genes = np.count_nonzero(all_valid)
    spinner.stop()
    print('\t%d Genes considered significant' % n_genes)
    spinner.start()

    if metric == 'fc_pair':
        n_samples = n_pairs
    elif metric == 'mean2mean':
        n_samples = n_datasets
    else:
        if samples == 'normal':
            n_samples = n_normal_samples
        elif samples == 'tumor':
            n_samples = n_cancer_samples
        elif samples == 'all':
            n_samples = n_cancer_samples + n_normal_samples
        else:
            raise ValueError('Value for sample was unexpected.')

    combined_set = np.zeros((n_genes, n_samples))
    group_labels = np.zeros(n_samples, dtype='|S20')
    label_list = []
    group_numbers = np.zeros(n_samples, dtype=int)
    insert_idx = 0

    for i, ds in enumerate(datasets):
        if metric == 'fc_pair':
            end_idx = insert_idx + ds.n_pairs
        elif metric == 'mean2mean':
            end_idx = insert_idx + 1
        else:
            if samples == 'normal':
                end_idx = insert_idx + ds.n_normal_samples
            elif samples == 'tumor':
                end_idx = insert_idx + ds.n_tumor_samples
            else:
                end_idx = insert_idx + ds.n_samples

        insert_range = np.arange(insert_idx, end_idx, dtype=int)
        insert_idx = end_idx

        if metric == 'rpkm':
            value_matrix = ds.exprdata[all_valid]
        elif metric == 'mean2mean' or metric == 'fc_mean':
            norm_expr_val = np.mean(ds.normal_samples[all_valid], 1)
            norm_expr_val += fold_change_fudge
            value_matrix = ((ds.exprdata[all_valid].T +
                            fold_change_fudge) / norm_expr_val).T
            value_matrix = np.log2(value_matrix)
        elif metric = 'fc_pair':
            norm_sel, tumor_sel = ds.sample_pairs
            # FIXME - this could be implemented better, doing excess operations
            normal = ds.exprdata[:, norm_sel] + fold_change_fudge
            tumor = ds.exprdata[:, tumor_sel] + fold_change_fudge
            value_matrix = tumor[all_valid] / normal[all_valid]
            value_matrix = np.log2(value_matrix)

        for idx in insert_range:
            group_labels[idx] = ds.cancer_type
        group_numbers[insert_range] = i
    spinner.stop()
    print('\tDone.')
    return (combined_set, group_labels, group_numbers, comb_genes)


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


if __name__ == "__main__":
    settings = {
        'min_norm_samples': 5,
        'version': '3.1',  # TODO - some type of automatic versioning
        'expression_cutoff': 0.2,  # 0.3 used in TANRIC paper
        'filter_method': 't_test',  # t_test, threshold, none
        'multi_hyp_procedure': 'bonferoni',  #bonferoni, crit, ben_hoch
        'alpha_crit': 0.05,
        'metric': 'fc',  # fc_mean, fc_pair, rpkm, mean2mean
        'samples': 'tumor',  # tumor, normal, all
        'fold_change_fudge': 5e-4,
        'analysis_date': str(datetime.datetime.now())
    }

    datasets = import_all_data(settings['min_norm_samples'])
    TanricDataset.get_gene_info()

    assess_validity(datasets,
                    settings['expression_cutoff'])

    perform_t_test(datasets,
                   settings['multi_hyp_procedure'],
                   a=settings['alpha_crit'])

    comp_ds = make_composite_dataset(datasets,
                                     settings['expression_cutoff'],
                                     settings['filter_method'],
                                     settings['metric'],
                                     settings['samples'],
                                     settings['fold_change_fudge'])

    save_for_matlab(datasets, comp_ds, settings)
