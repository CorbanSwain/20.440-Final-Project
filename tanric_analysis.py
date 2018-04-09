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

        cachefile = os.path.join(tanric_dir, 'np_cache', name + '.npy')
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


def perform_t_test(datasets, expr_cutoff, procedure, **kwargs):
    print('\nPerforming t-tests ...')
    # FIXME - May need to do some memory management here
    # FIXME - Validity testing should be done in its own function
    find_signif = mh_tests[procedure]
    n_counts = np.zeros((TanricDataset.n_genes,), dtype=int)
    all_valid = np.zeros((TanricDataset.n_genes,), dtype=bool)
    for ds in datasets:
        is_valid = np.mean(ds.exprdata, 1) > expr_cutoff
        all_valid = np.logical_or(all_valid, is_valid)
        norm_valid = ds.normal_samples[is_valid]
        tumor_valid = ds.tumor_samples[is_valid]

        t = np.zeros((ds.n_genes,))
        p = np.zeros((ds.n_genes,))
        t[is_valid], p[is_valid] = ttest_ind(norm_valid,
                                             tumor_valid, axis=1)

        is_signif = np.zeros((ds.n_genes,), dtype=bool)
        is_signif_valid = find_signif(p[is_valid], **kwargs)
        is_signif[is_valid] = is_signif_valid
        n_counts[is_signif] += 1
        stdout.write('\r\t%s: # implicated = %d'
                     % (ds.cancer_type, np.count_nonzero(is_signif)))
        stdout.flush()

        # fold change calculation
        fc = np.zeros((ds.n_genes,))
        # FIXME - parametrize fudge factor
        fc[is_valid] = (np.mean(tumor_valid, 1) + 1e-3) \
            / (np.mean(norm_valid, 1) + 1e-3)

        ds.results['t_test'] = (t, p, fc, is_valid, is_signif)
    stdout.write('\r\tDone.\n')

    print('\n\tCount Summary')
    print('\t\t%4d lncRNAs with significant expression in at least one '
          'cancer type.\n'
          % np.count_nonzero(all_valid))
    for i in range(1, max(n_counts) + 1):
        n = np.count_nonzero(n_counts == i)
        print('\t\t%4d lncRNAs implicated in %2d cancer types.' % (n, i))


def make_composite_dataset(datasets, expr_cutoff, pool_norm,
                           fold_change_fudge):
    print('\nMaking Composite Dataset ...')
    spinner.start()
    n_cancer_samples = 0
    n_normal_samples = 0

    # all_valid = np.ones((TanricDataset.n_genes, ), dtype=bool)
    expr_sum = np.zeros((TanricDataset.n_genes, ))

    for ds in datasets:
        n_cancer_samples += ds.n_tumor_samples
        # FIXME - might have too many ''if pool_norm:'' statements
        if pool_norm:
            n_normal_samples += ds.n_normal_samples
        expr_sum += np.sum(ds.tumor_samples, 1)
        # _, _, _, is_valid, _ = ds.results['t_test']
        # all_valid = np.logical_and(is_valid, all_valid)

    all_valid = (expr_sum / n_cancer_samples) >= expr_cutoff

    comb_genes = TanricDataset.gene_info['code'][all_valid]
    n_genes = np.count_nonzero(all_valid)
    comp_set = np.zeros((n_genes, n_cancer_samples))
    group_labels = np.zeros((n_cancer_samples, ), dtype=np.object)
    group_numbers = np.zeros((n_cancer_samples, ), dtype=int)
    insert_idx = 0

    if pool_norm:
        norm = np.zeros((n_genes, ))
        for ds in datasets:
            norm += np.sum(ds.normal_samples[all_valid], 1)
        norm /= n_normal_samples

    for i, ds in enumerate(datasets):
        finish = insert_idx + ds.n_tumor_samples
        insert_range = np.arange(insert_idx, finish, dtype=int)
        insert_idx = finish
        # FIXME - parametrize fudge factor
        # DONE - Try pooled normal samples
        if pool_norm:
            norm_expr_val = norm + fold_change_fudge
        else:
            norm_expr_val = np.mean(ds.normal_samples[all_valid],
                                    1) + fold_change_fudge
        expr_vals = ds.tumor_samples[all_valid] + fold_change_fudge
        fold_change = expr_vals.T / norm_expr_val
        comp_set[:, insert_range] = np.log2(fold_change.T)
        group_labels[insert_range] = ds.cancer_type
        group_numbers[insert_range] = i
    spinner.stop()
    print('\tDone.')
    return (comp_set, group_labels, group_numbers, comb_genes)


def save_for_matlab(datasets, comp_ds, settings):
    print('\nSaving results for MATLAB ...')
    spinner.start()
    n_sets = len(datasets)
    matlab_struct = {}
    cell_arr = lambda x: np.array(x, dtype=np.object) # TODO - add to utils
    col_vec = lambda x: x.reshape((-1, 1))  # TODO - add to utils
    row_vec = lambda x: x.reshape((1, -1))  # TODO - add to utils
    for ds in datasets:
        t, p, fc, is_valid, is_signif = ds.results['t_test']
        logfc = np.zeros(fc.shape)
        logfc[is_valid] = np.log2(fc[is_valid])
        sample_names = cell_arr(ds.sample_names)
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
    gene_ids = cell_arr(TanricDataset.gene_ids)

    # Using ... .copy(order='C') to ensure the slices are contiguous; see
    # https://stackoverflow.com/questions/26778079/valueerror-ndarray-is-not-c-contiguous-in-cython
    gene_codes = TanricDataset.gene_info['code'].copy(order='C')
    gene_codes = cell_arr(gene_codes)
    gene_descriptions = TanricDataset.gene_info['description'].copy(order='C')
    gene_descriptions = cell_arr(gene_descriptions)

    comp_set, labels, nums, comb_genes = comp_ds
    comb_genes = col_vec(cell_arr(comb_genes.copy(order='C')))

    matpath = os.path.join('data', 'matlab_io',
                           'part_1_analysis_v%s.mat' % settings['version'])
    spinner.stop()
    print('\t' + matpath)
    spinner.start()
    sio.savemat(matpath, {'S': matlab_struct,
                          'nGenes': n_genes,
                          'geneIDs': col_vec(gene_ids),
                          'geneCodes': col_vec(gene_codes),
                          'geneDescriptions': col_vec(gene_descriptions),
                          'analysisMetadata': settings,
                          'combGenes': comb_genes,
                          'combData': comp_set,
                          'combLabels': labels,
                          'combNumIDs': nums})
    spinner.stop()
    print('\tDone.')


if __name__ == "__main__":
    settings = {
        'min_norm_samples': 5,
        'version': '2.2',
        'expression_cutoff': 0.05,
        'multi_hyp_procedure': 'bonferoni',
        'alpha_crit': 0.05,
        'pool_normal_samples': True,
        'fold_change_fudge': 5e-4,
        'analysis_date': str(datetime.datetime.now())
    }

    datasets = import_all_data(settings['min_norm_samples'])
    TanricDataset.get_gene_info()

    perform_t_test(datasets,
                   expr_cutoff=settings['expression_cutoff'],
                   procedure=settings['multi_hyp_procedure'],
                   a=settings['alpha_crit'])

    comp_ds = make_composite_dataset(datasets,
                                     settings['expression_cutoff'],
                                     settings['pool_normal_samples'],
                                     settings['fold_change_fudge'])

    save_for_matlab(datasets, comp_ds, settings)
