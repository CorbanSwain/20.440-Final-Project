#!/bin/python3
# utils.py
# Corban Swain, 2018

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from utils import *
from scipy.stats import ttest_ind


def import_all_data(min_normal_samples):
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
        print('\tLoaded in dataset for %s (%3d, %3d)'
              % (tds.cancer_type, tds.n_normal_samples, tds.n_tumor_samples))
        datasets.append(tds)
    return datasets


def perform_t_test(datasets, expr_cutoff=0.3, procedure='crit', **kwargs):
    # FIXME - May need to do some memory management here
    find_signif = mh_tests[procedure]
    n_counts = np.zeros((datasets[0].n_genes,), dtype=int)
    for ds in datasets:
        is_valid = np.mean(ds.exprdata, 1) > expr_cutoff
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
        print('\t%s: # implicated = %d'
              % (ds.cancer_type, np.count_nonzero(is_signif)))

        # fold change calculation
        fc = np.zeros((ds.n_genes,))
        fc[is_valid] = (np.mean(tumor_valid, 1) + 1e-3) \
            / (np.mean(norm_valid, 1) + 1e-3)

        ds.results['t_test'] = (t, p, fc, is_valid, is_signif)

    print('\n\tCount Summary')
    for i in range(1, max(n_counts) + 1):
        n = np.count_nonzero(n_counts == i)
        print('\t\t%4d lncRNAs implicated in %2d cancer types.' % (n, i))


def save_for_matlab(datasets, version):
    n_sets = len(datasets)
    matlab_struct = {}
    col_vec = lambda x: x.reshape((-1, 1))
    row_vec = lambda x: x.reshape((1, -1))
    for ds in datasets:
        t, p, fc, is_valid, is_signif = ds.results['t_test']
        logfc = np.zeros(fc.shape)
        logfc[is_valid] = np.log2(fc[is_valid])
        sample_names = np.array(ds.sample_names, dtype=np.object)
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
    matpath = os.path.join('data', 'matlab_io',
                           'part_1_analysis_v%3.1f.mat' % version)
    n_genes = datasets[0].n_genes
    gene_ids = np.array(datasets[0].gene_ids, dtype=np.object)
    sio.savemat(matpath, {'S': matlab_struct, 'nGenes': n_genes,
                          'geneIDs': col_vec(gene_ids)})


if __name__ == "__main__":
    min_normal_samples = 100
    version = 1.1

    print('\n1-Beginning Data Import')
    datasets = import_all_data(min_normal_samples)

    print('\n2-Performing t-tests')
    perform_t_test(datasets, expr_cutoff=0.1, procedure='bonferoni')

    print('\nFetching Names ...')
    names = TanricDataset.geneid2name(datasets[0].gene_ids)
    gnamepath = os.path.join('data', 'tanric_data', 'np_cache',
                             'gene_names.npy')
    np.save(gnamepath, names)

    # print('\n3-Saving \'.mat\' File')
    # save_for_matlab(datasets, version)
