#!/bin/python3
# utils.py
# Corban Swain, 2018

import os
import numpy as np
import matplotlib.pyplot as plt
from utils import *
from scipy.stats import ttest_ind

MIN_NORMAL_SAMPLES = 5


def import_all_data():
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
        if tds.n_normal_samples < MIN_NORMAL_SAMPLES:
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
        valid_genes = np.mean(ds.exprdata, 1) > expr_cutoff
        norm_valid = ds.normal_samples[valid_genes]
        tumor_valid = ds.tumor_samples[valid_genes]

        t = np.zeros((ds.n_genes,))
        p = np.zeros((ds.n_genes,))
        t[valid_genes], p[valid_genes] = ttest_ind(norm_valid,
                                                   tumor_valid, axis=1)

        is_signif = np.zeros((ds.n_genes, ), dtype=bool)
        is_signif_valid = find_signif(p[valid_genes], **kwargs)
        is_signif[valid_genes] = is_signif_valid
        n_counts[is_signif] += 1
        print('\t%s: # implicated = %d'
               % (ds.cancer_type, np.count_nonzero(is_signif)))

        ds.results['t_test'] = (t, p, is_signif)

    print('\n\tCount Summary')
    for i in range(1, max(n_counts) + 1):
        n = np.count_nonzero(n_counts== i)
        print('\t\t%4d lncRNAs implicated in %2d cancer types.' % (n, i))


if __name__ == "__main__":
    print('\n1-Beginning Data Import')
    datasets = import_all_data()

    print('\n2-Performing t-tests')
    perform_t_test(datasets, procedure='bonferoni')







