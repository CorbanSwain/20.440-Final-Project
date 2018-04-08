#!/bin/python3
# utils.py
# Corban Swain, 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import eutils.exceptions
from eutils.client import Client as NCBI_Client
from sys import stdout
import warnings


def structarr2nparr(x):
    # TODO - check for max precision
    with warnings.catch_warnings():
        # Ignoring warning about viewing an array returned by selecting
        # multiple fields in a structured array.
        # see: https://github.com/numpy/numpy/issues/8383
        warnings.simplefilter('ignore')
        return x.view(np.float64).reshape((len(x), -1))


def file2dict(file, delimiter='\t'):
    data = np.genfromtxt(file, dtype=None, delimiter=delimiter, encoding=None)
    metadict = {}
    for row in data:
        key = row[0]
        value = row[1]
        try:
            value = int(value)
        except ValueError:  # can't cast to int, try float
            try:
                value = float(value)
            except ValueError:  # can't cast to float, keep as string
                pass
        metadict[key] = value
    return metadict


# benjamani-hochberg procedure for
def benhoch(p, q=0.1, plot=False):
    i = rankdata(p, method='ordinal')
    m = len(p)
    bh_crit = q * i / m
    try:
        cutoff_i = np.max(i[np.where(p <= bh_crit)[0]])
    except ValueError:
        return np.array([])
    is_signif = i <= cutoff_i

    if plot:
        s = np.argsort(i)
        plt.plot(i[s], bh_crit[s],  i[s], p[s])
        plt.plot(i[signif_idx], p[signif_idx], 'r.')
        plt.xlabel('Rank')
        plt.ylabel('P Value')
        plt.xlim([0, m])
        plt.ylim([0, 1])
        plt.show()

    return is_signif


mh_tests = {'ben-hoch': benhoch,
            'crit': lambda p, a=0.05: p < a,
            'bonferoni': lambda p, a=0.05: p < (a / len(p))}


class TanricDataset:

    # FIXME - gene IDs and lists should be global

    def __init__(self, metadict, expr_structarr=None):
        self.metadict = metadict

        self.cancer_type = metadict['Cancer_Type']

        self.n_normal_samples = metadict['Num_Normal_Samples']
        self.n_tumor_samples = metadict['Num_Tumor_Samples']
        self.n_samples = self.n_normal_samples + self.n_tumor_samples

        self.normal_sel = np.zeros((self.n_samples,), dtype=bool)
        self.normal_sel[np.arange(self.n_normal_samples)] = True

        self.tumor_sel = np.zeros((self.n_samples,), dtype=bool)
        self.tumor_sel[np.arange(self.n_normal_samples, self.n_samples)] = True

        self.results = {}
        if expr_structarr is not None:
            self.parse_exprdata(expr_structarr)
        else:
            self.gene_ids = None
            self.n_genes = None
            self.sample_names = None
            self.exprdata = None

    def parse_exprdata(self, expr_structarr):
        self.gene_ids = [s.strip('\"') for s in expr_structarr['Gene_ID']]
        self.n_genes = len(self.gene_ids)
        self.sample_names = list(expr_structarr.dtype.names[1:])
        self.exprdata = structarr2nparr(expr_structarr[self.sample_names])
        if not (self.n_samples - self.exprdata.shape[1]) < 1:
            raise ValueError('Metadata value for n_samples (%d) is not '
                             'equal to the number of columns found in the '
                             'expression data file (%d).'
                             % (self.n_samples, self.exprdata.shape[1]))


    @property
    def normal_samples(self):
        return self.exprdata[:, self.normal_sel]

    @property
    def tumor_samples(self):
        return self.exprdata[:, self.tumor_sel]

    @staticmethod
    def geneid2name(ids):
        ec = NCBI_Client()
        n_ids = len(ids)
        name_arr = np.zeros(n_ids, dtype={'names': ['code', 'description'],
                                          'formats': ['|S20', '|S100']})
        for i, gid in enumerate(ids):
            id_clean = gid.split('.')[0]
            search_result = ec.esearch(db='gene', term=id_clean)
            try:
                if not search_result.ids:
                    raise IndexError
                gene = ec.efetch(db='gene',
                                 id=search_result.ids[0]).entrezgenes[0]
                name_tuple = (gene.hgnc, gene.description)
            except (IndexError,
                    eutils.exceptions.EutilsNCBIError,
                    eutils.exceptions.EutilsError,
                    eutils.exceptions.EutilsRequestError,
                    eutils.exceptions.EutilsNotFoundError):
                name_tuple = (gid, gid)
            name_arr[i] = name_tuple
            stdout.write('\r\t%05d - %05.1f %% - \"%s\"'
                         % (i, 100*i/n_ids, name_tuple[1]))
            stdout.flush()
        stdout.write('\n')
        return name_arr
