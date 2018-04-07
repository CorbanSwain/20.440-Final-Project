#!/bin/python3
# tanric_utils.py
# Corban Swain, 2018

import numpy as np
import warnings

def structarr2nparr(x):
    # TODO - check for max precision
    with warnings.catch_warnings():
        # Ignoring warning about viewing an array returned by selecting
        # multiple fields in a structures array.
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


class TanricDataset:
    def __init__(self, metadict, expr_structarr=None):
        self.metadict = metadict
        self.cancer_type = metadict['Cancer_Type']
        self.n_normal_samples = metadict['Num_Normal_Samples']
        self.normal_sel = np.arange(self.n_normal_samples)
        self.n_tumor_samples = metadict['Num_Tumor_Samples']
        self.n_samples = self.n_normal_samples + self.n_tumor_samples
        self.tumor_sel = np.arange(self.n_normal_samples, self.n_samples)
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