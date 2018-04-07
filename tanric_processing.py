#!/bin/python3
# tanric_processing.py
# Corban Swain, 2018

import os
import matplotlib.pyplot as plt
import numpy as np


def file2dict(file, delimiter='\t'):
    data = np.genfromtxt(file, dtype=None,
                         delimiter=delimiter,
                         encoding=None)
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


class Dataset:
    def __init__(self, metadict, exprdata=None):
        self.metadict = metadict
        self.cancer_type = metadict['Cancer_Type']
        self.n_normal_samples = metadict['Num_Normal_Samples']
        self.normal_sel = np.arange(self.n_normal_samples)
        self.n_tumor_samples = metadict['Num_Tumor_Samples']
        self.n_samples = self.n_normal_samples + self.n_tumor_samples
        self.tumor_sel = np.arange(self.n_normal_samples, self.n_samples)
        if exprdata is not None:
            self.parse_exprdata(exprdata)
        else:
            self.gene_ids = None
            self.n_genes = None
            self.sample_names = None
            self.values = None

    def parse_exprdata(self, exprdata):
        self.gene_ids = [s.strip('\"') for s in exprdata['Gene_ID']]
        self.n_genes = len(self.gene_ids)
        self.sample_names = list(exprdata.dtype.names[1:])
        self.values = exprdata[self.sample_names]
        self.values = self.values.view(np.float64).reshape((self.n_genes, -1))

        if (self.n_samples - self.values.shape[1]) < 1:
            raise ValueError('Metadata value for n_samples (%d) is not '
                             'equal to the number of columns found in the '
                             'expression data file (%d).'
                             % (self.n_samples, self.values.shape[1]))

    @property
    def normal_samples(self):
        return self.values[:, self.normal_sel]

    @property
    def tumor_samples(self):
        return self.values[:, self.tumor_sel]



m_dict = file2dict('data/tanric_data/TCGA-BLCA-rnaexpr-META.tsv')
data = np.genfromtxt('data/tanric_data/TCGA-BLCA-rnaexpr.tsv', delimiter='\t',
                  dtype=None, encoding=None, names=True, deletechars='')
dataset = Dataset(m_dict, data)