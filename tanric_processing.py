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
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass
        metadict[key] = value
    return metadict


m_dict = file2dict('data/tanric_data/TCGA-BLCA-rnaexpr-META.tsv')

data = np.genfromtxt('data/tanric_data/TCGA-BLCA-rnaexpr.tsv', delimiter='\t',
                  dtype=None, encoding=None, names=True, deletechars='')


class Dataset:
    def __init__(self, metadict, exprdata=None):
        self.metadict = metadict
        self.n_normal_samples = metadict['Num_Normal_Samples']
        self.n_tumor_samples = metadict['Num_Tumor_Samples']
        self.cancer_type = metadict['Cancer_Type']
        self.gene_ids = None
        self.n_genes = None
        self.sample_names = None
        self.values = None

        if exprdata is not None:
            self.parse_exprdata(exprdata)

    def parse_exprdata(self, exprdata):
        self.gene_ids = [s.strip('\"') for s in exprdata['Gene_ID']]
        self.n_genes = len(self.gene_ids)
        self.sample_names = list(exprdata.dtype.names[1:])
        self.values = exprdata[self.sample_names]
        self.values = self.values.view(np.float64).reshape((self.n_genes, -1))


dataset = Dataset(m_dict, data)