#!/bin/python3
# tanric_processing.py
# Corban Swain, 2018

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


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


d = file2dict('data/tanric_data/TCGA-BLCA-rnaexpr-META.tsv')


class Dataset:
    def __init__(self, metadict, exprdata):
        pass

