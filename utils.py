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
import os
import time
import threading
from enum import Enum
import requests
import pickle



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


def strlist2np(x):
    max_len = max(len(s) for s in x)
    dtype_str = '|S%d' % max_len
    np_arr = np.zeros(len(x), dtype=dtype_str)
    for i, s in enumerate(x):
        np_arr[i] = s
    return np_arr


def matlab_cell_arr(x):
    return np.array(x, dtype=np.object)


def col_vec(x):
    return x.reshape((-1, 1))


def row_vec (x):
    return x.reshape((1, -1))


def benhoch(p, q, plot=False):
    """benjamani-hochberg procedure for controlling the false discovery rate."""
    i = rankdata(p, method='ordinal')
    m = len(p)
    bh_crit = q * i / m
    try:
        cutoff_i = np.max(i[p <= bh_crit])
    except ValueError:
        return np.array([])
    is_signif = i <= cutoff_i

    if plot:
        s = np.argsort(i)
        plt.plot(i[s], bh_crit[s], i[s], p[s])
        plt.plot(i[is_signif], p[is_signif], 'r.')
        plt.xlabel('Rank')
        plt.ylabel('P Value')
        plt.xlim([0, m])
        plt.ylim([0, 1])
        plt.show()

    return is_signif


class StrEnum(str, Enum):

    def __str__(self):
        return self.value


class Filter(StrEnum):
    T_TEST = 't_test'
    THRESHOLD = 'threshold'
    NONE = 'no_filter'


class MultiHypProc(StrEnum):
    BONFERONI = 'bonferoni'
    CRITICAL_VALUE = 'crit'
    BEN_HOCH = 'ben_hoch'

    def signif_func(self):
        return {
            'ben_hoch': benhoch,
            'crit': lambda p, a: p <= a,
            'bonferoni': lambda p, a: p <= (a / len(p))
        }[self.value]


class Metric(StrEnum):
    FC_MEAN = 'fold_change_mean'
    FC_PAIR = 'fold_change_pairwise'
    RPKM = 'rpkm'
    MEAN2MEAN = 'mean_to_mean'


class Samples(StrEnum):
    TUMOR = 'tumor'
    NORMAL = 'normal'
    ALL = 'all'


ec = NCBI_Client()


def clean_id(gid):
    return gid.split('.')[0]

def geneid2name(gid):
    """Queries the NCBI database for gene HGNC codes and descriptions.

    :param gids: a string contaning a gene id
    :return:
    """
    cid = clean_id(gid)
    search_result = ec.esearch(db='gene', term=cid)
    try:
        if not search_result.ids:
            raise IndexError
        egene = ec.efetch(db='gene',
                          id=search_result.ids[0]).entrezgenes[0]
        name_tuple = (egene.hgnc, egene.description)
    except (IndexError,
            eutils.exceptions.EutilsNCBIError,
            eutils.exceptions.EutilsError,
            eutils.exceptions.EutilsRequestError,
            eutils.exceptions.EutilsNotFoundError):
        name_tuple = (cid, gid)
    return name_tuple


class Spinner:
    """
    Spinning cursor for console.
    """
    busy = False
    delay = 0.2

    @staticmethod
    def spinning_cursor():
        while True:
            for cursor in '|/-\\':
                yield cursor

    def __init__(self, delay=None):
        self.spinner_generator = self.spinning_cursor()
        if delay and float(delay): self.delay = delay

    def spinner_task(self):
        while self.busy:
            stdout.write('\r\t' + next(self.spinner_generator))
            time.sleep(self.delay)
            stdout.flush()
        stdout.write('\r')

    def start(self):
        self.busy = True
        threading.Thread(target=self.spinner_task).start()

    def stop(self):
        self.busy = False
        time.sleep(self.delay)


spinner = Spinner()


class TanricDataset:
    # DONE - gene IDs and lists should be class - level atributtes
    # FIXME - there should be a gene class to handle all of this
    gene_ids = None
    n_genes = None
    gene_info = None
    transcripts = []

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

        self.gene_ids = None
        self.n_genes = None
        self.sample_names = None
        self.exprdata = None
        self.sample_pairs = None
        self.n_pairs = None

        self.results = {}

        if expr_structarr is not None:
            self.parse_exprdata(expr_structarr)

    def parse_exprdata(self, expr_structarr):
        self.gene_ids = [s.strip('\"') for s in expr_structarr['Gene_ID']]
        self.n_genes = len(self.gene_ids)

        if TanricDataset.n_genes is None:
            TanricDataset.n_genes = self.n_genes
        else:
            assert self.n_genes == TanricDataset.n_genes
        if TanricDataset.gene_ids is None:
            TanricDataset.gene_ids = self.gene_ids
        else:
            assert np.all(self.gene_ids == TanricDataset.gene_ids)

        self.sample_names = list(expr_structarr.dtype.names[1:])
        self.exprdata = structarr2nparr(expr_structarr[self.sample_names])
        if not (self.n_samples - self.exprdata.shape[1]) < 1:
            raise ValueError('Metadata value for n_samples (%d) is not '
                             'equal to the number of columns found in the '
                             'expression data file (%d).'
                             % (self.n_samples, self.exprdata.shape[1]))

        self.pair_samples()

    @staticmethod
    def get_patient_id(name):
        return '-'.join(name.split('-')[-2:])

    def pair_samples(self):
        normal_ids = [self.get_patient_id(n) for n in self.normal_names]
        tumor_ids = [self.get_patient_id(n) for n in self.tumor_names]

        pairs = []
        for i_normal, n_id in enumerate(normal_ids):
            try:
                i_tumor = tumor_ids.index(n_id) + self.n_normal_samples
            except ValueError:
                continue
            pairs.append((i_normal, i_tumor))

        if pairs:
            self.n_pairs = len(pairs)
            norm_pair_idx, tumor_pair_idx = zip(*pairs)
            self.sample_pairs = (np.array(norm_pair_idx, dtype=int),
                                 np.array(tumor_pair_idx, dtype=int))
        else:
            self.n_pairs = 0

    @property
    def normal_names(self):
        return [self.sample_names[i] for i in range(self.n_samples)
                if self.normal_sel[i]]

    @property
    def tumor_names(self):
        return [self.sample_names[i] for i in range(self.n_samples)
                if self.tumor_sel[i]]

    @property
    def normal_samples(self):
        return self.exprdata[:, self.normal_sel]

    @property
    def tumor_samples(self):
        return self.exprdata[:, self.tumor_sel]

    @classmethod
    def get_gene_info(cls):
        print('\nFetching Gene Info ...')
        if cls.gene_info is None:
            gnamepath = os.path.join('data', 'np_cache', 'gene_names.npy')
            try:
                cls.gene_info = np.load(gnamepath)
            except FileNotFoundError:
                name_arr = np.zeros(cls.n_genes,
                                    dtype={'names': ['code', 'description'],
                                           'formats': ['|S20', '|S200']})
                for i, gid in enumerate(cls.gene_ids):
                    name_tuple = geneid2name(gid)
                    name_arr[i] = name_tuple
                    stdout.write('\r\t%05d - %05.1f %% - \"%s\"'
                                 % (i, 100 * i / cls.n_genes, name_tuple[1]))
                    stdout.flush()
                stdout.write('\n')
                cls.gene_info = name_arr
                np.save(gnamepath, name_arr)
        print('\tDone.')
        return cls.gene_info

    @classmethod
    def get_transcript_info(cls):
        print('\nLoading all Transcript Info ...')
        cache_path = os.path.join('data', 'np_cache', 'lncpedia_info.p')
        try:
            cls.transcripts = pickle.load(open(cache_path, 'rb'))
        except FileNotFoundError:
            for i, gid in enumerate(cls.gene_ids):
                response_dict = query_lncpedia(gid)
                if response_dict['count'] < 1:
                    cls.transcripts.append([])
                    continue
                ref_genome = response_dict['refGenome']
                ts_list = []
                for ts in response_dict['transcripts']:
                    ts_list.append(Transcript(ref_genome, ts))
                cls.transcripts.append(ts_list)
                stdout.write('\r\t%05d - %05.1f%% - %s'
                             % (i,
                                i / cls.n_genes * 100,
                                cls.transcripts[-1][0].lncpedia_id))
                stdout.flush()
            stdout.write('\r\tDone.\n')
            pickle.dump(cls.transcripts, open(cache_path, 'wb'))



lncpedia_url = 'https://lncipedia.org/api/search'


def query_lncpedia(gid):
    resp = requests.get(lncpedia_url, params={'id': gid})
    return resp.json()


class Transcript:

    def __init__(self, ref_genome, transcript_dict):
        self.ref_genome = ref_genome
        self.sequence = transcript_dict['sequence']
        self.chromosome = transcript_dict['chromosome']
        self.start = transcript_dict['start']
        self.end = transcript_dict['end']
        self.strand = transcript_dict['strand']
        self.size = transcript_dict['transcriptSize']
        self.lncpedia_id = transcript_dict['lncipediaTranscriptID']
        self.n_aliases = len(transcript_dict['transcriptAliases'])