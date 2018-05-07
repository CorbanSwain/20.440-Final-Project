#!/bin/python3
# utils.py
# Corban Swain, 2018

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import rankdata, mannwhitneyu
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
import multiprocessing


# parmap implementation for functions with non-trivial execution time
def funmap(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))
    return


def parmap(f, X, nprocs=(multiprocessing.cpu_count())):
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()
    proc = [multiprocessing.Process(
        target=funmap,
        args=(f, q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(X))]
    [p.join() for p in proc]
    return [x for i, x in sorted(res)]


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


def mwu(tpl):
    a, b = tpl
    return mannwhitneyu(a, b, alternative='two-sided')

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

class Results(StrEnum):
    SIGNIF_TEST = 'significance_test'
    VALIDITY = 'validity'


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
    # FIXME - there should be a gene class to handle all of this
    fcf = None
    gene_ids = None
    n_genes = None
    gene_info = None
    transcripts = []
    cancer_dict = {}

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

        self._normal_mean = None
        self._tumor_mean = None
        self._normal_pair_samples = None
        self._tumor_pair_samples = None
        self._paired_fc = None
        self._paired_fc_mean = None
        self._paired_fc_median = None
        self._fc = None
        self._fc_mean = None
        self._normal_names = None
        self._tumor_names = None
        self._normal_samples = None
        self._tumor_samples = None

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
    def normal_mean(self):
        if self._normal_mean is None:
            self._normal_mean = np.mean(self.normal_samples, 1)
        return self._normal_mean

    @property
    def tumor_mean(self):
        if self._tumor_mean is None:
            self._tumor_mean = np.mean(self.tumor_samples, 1)
        return self._tumor_mean

    @property
    def normal_pair_samples(self):
        if self._normal_pair_samples is None:
            idxs, _ = self.sample_pairs
            self._normal_pair_samples = self.exprdata[:, idxs]
        return self._normal_pair_samples

    @property
    def tumor_pair_samples(self):
        if self._tumor_pair_samples is None:
            _, idxs = self.sample_pairs
            self._tumor_pair_samples = self.exprdata[:, idxs]
        return self._tumor_pair_samples

    @property
    def paired_fc(self):
        if self._paired_fc is None:
            self._paired_fc = ((self.tumor_pair_samples + self.fcf) /
                               (self.normal_pair_samples + self.fcf))
            self._paired_fc = np.log2(self._paired_fc)
        return self._paired_fc

    @property
    def paired_fc_mean(self):
        if self._paired_fc_mean is None:
            self._paired_fc_mean = np.mean(self.paired_fc, 1)
        return self._paired_fc_mean


    @property
    def paired_fc_median(self):
        if self._paired_fc_median is None:
            self._paired_fc_median = np.median(self.paired_fc, 1)
        return self._paired_fc_median

    @property
    def fc(self):
        if self._fc is None:
            self._fc = ((self.tumor_samples.T + self.fcf) /
                        self.normal_mean + self.fcf).T
            self._fc = np.log2(self._fc)
        return self._fc

    @property
    def fc_mean(self):
        if self._fc_mean is None:
            self._fc_mean = ((self.tumor_mean + self.fcf) /
                             (self.normal_mean + self.fcf))
            self._fc_mean = np.log2(self._fc_mean)
        return self._fc_mean

    @property
    def normal_names(self):
        if self._normal_names is None:
            self._normal_names = [self.sample_names[i]
                                  for i in range(self.n_samples)
                                  if self.normal_sel[i]]
        return self._normal_names

    @property
    def tumor_names(self):
        if self._tumor_names is None:
            self._tumor_names = [self.sample_names[i]
                                 for i in range(self.n_samples)
                                 if self.tumor_sel[i]]
        return self._tumor_names

    @property
    def normal_samples(self):
        if self._normal_samples is None:
            self._normal_samples = self.exprdata[:, self.normal_sel]
        return self._normal_samples

    @property
    def tumor_samples(self):
        if self._tumor_samples is None:
            self._tumor_samples = self.exprdata[:, self.tumor_sel]
        return self._tumor_samples

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

        if cls.transcripts:
            for i, ts_list in enumerate(cls.transcripts):
                if ts_list:
                    new_code = ts_list[0].lncpedia_name
                    if 'ENSG' in str(cls.gene_info['code'][i]):
                        cls.gene_info['code'][i] = new_code

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
        return cls.transcripts

    @classmethod
    def sampleid2name(cls, cid):
        if not cls.cancer_dict:
            fname = os.path.join('data', 'tcga_codes', 'diseaseStudy.tsv')
            cls.cancer_dict = file2dict(fname)
        return cls.cancer_dict[cid]


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

    @property
    def lncpedia_name(self):
        return self.lncpedia_id.split(':')[0]


class CompositeDataset:
    def __init__(self, datasets):
        self.dss = datasets

        self._n_dss = None
        self._n_samples = None
        self._dss_names = None
        self._n_genes = None
        self._long_panda = None
        self._gene_names = None
        self._num_signif = None
        self._all_expressed = None
        self._any_expressed = None
        self._all_nonzero = None
        self._any_nonzero = None

    @property
    def n_dss(self):
        if self._n_dss is None:
            self._n_dss = len(self.dss)
        return self._n_dss

    @property
    def n_samples(self):
        if self._n_samples is None:
            self._n_samples = sum(ds.n_pairs for ds in self.dss)
        return self._n_samples

    @property
    def dss_names(self):
        if self._dss_names is None:
            self._dss_names = [ds.cancer_type for ds in self.dss]
        return self._dss_names

    @property
    def n_genes(self):
        if self._n_genes is None:
            self._n_genes = TanricDataset.n_genes
        return self._n_genes

    @property
    def gene_names(self):
        if self._gene_names is None:
            self._gene_names = [bytes.decode(v, 'utf8')
                                for v in TanricDataset.gene_info['code']]
        return self._gene_names

    @property
    def num_signif(self):
        if self._num_signif is None:
            self._num_signif = np.sum(np.column_stack(
                [ds.results['t_test'][2] for ds in self.dss]), 1)
        return self._num_signif

    @property
    def all_expressed(self):
        if self._all_expressed is None:
            self._all_expressed = np.all(np.column_stack(
                [ds.results['is_expressed'] for ds in self.dss]), 1)
        return self._all_expressed

    @property
    def any_expressed(self):
        if self._any_expressed is None:
            self._any_expressed = np.any(np.column_stack(
                [ds.results['is_expressed'] for ds in self.dss]), 1)
        return self._any_expressed

    @property
    def all_nonzero(self):
        if self._all_nonzero is None:
            self._all_nonzero = np.all(np.column_stack(
                [ds.results['is_nonzero'] for ds in self.dss]), 1)
        return self._all_nonzero

    @property
    def any_nonzero(self):
        if self._any_nonzero is None:
            self._any_nonzero = np.any(np.column_stack(
                [ds.results['is_nonzero'] for ds in self.dss]), 1)
        return self._any_nonzero

    @property
    def long_panda(self):
        if self._long_panda is None:
            n_points = self.n_dss * self.n_genes
            gene_name = sum([self.gene_names for _ in self.dss], [])
            cancer_type = [t for t in self.dss_names
                           for _ in range(self.n_genes)]
            cancer_num_id = np.repeat(range(self.n_dss), self.n_genes)
            is_nonzero = np.concatenate([ds.results['is_nonzero'] for ds in
                                         self.dss])
            is_expressed = np.concatenate([ds.results['is_expressed'] for ds in
                                          self.dss])
            t_stat = np.concatenate([ds.results['t_test'][0] for ds in
                                     self.dss])
            q_values = np.concatenate([ds.results['q_values'] for ds in
                                       self.dss])
            is_signif = np.concatenate([ds.results['t_test'][2] for ds in
                                       self.dss])

            num_signif = np.zeros(n_points, dtype=int)
            for ds in self.dss:
                this_signif = np.tile(ds.results['t_test'][2], self.n_dss)
                num_signif[np.logical_and(is_signif, this_signif)] += 1

            expression = np.concatenate([ds.tumor_mean for ds in self.dss])

            fc = np.concatenate([ds.paired_fc_median for ds in self.dss])

            df_dict = {
                'gene_name': gene_name,
                'cancer_type': cancer_type,
                'cancer_num_id': cancer_num_id,
                'is_nonzero': is_nonzero,
                'is_expressed': is_expressed,
                't_stat': t_stat,
                'q_values': q_values,
                'is_signif': is_signif,
                'num_signif': num_signif,
                'expression': expression,
                'fc': fc
            }
            self._long_panda = pd.DataFrame.from_dict(df_dict)
        return self._long_panda

    def profile_panda(self, filt=None):
        if filt is None:
            filt = self.all_expressed

        mi_arrs = [[], []]
        arr = np.zeros((np.count_nonzero(filt), self.n_samples))
        col_idx = 0
        for ds in self.dss:
            n_samps = ds.paired_fc.shape[1]
            for iSamp in range(n_samps):
                mi_arrs[0].append(ds.cancer_type)
                mi_arrs[1].append(ds.sample_names[iSamp])
                arr[:, col_idx] = ds.paired_fc[filt, iSamp]
                col_idx += 1

        index = pd.MultiIndex.from_arrays(mi_arrs, names=['cancer_type',
                                                          'sample_name'])
        df = pd.DataFrame(arr,
                          index=[self.gene_names[i] for i in
                                 range(self.n_genes) if filt[i]],
                          columns=index)
        return df

    def median_profile_panda(self, filt=None):
        if filt is None:
            filt = self.all_expressed
        df = pd.DataFrame([], index=[self.gene_names[i] for i in
                                     range(self.n_genes) if filt[i]])
        for ds in self.dss:
            df[ds.cancer_type] = ds.paired_fc_median[filt]
        return df

