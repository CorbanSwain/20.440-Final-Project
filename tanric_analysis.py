#!/bin/python3
# tanric_analysis.py
# Corban Swain, 2018


from utils import *
import datetime
import scipy.io as sio
from scipy.stats import ttest_ind
from paper_figures import *


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

        cachefile = os.path.join('data', 'np_cache', name + '.npy')
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


def assess_validity(datasets, expr_cutoff):
    for ds in datasets:
        ds.results['is_nonzero'] = np.logical_or(
            np.any(ds.normal_samples > 0, 1),
            np.any(ds.tumor_samples > 0, 1))
        ds.results['is_expressed'] = np.logical_or(
            np.mean(ds.normal_samples, 1) > expr_cutoff,
            np.mean(ds.tumor_samples, 1) > expr_cutoff)


def perform_t_test(datasets, test, t_filter, procedure, **kwargs):
    print('\nPerforming t-tests ...')
    # FIXME - May need to do some memory management here
    # FIXME - Validity testing should be done in its own function
    find_signif = procedure.signif_func()
    n_counts = np.zeros((TanricDataset.n_genes,), dtype=int)
    all_valid = np.zeros((TanricDataset.n_genes,), dtype=bool)
    signif_mat = np.zeros((TanricDataset.n_genes, len(datasets)), dtype=bool)
    for i, ds in enumerate(datasets):
        is_valid = ds.results[t_filter]
        all_valid = np.logical_or(all_valid, is_valid)
        norm_valid = ds.normal_samples[is_valid]
        tumor_valid = ds.tumor_samples[is_valid]

        t = np.zeros((ds.n_genes,))
        p = np.zeros((ds.n_genes,))

        if test is 'mwu':
            idxs = list(np.where(is_valid)[0])
            for j in range(np.count_nonzero(is_valid)):
                t[idxs[j]], p[idxs[j]] = mwu((tumor_valid[j, :],
                                              norm_valid[j, :]))
        elif test is 't_test':
            t[is_valid], p[is_valid] = ttest_ind(tumor_valid,
                                                 norm_valid,
                                                 axis=1)
        else:
            raise ValueError('Unexpected test specification -> \'%s\'.' % test)

        is_signif = np.zeros((ds.n_genes,), dtype=bool)
        is_signif_valid = find_signif(p[is_valid], **kwargs)
        is_signif[is_valid] = is_signif_valid
        signif_mat[:, i] = is_signif
        n_counts[is_signif] += 1
        stdout.write('\r\t%s: # implicated = %d'
                     % (ds.cancer_type, np.count_nonzero(is_signif)))
        stdout.flush()
        ds.results['t_test'] = (t, p, is_signif)
    stdout.write('\r\tDone.\n')

    signif_path = os.path.join('data', 'matlab_io', 'signif_matrix.m')
    sio.savemat(signif_path, {'isSignif': signif_mat,
                              'geneNames': matlab_cell_arr(
                                  TanricDataset.gene_info['code'])})

    print('\n\tCount Summary')
    print('\t\t%4d lncRNAs with significant expression in at least one '
          'cancer type.\n' % np.count_nonzero(all_valid))
    for i in range(1, max(n_counts) + 1):
        n = np.count_nonzero(n_counts == i)
        print('\t\t%4d lncRNAs implicated in %2d cancer types.' % (n, i))


def make_composite_dataset(datasets, filter_method, metric,
                           samples, fold_change_fudge):
    print('\tMaking Composite Dataset ...')

    # FIXME - There is a better way than all these if statements,
    # maybe implement a set of lambda functions and call them for all of the
    # operations.
    if metric in [Metric.FC_PAIR, Metric.MEAN2MEAN]:
        assert samples is Samples.TUMOR

    # Filtering and counting
    n_datasets = len(datasets)
    n_cancer_samples = 0
    n_normal_samples = 0
    n_pairs = 0
    all_valid = np.ones(TanricDataset.n_genes, dtype=bool)
    num_signif = np.zeros(TanricDataset.n_genes, dtype=int)
    signif_matrix = np.zeros((TanricDataset.n_genes, n_datasets))
    if filter_method is Filter.T_TEST:
        any_signif = np.zeros(TanricDataset.n_genes, dtype=bool)

    num_annotations_v1 = np.zeros(TanricDataset.n_genes, dtype=int)
    num_annotations_v2 = np.zeros(TanricDataset.n_genes, dtype=int)
    for i_gene in range(TanricDataset.n_genes):
        num_annotations_v1[i_gene] += len(TanricDataset.transcripts[i_gene])
        # FIXME - need to try except here
        try:
            num_annotations_v2[i_gene] += max(
                ts.n_aliases for ts in TanricDataset.transcripts[i_gene])
        except ValueError:
            num_annotations_v2[i_gene] = 0

    for i_ds, ds in enumerate(datasets):
        n_cancer_samples += ds.n_tumor_samples
        n_normal_samples += ds.n_normal_samples
        n_pairs += ds.n_pairs
        _, _, is_signif = ds.results['t_test']
        signif_matrix[:, i_ds] = is_signif
        num_signif += is_signif.astype(int)
        if filter_method is Filter.NONE:
            is_valid = ds.results['is_nonzero']
        elif filter_method in [Filter.THRESHOLD, Filter.T_TEST]:
            is_valid = ds.results['is_expressed']
        all_valid = np.logical_and(all_valid, is_valid)
        if filter_method is Filter.T_TEST:
            any_signif = np.logical_or(any_signif, is_signif)

    if filter_method is Filter.T_TEST:
        all_valid = np.logical_and(all_valid, any_signif)


    # setting up combined dataset
    genes_names = TanricDataset.gene_info['code'].copy(order='C')
    gene_descr = TanricDataset.gene_info['description'].copy(order='C')
    signif_matrix = signif_matrix[all_valid, :]
    genes_names = genes_names[all_valid]
    gene_descr = gene_descr[all_valid]
    num_signif = num_signif[all_valid]
    num_annotations_v1 = num_annotations_v1[all_valid]
    num_annotations_v2 = num_annotations_v2[all_valid]
    n_genes = np.count_nonzero(all_valid)
    print('\t%d Genes considered significant' % n_genes)

    if metric is Metric.FC_PAIR:
        n_samples = n_pairs
    elif metric is Metric.MEAN2MEAN:
        n_samples = n_datasets
    else:
        if samples is Samples.NORMAL:
            n_samples = n_normal_samples
        elif samples is Samples.TUMOR:
            n_samples = n_cancer_samples
        elif samples is Samples.ALL:
            n_samples = n_cancer_samples + n_normal_samples

    combined_values = np.zeros((n_genes, n_samples), dtype=np.float64)
    group_labels = np.zeros(n_samples, dtype='|S20')
    label_list = []
    label_counter = 0
    group_numbers = np.zeros(n_samples, dtype=int)
    insert_idx = 0

    for ds in datasets:
        if metric == Metric.FC_PAIR:
            end_idx = insert_idx + ds.n_pairs
        elif metric == Metric.MEAN2MEAN:
            end_idx = insert_idx + 1
        else:
            if samples is Samples.NORMAL:
                end_idx = insert_idx + ds.n_normal_samples
            elif samples is Samples.TUMOR:
                end_idx = insert_idx + ds.n_tumor_samples
            elif samples is Samples.ALL:
                end_idx = insert_idx + ds.n_samples

        insert_range = np.arange(insert_idx, end_idx, dtype=int)
        insert_idx = end_idx

        if metric is Metric.RPKM:
            value_matrix = ds.exprdata[all_valid]
        elif metric in [Metric.MEAN2MEAN, Metric.FC_MEAN]:
            norm_expr_val = np.mean(ds.normal_samples[all_valid], 1)
            norm_expr_val += fold_change_fudge
            value_matrix = ((ds.exprdata[all_valid].T +
                            fold_change_fudge) / norm_expr_val).T
            value_matrix = np.log2(value_matrix)
        elif metric is Metric.FC_PAIR:
            norm_sel, tumor_sel = ds.sample_pairs
            # FIXME - this could be implemented better, doing excess operations
            normal = ds.exprdata[:, norm_sel] + fold_change_fudge
            tumor = ds.exprdata[:, tumor_sel] + fold_change_fudge
            value_matrix = tumor[all_valid, :] / normal[all_valid, :]
            value_matrix = np.log2(value_matrix)

        if samples is Samples.NORMAL:
            combined_values[:, insert_range] = value_matrix[:, ds.normal_sel]
        elif samples is Samples.TUMOR:
            if metric is Metric.MEAN2MEAN:
                combined_values[:, insert_range] \
                    = np.mean(value_matrix[:, ds.tumor_sel], 1, keepdims=True)
            elif metric in [Metric.FC_MEAN, Metric.RPKM]:
                combined_values[:, insert_range] = value_matrix[:, ds.tumor_sel]
            elif metric is Metric.FC_PAIR:
                combined_values[:, insert_range] = value_matrix
        elif samples is Samples.ALL:
            combined_values[:, insert_range] = value_matrix

        cancer_name = ds.cancer_type
        norm_name = 'Normal-' + cancer_name
        for i, idx in enumerate(insert_range):
            if samples is Samples.NORMAL:
                group_labels[idx] = norm_name
                group_numbers[idx] = label_counter
            elif samples is Samples.TUMOR:
                group_labels[idx] = cancer_name
                group_numbers[idx] = label_counter
            elif samples is Samples.ALL:
                if i < ds.n_normal_samples:
                    group_labels[idx] = norm_name
                    group_numbers[idx] = label_counter
                else:
                    group_labels[idx] = cancer_name
                    group_numbers[idx] = label_counter + 1

        if samples is Samples.ALL:
            label_counter = label_counter + 2
            label_list.append(norm_name)
            label_list.append(cancer_name)
        else:
            label_counter = label_counter + 1
            if samples is Samples.NORMAL:
                label_list.append(norm_name)
            elif samples is Samples.TUMOR:
                label_list.append(cancer_name)

    print('\tDone.')
    return (combined_values,
            group_labels,
            group_numbers,
            strlist2np(label_list),
            genes_names,
            gene_descr,
            num_signif,
            signif_matrix,
            num_annotations_v1,
            num_annotations_v2,
            n_genes)


def save_for_matlab(datasets, comp_ds, settings):
    print('\nSaving results for MATLAB ...')
    spinner.start()
    n_sets = len(datasets)
    matlab_struct = {}
    for ds in datasets:
        t, p, fc, is_valid, is_signif = ds.results['t_test']
        logfc = np.zeros(fc.shape)
        logfc[is_valid] = np.log2(fc[is_valid])
        sample_names = matlab_cell_arr(ds.sample_names)
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
    gene_ids = matlab_cell_arr(TanricDataset.gene_ids)

    # Using ... .copy(order='C') to ensure the slices are contiguous; see
    # https://stackoverflow.com/questions/26778079/valueerror-ndarray-is-not
    # -c-contiguous-in-cython
    gene_codes = TanricDataset.gene_info['code'].copy(order='C')
    gene_codes = matlab_cell_arr(gene_codes)
    gene_descriptions = TanricDataset.gene_info['description'].copy(order='C')
    gene_descriptions = matlab_cell_arr(gene_descriptions)

    comp_set, labels, nums, comb_genes = comp_ds
    labels = matlab_cell_arr(labels)
    comb_genes = matlab_cell_arr(comb_genes)

    matpath = os.path.join('data', 'matlab_io', '%s_analysis_v%s.mat'
                           % ('%s', settings['version']))
    spinner.stop()
    print('\t' + matpath % 't_test')
    spinner.start()
    sio.savemat(matpath % 't_test',
                {'S': matlab_struct,
                 'nGenes': n_genes,
                 'geneIDs': col_vec(gene_ids),
                 'geneCodes': col_vec(gene_codes),
                 'geneDescriptions': col_vec(gene_descriptions),
                 'analysisMetadata': settings})
    spinner.stop()
    print('\t' + matpath % 'combined')
    spinner.start()
    sio.savemat(matpath % 'combined',
                {'analysisMetadata': settings,
                 'combGenes': col_vec(comb_genes),
                 'combData': comp_set,
                 'combLabels': labels,
                 'combNumIDs': nums})
    spinner.stop()
    print('\tDone.')


def save_for_matlab_2(filename, combined_data, settings):
    print('\tSaving results for MATLAB ...')
    spinner.start()

    ml_settings = {}
    for k, v in settings.items():
        if isinstance(v, Enum):
            ml_settings[k] = v.value
        else:
            ml_settings[k] = v

    values, g_labels, g_numbers, label_list, gene_names, gene_descr, \
        num_signif, signif_mat, n_anat1, n_anat2, n_genes = combined_data
    n_samples = len(g_labels)
    g_labels = matlab_cell_arr(g_labels)
    label_list = matlab_cell_arr(label_list)
    gene_names = matlab_cell_arr(gene_names)
    gene_descr = matlab_cell_arr(gene_descr)

    matpath = os.path.join('data', 'matlab_io',
                           'multi_analysis', filename)

    sio.savemat(matpath,
                {'analysisMetadata': ml_settings,
                 'values': values,
                 'numSamples': n_samples,
                 'sampleNames': row_vec(g_labels),
                 'sampleGroupNumbers': row_vec(g_numbers),
                 'sampleList': label_list,
                 'numGenes': n_genes,
                 'geneNames': col_vec(gene_names),
                 'geneDescrips': col_vec(gene_descr),
                 'geneNumSignif': col_vec(num_signif),
                 'signifMatrix': signif_mat,
                 'geneNumAnat1': col_vec(n_anat1),
                 'geneNumAnat2': col_vec(n_anat2)})
    spinner.stop()
    print('\tDone.')


def make_multi_analysis(datasets, settings):
    file_list = []

    for filt in [Filter.NONE]:
        for metric in [Metric.FC_MEAN, Metric.FC_PAIR]:
            for samples in [Samples.TUMOR]:
                print('\nANALYSIS - %s - %s - %s' % (filt, metric, samples))

                settings['filter_method'] = filt
                settings['metric'] = metric
                settings['samples'] = samples

                try:
                    data = make_composite_dataset(
                        datasets,
                        filt,
                        metric,
                        samples,
                        settings['fold_change_fudge'])
                except AssertionError:
                    print('\tEXCEPTION: Couldn\'t Process')
                    continue

                if settings['do_save']:
                    filename = '%s-%s-%s_v%s.mat' % (filt,
                                                     metric,
                                                     samples,
                                                     settings['version'])

                    save_for_matlab_2(filename, data, settings)
                    file_list.append(filename)

    if settings['do_save']:
        path = os.path.join('data', 'matlab_io',
                            'multi_analysis', 'file_list.mat')
        sio.savemat(path, {'fileNames': matlab_cell_arr(file_list)})


if __name__ == "__main__":
    settings = {
        'min_norm_samples': 5,
        'test': 'mwu',  # mwu, t_test
        'version': '6.2',  # TODO - some type of automatic versioning
        'expression_cutoff': 0.3,  # 0.3 used in TANRIC paper
        'filter_method': None,
        't_filter': 'is_expressed', # is_nonzero, is_expressed
        'multi_hyp_procedure': MultiHypProc.BONFERONI,
        'alpha_crit': 5e-3,
        'metric': None,
        'samples': None,
        'fold_change_fudge': 1e-4,
        'do_save': False,
        'do_plot': False,
        'analysis_date': str(datetime.datetime.now())
    }

    if settings['multi_hyp_procedure'] is MultiHypProc.BEN_HOCH:
        add_args = {'q': settings['alpha_crit'],
                    'plot': False}
    else:
        add_args = {'a': settings['alpha_crit']}

    datasets = import_all_data(settings['min_norm_samples'])

    TanricDataset.get_transcript_info()

    TanricDataset.get_gene_info()

    assess_validity(datasets, settings['expression_cutoff'])

    perform_t_test(datasets,
                   settings['test'],
                   settings['t_filter'],
                   settings['multi_hyp_procedure'],
                   **add_args)

    make_multi_analysis(datasets, settings)

    plt.style.use('seaborn-notebook')

    # make_simple_charts(datasets)

    make_ma_plots(datasets,
                 settings['fold_change_fudge'])



