from seqdataloader import *
classification_params={
    'task_list':"tasks.tsv",
    'outf':"classificationlabels.SummitWithin200bpCenter.tsv.gz",
    'output_type':'gzip',
    'chrom_sizes':'hg38.chrom.sizes',
    'chroms_to_keep':['chr21'],
    "store_positives_only":True,
    'bin_stride':50,
    'left_flank':400,
    'right_flank':400,
    'bin_size':200,
    'threads':10,
    'subthreads':4,
    'allow_ambiguous':True,
    'labeling_approach':'peak_summit_in_bin_classification'
    }
genomewide_labels(classification_params)

regression_params={
    'task_list':"tasks.tsv",
    'outf':"regressionlabels.all_genome_bins_regression.hdf5",
    'output_type':'hdf5',
    'chrom_sizes':'hg38.chrom.sizes',
    'store_values_above_thresh': 0,
    'chroms_to_keep':['chr21'],
    'bin_stride':50,
    'left_flank':400,
    'right_flank':400,
    'bin_size':200,
    'threads':10,
    'subthreads':4,
    'labeling_approach':'all_genome_bins_regression'
    }
genomewide_labels(regression_params)
