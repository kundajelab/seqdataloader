from seqdataloader import *
params={
    'task_list':"/users/annashch/seqdataloader/examples/tasks.tsv",
    'outf':"classificationlabels.SummitWithin200bpCenter.tsv.gz",
    'output_type':'gzip',
    'chrom_sizes':'hg38.chrom.sizes',
    'chroms_to_keep':['chr21'],
    'bin_stride':50,
    'left_flank':400,
    'right_flank':400,
    'bin_size':200,
    'threads':10,
    'subthreads':4,
    'allow_ambiguous':True,
    'labeling_approach':'peak_summit_in_bin_classification'
    }

genomewide_labels(params)
