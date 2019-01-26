# seqdataloader
Sequence data label generation and ingestion into deep learning models

## Installation
`pip install seqdataloader`

If you plan to modify the code, you can install it in development mode: 
`pip install -e seqdataloader` 

## Quick Start

tasks.tsv is a 3 column tab-delimited file

Column 1 -- User-specified task name

Column 2 -- path to narrowPeak file (or empty)

Column 3 -- peath to bigwig file (or empty)

```
genomewide_labels --task_list tasks.tsv \
		  --outf classificationlabels.SummitWithin200bpCenter.tsv.gz \
		  --output_type gzip \ # (one of gzip, bz2, hdf5, pkl) 
		  --chrom_sizes hg38.chrom.sizes \
		  --bin_stride 50 \
		  --left_flank 400 \
		  --right_flank 400 \
		  --bin_size 200 \
		  --threads 10 \
		  --subthreads 4 \
		  --allow_ambiguous \
		  --labeling_approach peak_summit_in_bin_classification 
```
labeling_approach can be one of:

    "peak_summit_in_bin_classification"

    "peak_percent_overlap_with_bin_classification"

    "peak_summit_in_bin_regression"

    "peak_percent_overlap_with_bin_regression"
    
    "all_genome_bins_regression"
    

## How to run 
Sample datasets are included in the folder `examples/peak_files_from_encode_for_label_comparison` and `examples/bigwig_files_from_encode_for_label_comparison`

Execute the script:

`examples/genomewide_labels.sh` for examples on how to generate classification and regression labels on sample datasets.
The script generates binary classification labels (1,0,-1 for ambiguous) or continuous regression labels reflective of bigWig coverage in a bin  in bed file format:

http://mitra.stanford.edu/kundaje/seqdataloader/classificationlabels.50PercentOverlap.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/classificationlabels.SummitWithin200bpCenter.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/regressionlabels.50PercentOverlap.tsv.gz

http://mitra.stanford.edu/kundaje/seqdataloader/regressionlabels.SummitWithin200bpCenter.tsv.gz

Corresponding WashU Browser Tracks with optimal narrowPeak and associated bin labels are here:
http://epigenomegateway.wustl.edu/legacy/?genome=hg38&session=GDB2BTMGnB&statusId=1154897038


## Dependencies

Please make sure the following dependencies are installed on your system to use SeqDataLoader:
* pybedtools
* pyBigWig 
* pandas
* numpy
* multiprocessing


## A note on file outputs

The code supports several output types: `hdf5`, `gzip`, `pkl`, `bz2`.
Specify your desired output type with the flag `--output_type`. The default setting for this flag is `gzip`
Please note that the large bottleneck in the code is writing the files to disk. `hdf5` has negligible overhead, but using `gzip` or `bz2` may increase runtime. Timining benchmarks are provided in `examples/genomewide_labels.sh`

You may speed up i/o by writing chromosome outputs to separate files in parallel. This is currently only supported for the `gzip` and `bz2` output types, as i/o is less of a bottleneck for `hdf5` and `pkl` output formats. Use the flag `--split_output_by_chrom` to invoke this parallelized saving of chromosomes.

## Documentation and benchmarks

Testing, benchmarks, and documentation can be found in the `docs` folder
