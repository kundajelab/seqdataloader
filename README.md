# seqdataloader
Sequence data label generation and ingestion into deep learning models

## Installation
`pip install seqdataloader`


## How to run 
Sample datasets are included in the folder `peak_files_from_encode_for_label_comparison` and `bigwig_files_from_encode_for_label_comparison`

Execute the script:

`genomewide_labels.sh` for examples on how to generate classification and regression labels on sample datasets.
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
Please note that the large bottleneck in the code is writing the files to disk. `hdf5` has negligible overhead, but using `gzip` or `bz2` may increase runtime. Timining benchmarks are provided in `genomewide_labels.sh`

You may speed up i/o by writing chromosome outputs to separate files in parallel. This is currently only supported for the `gzip` and `bz2` output types, as i/o is less of a bottleneck for `hdf5` and `pkl` output formats. Use the flag `--split_output_by_chrom` to invoke this parallelized saving of chromosomes.

## Documentation and benchmarks

Testing, benchmarks, and documentation can be found in the `docs` folder
