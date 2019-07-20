from .utils import *

def get_attribute_info():
    required_attribs=['fc_bigwig','pval_bigwig','count_bigwig','idr_peak','overlap_peak','ambig_peak']

    attrib_info=dict()
    attrib_info['pval_bigwig']={'dtype':'float32',
                                'opener':open_bigwig_for_parsing,
                                'parser':parse_bigwig_chrom_vals}
    attrib_info['fc_bigwig']={'dtype':'float32',
                              'opener':open_bigwig_for_parsing,
                              'parser':parse_bigwig_chrom_vals}
    attrib_info['count_bigwig']={'dtype':'float32',
                                 'opener':open_bigwig_for_parsing,
                                 'parser':parse_bigwig_chrom_vals}
    attrib_info['idr_peak']={'dtype':'int',
                             'opener':open_csv_for_parsing,
                             'parser':parse_narrowPeak_chrom_vals}
    attrib_info['overlap_peak']={'dtype':'int',
                                 'opener':open_csv_for_parsing,
                                 'parser':parse_narrowPeak_chrom_vals}
    attrib_info['ambig_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals}
    return attrib_info 
