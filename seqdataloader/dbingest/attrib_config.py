from .utils import *

def get_attribute_info():
    required_attribs=['fc_bigwig','pval_bigwig','count_bigwig','idr_peak','overlap_peak','ambig_peak']

    attrib_info=dict()

    attrib_info['pval_bigwig']={'dtype':'float32',
                                'opener':open_bigwig_for_parsing,
                                'parser':parse_bigwig_chrom_vals,
                                'store_summits':False}

    attrib_info['fc_bigwig']={'dtype':'float32',
                              'opener':open_bigwig_for_parsing,
                              'parser':parse_bigwig_chrom_vals,
                              'store_summits':False}

    attrib_info['count_bigwig_plus_5p']={'dtype':'float32',
                                         'opener':open_bigwig_for_parsing,
                                         'parser':parse_bigwig_chrom_vals,
                                         'store_summits':False}
    
    attrib_info['count_bigwig_minux_5p']={'dtype':'float32',
                                          'opener':open_bigwig_for_parsing,
                                          'parser':parse_bigwig_chrom_vals,
                                          'store_summits':False}    
    attrib_info['idr_peak']={'dtype':'int',
                             'opener':open_csv_for_parsing,
                             'parser':parse_narrowPeak_chrom_vals,
                             'store_summits':True,
                             'summit_indicator':2}

    
    attrib_info['overlap_peak']={'dtype':'int',
                                 'opener':open_csv_for_parsing,
                                 'parser':parse_narrowPeak_chrom_vals,
                                 'store_summits':True,
                                 'summit_indicator':2}
    
    attrib_info['ambig_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':False}
    return attrib_info 
