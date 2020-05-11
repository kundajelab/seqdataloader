from .utils import *

def get_multi_peak():
    required_attribs=['count_bigwig_plus_5p','count_bigwig_minus_5p','count_bigwig_unstranded_5p','100bp_peak','200bp_peak','500bp_peak','ambig_peak']
    attrib_info={}
    attrib_info['count_bigwig_plus_5p']={'dtype':'float32',
                                         'opener':open_bigwig_for_parsing,
                                         'parser':parse_bigwig_chrom_vals,
                                         'store_summits':False}

    attrib_info['count_bigwig_minus_5p']={'dtype':'float32',
                                          'opener':open_bigwig_for_parsing,
                                          'parser':parse_bigwig_chrom_vals,
                                          'store_summits':False}
    attrib_info['count_bigwig_unstranded_5p']={'dtype':'float32',
                                               'opener':open_bigwig_for_parsing,
                                               'parser':parse_bigwig_chrom_vals,
                                               'store_summits':False}
    attrib_info['100bp_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':True,
                               'summit_indicator':2,
                               'summit_from_peak_center':True}
    attrib_info['200bp_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':True,
                               'summit_indicator':2,
                               'summit_from_peak_center':True}
    attrib_info['500bp_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':True,
                               'summit_indicator':2,
                               'summit_from_peak_center':True}
    attrib_info['ambig_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':False,
                               'summit_from_peak_center':False}

    return attrib_info

def get_generic_bigwig_config():
    required_attribs=['bigwig_track'],

    attrib_info=dict()
    attrib_info['bigwig_track']={'dtype':'float32',
                                  'opener':open_bigwig_for_parsing,
                                  'parser':parse_bigwig_chrom_vals,
                                  'store_summits':False}
    return attrib_info 

def get_encode_config():
    required_attribs=['fc_bigwig','pval_bigwig','count_bigwig_plus_5p','count_bigwig_minus_5p','idr_peak','overlap_peak','ambig_peak']

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
    
    attrib_info['count_bigwig_minus_5p']={'dtype':'float32',
                                          'opener':open_bigwig_for_parsing,
                                          'parser':parse_bigwig_chrom_vals,
                                          'store_summits':False}
    attrib_info['count_bigwig_unstranded_5p']={'dtype':'float32',
                                          'opener':open_bigwig_for_parsing,
                                          'parser':parse_bigwig_chrom_vals,
                                          'store_summits':False}
    
    attrib_info['idr_peak']={'dtype':'int',
                             'opener':open_csv_for_parsing,
                             'parser':parse_narrowPeak_chrom_vals,
                             'store_summits':True,
                             'summit_indicator':2,
                             'summit_from_peak_center':False}

    
    attrib_info['overlap_peak']={'dtype':'int',
                                 'opener':open_csv_for_parsing,
                                 'parser':parse_narrowPeak_chrom_vals,
                                 'store_summits':True,
                                 'summit_indicator':2,
                                 'summit_from_peak_center':False}
    
    attrib_info['ambig_peak']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':False,
                               'summit_from_peak_center':False}
    return attrib_info 

def get_attribute_info(attribute_config):
    try:
        name_to_config=dict()
        name_to_config['encode_pipeline']=get_encode_config()
        name_to_config['generic_bigwig']=get_generic_bigwig_config()
        name_to_config['multi_peak']=get_multi_peak()
        attrib_info=name_to_config[attribute_config]
        return attrib_info
    except Exception as e:
        raise e