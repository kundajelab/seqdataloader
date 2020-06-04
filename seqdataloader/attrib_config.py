from .utils import *
allowed_attributes={}
allowed_attributes['bigwig']={'dtype':'float32',
                                         'opener':open_bigwig_for_parsing,
                                         'parser':parse_bigwig_chrom_vals,
                                         'store_summits':False}
allowed_attributes['bed_no_summit']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':False,
                               'summit_from_peak_center':False}
allowed_attributes['bed_summit_from_peak_center']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':True,
                               'summit_indicator':2,
                               'summit_from_peak_center':True}
allowed_attributes['bed_summit_from_last_col']={'dtype':'int',
                               'opener':open_csv_for_parsing,
                               'parser':parse_narrowPeak_chrom_vals,
                               'store_summits':True,
                               'summit_indicator':2,
                               'summit_from_peak_center':False}

def get_generic_bigwig_config():
    attrib_info=dict()
    attrib_info['bigwig_track']=allowed_attributes['bigwig']
    attrib_info['ambig_peak']=allowed_attributes['bed_no_summit']
    return attrib_info 

def get_encode_with_controls_config():
    attrib_info=get_encode_config()
    #add the control count tracks 
    attrib_info['control_count_bigwig_unstranded_5p']=allowed_attributes['bigwig']
    attrib_info['control_count_bigwig_plus_5p']=allowed_attributes['bigwig']
    attrib_info['control_count_bigwig_minus_5p']=allowed_attributes['bigwig']
    return attrib_info
    

def get_encode_config():
    attrib_info=dict()

    attrib_info['pval_bigwig']=allowed_attributes['bigwig']
    attrib_info['fc_bigwig']=allowed_attributes['bigwig']
    attrib_info['count_bigwig_plus_5p']=allowed_attributes['bigwig']
    attrib_info['count_bigwig_minus_5p']=allowed_attributes['bigwig']
    attrib_info['count_bigwig_unstranded_5p']=allowed_attributes['bigwig']
    attrib_info['idr_peak']=allowed_attributes['bed_summit_from_last_col']
    attrib_info['overlap_peak']=allowed_attributes['bed_summit_from_last_col']
    attrib_info['ambig_peak']=allowed_attributes['bed_no_summit']
    return attrib_info 

def get_attribute_info_from_file(attribute_config_file):
    config_metadata=open(attribute_config_file,'r').read().strip().split('\n')
    attrib_info={}
    for line in config_metadata:
        tokens=line.split('\t')
        field_name=tokens[0]
        field_type=tokens[1]
        attrib_info[field_name]=allowed_attributes[field_type]
    return attrib_info

def get_attribute_info(attribute_config,attribute_config_file):
    assert (attribute_config is None) or (attribute_config_file is None)
    if attribute_config_file is not None:
        return get_attribute_info_from_file(attribute_config_file) 
    try:
        name_to_config=dict()
        name_to_config['encode_pipeline_with_controls']=get_encode_with_controls_config()
        name_to_config['encode_pipeline']=get_encode_config()
        name_to_config['generic_bigwig']=get_generic_bigwig_config()
        attrib_info=name_to_config[attribute_config]
        return attrib_info
    except Exception as e:
        raise e
