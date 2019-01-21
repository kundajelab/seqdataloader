import numpy as np

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def merge_dictionaries(main_dict,new_dict):
    for cur_bin in new_dict:
        if cur_bin not in main_dict:
            main_dict[cur_bin]=new_dict[cur_bin]
        else:
            for task_name in new_dict[cur_bin]:
                main_dict[cur_bin][task_name]=new_dict[cur_bin][task_name]
    return main_dict
