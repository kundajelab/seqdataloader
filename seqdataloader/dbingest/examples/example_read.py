import tiledb
import numpy as np
fname="/users/annashch/seqdataloader/seqdataloader/dbingest/encode_dnase_hg38/ENCFF001CTB.chr21"
with tiledb.DenseArray(fname,mode='r') as array:
    data=array[20000000:20001000]
    print(data)
    print(data['signal_value']) 

    
