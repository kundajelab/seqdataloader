#unit tests for class seqdataloader.batchproducers.coordbased.coordstovals.BasicTiledbProfileCoordsToVals
import pdb
from seqdataloader.batchproducers.coordbased.coordstovals.tiledb import *

#generate some test coords objects 
from collections import namedtuple
Coord=namedtuple('Coord','chrom start end isplusstrand')
coords=[Coord('chr1',1000000,2000000,True),
        Coord('chr2',1000000,2000000,True),
        Coord('chr3',1000000,2000000,True),
        Coord('chr4',1000000,2000000,True),
        Coord('chr5',1000000,2000000,True),
        Coord('chr6',1000000,2000000,True),
        Coord('chr7',1000000,2000000,True),
        Coord('chr1',1000000,2000000,False),
        Coord('chr2',1000000,2000000,False),
        Coord('chr3',1000000,2000000,False),
        Coord('chr4',1000000,2000000,False),
        Coord('chr5',1000000,2000000,False),
        Coord('chr6',1000000,2000000,False),
        Coord('chr7',1000000,2000000,False)]


pos_label_source_attribute="fc_bigwig"
neg_label_source_attribute="fc_bigwig"


#case 1: tiledb_paths is a string
tiledb_paths="/mnt/data/tiledb/encode/dnase/ENCSR000EOY"
ctov=BasicTiledbProfileCoordsToVals(tiledb_paths=tiledb_paths,
                                    pos_label_source_attribute=pos_label_source_attribute,
                                    neg_label_source_attribute=neg_label_source_attribute)
string_vals=ctov.__call__(coords)

