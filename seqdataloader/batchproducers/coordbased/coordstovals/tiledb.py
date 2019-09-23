import tiledb
import numpy as np
from .core import CoordsToVals

class BasicTiledbProfileCoordsToVals(CoordsToVals):
    def __init__(self, tiledb_paths, pos_label_source_attribute, neg_label_source_attribute=None, center_size_to_use=None, **kwargs):
        '''
        tiledb_paths can be a single string or a list of strings or a dictionary mapping from  mode name to string. 
        '''
        self.tiledb_paths=tiledb_paths
        #identify the data type of tiledb_paths 
        self.type_tiledb_paths=type(self.tiledb_paths)
        #identify the corresponding function to use for querying tiledb 
        self.call_function=self.get_call_function()
        #positive and negative strand values may correspond to different attirbutes of tiledb database
        self.pos_label_source_attribute=pos_label_source_attribute
        self.neg_label_source_attribute=neg_label_source_attribute

    def get_call_function(self):
        '''
        determines function to use for querying coord values based 
        on the data type of tiledb_paths attribute 
        '''
        if self.type_tiledb_paths == str:
            return self.__call__string
        elif self.type_tiledb_paths == list:
            return self.__call__list
        elif self.type_tiledb_paths == dict:
            return self.__call__dict
        else:
            raise Exception("Unsupported data type for BasicTiledbProfileCoordsToVals:"+str(self.type_tiledb_paths))
                
    def __call__dict(self,coords):
        '''
        self.tiledb_paths is a dictinary mapping from mode name to string
        '''
        vals={}
        for mode_name in self.tiledb_paths:
            cur_tiledb_path=self.tiledb_paths[mode_name]
            vals[mode_name]=self.query_tiledb(cur_tiledb_path,coords)
        return vals
    
    def __call__list(self,coords):
        '''
        self.tiledb_paths is a list of strings  
        '''
        vals=[self.query_tiledb(cur_tiled_path,coords) for cur_tiled_path in self.tiledb_paths]
        return vals 
    
    def __call__string(self,coords):
        '''
        self.tiledb_paths is a string 
        '''
        vals=self.query_tiledb(self.tiledb_paths,coords)
        return vals
    
    def __call__(self,coords):
        '''
        coords is a list of named tuples : .chrom, .start, .end, .isplusstrand    
        returns nparray of values associated with coordinates
        '''
        assert len(coords)>0        
        self.ctx = tiledb.Ctx()
        return self.call_function(coords)

    def query_tiledb(self,cur_tiledb_path,coords):
        '''
        queries tiledb database for a specific batch of coordinates for a single dataset/task. 
        '''
        labels=np.zeros((len(coords),coords[0].end-coords[0].start))
        for i in range(len(coords)):
            coord=coords[i]
            #open the tiledb for access in a pre-defined context 
            with tiledb.DenseArray('.'.join([cur_tiledb_path,coord.chrom]), mode='r',ctx=self.ctx) as cur_array:
                if coord.isplusstrand:
                    #query positive strand (or non-stranded entity)
                    cur_vals=cur_array[coord.start:coord.end][self.pos_label_source_attribute]
                else:
                    #query negative strand , make sure to reverse the values
                    cur_vals=cur_array[coord.start:coord.end][self.neg_label_source_attribute][::-1]
            labels[i]=cur_vals
        return labels
