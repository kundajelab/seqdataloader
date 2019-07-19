import numpy as np
import sys
import tiledb

# Name of the array to create.
array_name = "quickstart_dense"

def create_array():

    # Check if the array already exists.
    # The array will be 4x4 with dimensions "rows" and "cols", with domain [1,4].
    dom = tiledb.Domain(tiledb.Dim(name="rows", domain=(1, 4), tile=4, dtype=np.int32),
                        tiledb.Dim(name="cols",domain=(1, 4), tile=4, dtype=np.int32))

    # The array will be dense with a single attribute "a" so each (i,j) cell can store an integer.
    schema = tiledb.ArraySchema(domain=dom, sparse=False,
                                attrs=[tiledb.Attr( name="a", dtype=np.int32)])

    # Create the (empty) array on disk.
    tiledb.DenseArray.create(array_name, schema)

def write_array():
    # Open the array and write to it.
    with tiledb.DenseArray(array_name, mode='w') as A:
        A[1:2] = np.array([1,1,1,1])

def read_array():
        # Open the array and read from it.
        with tiledb.DenseArray(array_name, mode='r') as A:
            # Slice only rows 1, 2 and cols 2, 3, 4.
            data = A[1:3, 2:5]
            print(data["a"])

if tiledb.object_type(array_name)!="array": 
    create_array()
write_array()
read_array()
