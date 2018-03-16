#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Helper functions (HDF5 files)
MIT License, copyright (c) 2018 Thomas Euler
"""
# ---------------------------------------------------------------------
__author__ = "code@eulerlab.de"

import h5py

# ---------------------------------------------------------------------
def read_hdf5_data(path, folder):
    return {key:item[:] for key, item in h5py.File(path, 'r')[folder].items()}

def nan2zero(data):
    """ Replace NaNs and negative values by zeros
    """
    return [val if val>0 else 0 for val in data]
    
# ---------------------------------------------------------------------
preLev = 0

def __print(path, obj):
    global preLev
    
    curLev = path.count("/")
    minInd = "   "
    indStr = minInd *curLev
    if isinstance(obj, h5py.Group):
        print("{0}{1}".format(indStr, obj.name))
        
    elif isinstance(obj, h5py.Dataset):    
        dataSetName = obj.name.split("/")[-1]
        if dataSetName == "_info":
            text = obj[()].decode("utf-8")
            print("{0}|  info='{1}'".format(indStr, text))
        else:
            print("{0}|--{1}{2}".format(indStr, dataSetName, obj.shape))
        
    if len(obj.attrs) > 0:
        text = obj.attrs["desc"][0].decode("utf-8")
        print("{0}desc='{1}'".format(indStr +minInd, text))
    
def print_hfd5_directory(path):    
    global preLev 
    
    preLev = 0    
    f = h5py.File(path, 'r')
    f.visititems(__print)
    
# ---------------------------------------------------------------------
