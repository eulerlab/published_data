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
curLev = 0
nLine  = 0
sNew   = ""
sPrev  = ""

def __print(path, obj):
    global preLev, curLev, nLine, sNew, sPrev
    
    preLev = curLev
    curLev = path.count("/")
    minInd = "   "
    indStr = minInd *curLev
    sNew   = ""
    if isinstance(obj, h5py.Group):
        sNew = "{0}{1}".format(indStr, obj.name.split("/")[-1])
        
    elif isinstance(obj, h5py.Dataset):    
        dataSetName = obj.name.split("/")[-1]
        if dataSetName == "_info":
            text = obj[()].decode("utf-8")
            sNew = "{0}|  info='{1}'".format(indStr, text)
        else:
            sNew = "{0}|--{1}{2}".format(indStr, dataSetName, obj.shape)
        
    if len(obj.attrs) > 0:
        text = obj.attrs["desc"][0].decode("utf-8")
        sNew = sNew +"\r{0}desc='{1}'".format(indStr +minInd, text)
        
    if nLine > 0:
        if (curLev < preLev):
            sPrev = sPrev.replace("|", "`")
        print(sPrev)
    sPrev = sNew
    nLine += 1
        
    
def print_hfd5_directory(path):    
    """ Prints a simple directory of the file
    """
    global preLev, nLine, sPrev, sNew 
    
    preLev = 0
    curLev = 0
    nLine  = 0
    sNew   = ""
    sPrev  = ""
    f = h5py.File(path, 'r')
    f.visititems(__print)
    print(sNew.replace("|", "`"))
    
# ---------------------------------------------------------------------
