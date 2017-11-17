#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 08:13:02 2017

@author: pgweb
"""

import subprocess
import os
import pandas as pd
import numpy as np

# print(os.getcwd())
for array in os.listdir('test_data/sample_series'):
    if array.startswith('GSM'):
        # subprocess.run('cp {0}/probes,cn.tsv {0}/probes_copy,cn.tsv'.format(os.path.join('test_data/sample_series',array)) , shell=True) #save a copy
        # subprocess.run('cp {0}/probes_copy,cn.tsv {0}/probes,cn.tsv'.format(os.path.join('test_data/sample_series',array)) , shell=True) #revert
        f = pd.read_csv(os.path.join('test_data/sample_series',array,'probes,cn.tsv'), sep='\t',skiprows=[0],header=None,dtype={0:np.int32,1:np.character,2:np.int32,3:np.float64})
        f.columns = ['ID','chro','pos','log2']
        f['ID'] = list(map(lambda x: 'ID_'+str(x),f['ID']))
        f.to_csv(os.path.join('test_data/sample_series',array,'probes,cn.tsv'),sep='\t',index=False)
