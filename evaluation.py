#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 14:53:48 2017

@author: pgweb
"""
####

import scipy.stats as st
import pandas as pd
import numpy as np
import subprocess
import click
#from rpy2.robjects.packages import importr
#from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
#with open('/Users/pgweb/arraydata/aroma/AromaPack/Rpipeline/final_doCRMAv2.R','r') as f:
#    content = f.readlines()
#funcstring = '\n'.join([i.strip() for i in content])
#aromaR = SignatureTranslatedAnonymousPackage(funcstring, "aromaR")

@click.command()
@click.options('-s','--series', default = "GSE59150")
def cli(series):
    cnseg = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE59150/GSM1442240/segments,cn.tsv',sep='\t')
    clen = len(cnseg)
    fbseg = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE59150/GSM1442240/segments,fracb.tsv',sep='\t')
    blen = len(fbseg)
    
    series = 'GSE59150'
    array = 'GSM1442240'
    
    if clen > 1000:
        args = ['cnseg','4', '/Users/pgweb/arraydata/aroma/hg19', series, array, '/Users/pgweb/arraydata/aroma/AromaPack/Rpipeline/']
        subprocess.run('/usr/local/bin/r --vanilla </Users/pgweb/arraydata/aroma/AromaPack/RchooseFunc.r --args {0}'.format(' '.join(args)), shell = True)
    if blen > 1000: 
        args = ['fracBseg','4', '/Users/pgweb/arraydata/aroma/hg19', series, array, '/Users/pgweb/arraydata/aroma/AromaPack/Rpipeline/']
        subprocess.run('/usr/local/bin/r --vanilla </Users/pgweb/arraydata/aroma/AromaPack/RchooseFunc.r --args {0}'.format(' '.join(args)), shell = True)
        
    ### adjust by median ###
    cnprob = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE59150/GSM1442240/probes,cn.tsv',sep='\t')
    adjMed = np.average(cnprob['value'])
    cnprob['value'] = cnprob['value'] - adjMed
    
    
    cnseg['value'] = cnseg['value'] - adjMed
    st.skew(cnseg['value'])
    st.kurtosis(cnseg['value'])
    
    
    
    p = cnseg [cnseg['value'] >0]
    np.average(p['value'],weights= p['end']-p['start'])
    n = cnseg [cnseg['value'] <0]
    np.average(n['value'],weights= n['end']-n['start'])
    
