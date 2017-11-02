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
import os
from yaml import dump#, load 
from time import strftime
#try:
#    from yaml import CLoader as Loader, CDumper as Dumper
#except ImportError:
#    from yaml import Loader, Dumper


@click.command()
@click.option('-s','--series', default = "GSE40546")
@click.option('-t','--test', default = 0)
def cli(series,test):
    print(os.getcwd())
#    array = 'GSM996044'

    if test:
        data_dir = 'test_data/'
    else:
        data_dir = '/Volumes/arraymapMirror/arraymap/hg19/'

### run CNARA, per series
    arguments = [series, data_dir]
    subprocess.run('/usr/local/bin/r --vanilla </Users/pgweb/arraydata/aroma/EvaluationPack/run_CNARA.r --args {0} '.format(' '.join(arguments)), shell = True)  

### run analysis per array
#    for array in os.listdir(os.path.join(data_dir,series) ):
#        if not array.startswith('GSM'): next
#        try:
#            cnseg = pd.read_csv(os.path.join(data_dir,series,array,'segments,cn.tsv'),sep='\t')
#        except:
#            with open('err.txt','a+') as f:
#                f.write(' '.join([strftime("%Y-%m-%d %H:%M"),series,array,'no cn segments']) + '\n' )
#            next   
#        try:
#            fbseg = pd.read_csv(os.path.join(data_dir,series,array,'segments,fracb.tsv'),sep='\t')
#        except:
#            with open('err.txt','a+') as f:
#                f.write(' '.join([strftime("%Y-%m-%d %H:%M"),series,array,'no fracb segments']) + '\n')
#            next
#        
#        clen = len(cnseg)
#        blen = len(fbseg)
#    
#    
#        ### adjust by mean ###
#        cnprob = pd.read_csv(os.path.join(data_dir,series,array,'probes,cn.tsv'),sep='\t')
#        cnprob = cnprob.convert_objects(convert_numeric=True)
#        adjMean = (cnprob['VALUE'][~np.isnan(cnprob['VALUE'])]).mean()
#
#        cnprob['VALUE'] -= adjMean
##        with open(os.path.join(data_dir,series,array,'adjusted_parameter.tsv'), 'a+') as f:
##            f.write('\t'.join(['Mean', str(-adjMean)]))
#        
#        cnseg['value'] -= adjMean
#        sk = st.skew(cnseg['value'])
#        kt = st.kurtosis(cnseg['value'])
#        
#        poslen = sum(cnseg[cnseg['value'] > 0]['end'] - cnseg[cnseg['value'] > 0]['start'])
#        neglen = sum(cnseg[cnseg['value'] < 0]['end'] - cnseg[cnseg['value'] < 0]['start'])
#        posratio = poslen/(poslen + neglen)
#        
#        names = ['CNsegments','fracbsegments','skewness','kurtosis','posratio']
#        values = [clen,blen,sk,kt,posratio]
#        with open(os.path.join(data_dir,series,array,'sample_evaluation.tsv'), 'a+') as f:
#            for i in range(len(names)):
#                f.write('\t'.join([names[i], str(round(values[i],4))])+ '\n')
#    
#        
#        p = cnseg [cnseg['value'] >0]
#        np.average(p['value'],weights= p['end']-p['start'])
#        n = cnseg [cnseg['value'] <0]
#        np.average(n['value'],weights= n['end']-n['start'])
#        
#        default_param = ['GAINTHRESH', 'BASELINECORR', 'MINSEGSIZE', 'MAXSEGSIZE', 'MAXY', 
#                         'MINPROBES', 'IMGW', 'PLOTAREAH', 'DOTSCALE', 'FONTPX', 
#                         'GENEFONTSIZE', 'SEGSTROKE', 'CHROPLOT', 'PROBEMEAN', '-center', '-scale']
#        default_value = ['0.15', '0', '0', '250000000', '3.5', 
#                         '2', '780', '180', '1', '11',
#                         '11', '3', '1', '%.4f'%adjMean, 'n', 'n']
#        with open(os.path.join(data_dir,series,array,'defaults.yaml'),'w') as defaultfile:
#            dump(dict(zip(default_param,default_value)), defaultfile, default_flow_style=False)
        
        

####
#    if clen > 1000:
#        args = ['cnseg','4', '/Users/pgweb/arraydata/aroma/hg19', series, array, '/Users/pgweb/arraydata/aroma/AromaPack/Rpipeline/']
#        subprocess.run('/usr/local/bin/r --vanilla </Users/pgweb/arraydata/aroma/AromaPack/RchooseFunc.r --args {0}'.format(' '.join(args)), shell = True)
#    if blen > 1000: 
#        args = ['fracBseg','4', '/Users/pgweb/arraydata/aroma/hg19', series, array, '/Users/pgweb/arraydata/aroma/AromaPack/Rpipeline/']
#        subprocess.run('/usr/local/bin/r --vanilla </Users/pgweb/arraydata/aroma/AromaPack/RchooseFunc.r --args {0}'.format(' '.join(args)), shell = True)

    
def test_example():
    cli('sample_series',1)
    assert()

if __name__ == '__main__':
    print()
    cli()
