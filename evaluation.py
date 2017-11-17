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
#import subprocess
import click
import os
#from yaml import dump#, load
from time import strftime
import csv

#try:
#    from yaml import CLoader as Loader, CDumper as Dumper
#except ImportError:
#    from yaml import Loader, Dumper


@click.command()
@click.option('-s','--series', default = "GSE40546")
@click.option('-t','--test', default = 0)
@click.option('-u','--update',default = 0)
def cli(series,test,update):
#    array = 'GSM996044'

    if test:
        data_dir = 'test_data'
    else:
        data_dir = '/Volumes/arraymapMirror/arraymap/hg19/'

#### run CNARA, per series
#    arguments = [series, data_dir, str(update)]
#    subprocess.run('/usr/local/bin/r --vanilla </Users/pgweb/arraydata/aroma/EvaluationPack/run_CNARA.r --args {0} '.format(' '.join(arguments)), shell = True)

## run analysis per array
    arrayPlatform = {}
    arrayTumor = {}
    with open('blockfile_PLATFORM_knownser.tsv', 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        for i in reader:
            arrayPlatform.update({i[1]:i[2]})
            arrayTumor.update({i[1]:i[3]})

    names = ['platform','TumorOrNormal','CNsegments','fracbsegments','skewness_seg','kurtosis_seg','skewness_prob','kurtosis_prob','probeMean','segMedian','positive_segment_ratio_mean','positive_segment_ratio_med']
    if not os.path.isfile('sample_evalutaion_summary.tsv'):
        with open('sample_evalutaion_summary.tsv','w') as evaluationfile:
            evaluationfile.write('\t'.join(['Series', 'Array'] + names) + '\n')

    ###to check if the array is new in file
    donearrays = []
    with open('sample_evalutaion_summary.tsv','r') as evaluationfile:
        reader = csv.reader(evaluationfile,delimiter='\t')
        for i in reader:
            donearrays.append(i[1])
    donearrays = donearrays[1:]

    for array in os.listdir(os.path.join(data_dir,series) ):

        if not array.startswith('GSM'): continue
        if array in donearrays: continue
        try:
            cnseg = pd.read_csv(os.path.join(data_dir,series,array,'segments,cn.tsv'),sep='\t')
        except:
            with open('err.txt','a+') as f:
                f.write(' '.join([strftime("%Y-%m-%d %H:%M"),series,array,'no data or problem parsing cn segments']) + '\n' )
            continue
        try:
            fbseg = pd.read_csv(os.path.join(data_dir,series,array,'segments,fracb.tsv'),sep='\t')
        except:
            with open('err.txt','a+') as f:
                f.write(' '.join([strftime("%Y-%m-%d %H:%M"),series,array,'no data or problem parsing fracb segments']) + '\n')
            continue

        clen = len(cnseg)
        blen = len(fbseg)


        ### adjust by mean ###
        try:
            cnprob = pd.read_csv(os.path.join(data_dir,series,array,'probes,cn.tsv'),sep='\t')
        except:
            with open('err.txt','a+') as f:
                f.write(' '.join([strftime("%Y-%m-%d %H:%M"),series,array,'no data or problem parsing cn probes']) + '\n' )
            continue

        cnprob['VALUE'] = pd.to_numeric(cnprob['VALUE'])
        adjMean = (cnprob['VALUE'][~np.isnan(cnprob['VALUE'])]).mean()
        adjMedian = weighted_median(cnseg,'value',('start','end'))

#        cnprob['VALUE'] -= adjMean
#        with open(os.path.join(data_dir,series,array,'adjusted_parameter.tsv'), 'a+') as f:
#            f.write('\t'.join(['Mean', str(-adjMean)]))


        skseg = st.skew(cnseg['value'])
        ktseg = st.kurtosis(cnseg['value'])

        skprob = float(np.ma.getdata(st.skew(cnprob['VALUE'],nan_policy='omit')))
        ktprob = float(np.ma.getdata(st.kurtosis(cnprob['VALUE'],nan_policy='omit')))

        cnseg['value'] -= adjMean
        poslen = sum(cnseg[cnseg['value'] > 0]['end'] - cnseg[cnseg['value'] > 0]['start'])
        neglen = sum(cnseg[cnseg['value'] < 0]['end'] - cnseg[cnseg['value'] < 0]['start'])
        posratio_mean = poslen/(poslen + neglen)

        cnseg['value'] = cnseg['value'] - adjMedian + adjMean
        poslen = sum(cnseg[cnseg['value'] > 0]['end'] - cnseg[cnseg['value'] > 0]['start'])
        neglen = sum(cnseg[cnseg['value'] < 0]['end'] - cnseg[cnseg['value'] < 0]['start'])
        posratio_med = poslen/(poslen + neglen)

        platform = arrayPlatform[array]
        tumorOrNormal =arrayTumor[array]
        values = [platform,tumorOrNormal] + list(map(str,[clen,blen])) + list(map(lambda x: '%.4f'%(x) ,[skseg,ktseg,skprob,ktprob,adjMean,adjMedian,posratio_mean,posratio_med]))
#        with open(os.path.join(data_dir,series,array,'sample_evaluation.tsv'), 'a') as f:
#            for i in range(len(names)):
#                f.write('\t'.join([names[i], values[i]]) + '\n')

        with open('sample_evalutaion_summary.tsv','a') as evaluationfile:
            evaluationfile.write('\t'.join([series, array] + values) + '\n')

        p = cnseg [cnseg['value'] >0]
        np.average(p['value'],weights= p['end']-p['start'])
        n = cnseg [cnseg['value'] <0]
        np.average(n['value'],weights= n['end']-n['start'])

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

def weighted_median(pddata,value_colname, weights_StartEndinTuple):
    pddata.sort_values(value_colname,inplace=True)
    cumsum = (pddata[weights_StartEndinTuple[1]] - pddata[weights_StartEndinTuple[0]]).cumsum()
    cutoff = (pddata[weights_StartEndinTuple[1]] - pddata[weights_StartEndinTuple[0]]).sum() / 2.0
    return (pddata[value_colname][cumsum >= cutoff].iloc[0])


##for tests:
cnseg = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE40546/GSM996083/segments,cn.tsv',sep='\t')
#pddata= cnseg
#value_colname= 'value'
#weights_StartEndinTuple= ('start','end')
cnprob = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE40546/GSM996083/probes,cn.tsv',sep='\t')
if __name__ == '__main__':
    print()
    cli()
