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
import click
import os
from time import strftime
import csv
from multiprocessing.dummy import Pool
import random


@click.command()
@click.option('-x','--test', default = 0)
@click.option('-u','--update',default = 0)
@click.option('-w','--workingdir',default = '~/aroma/EvaluationPack')
@click.option('-t','--threads',default = 5)
@click.option('-b','--block', default=0, help='Block number out of 5 blocks: 0--4')
@click.option('-a','--allblocks',default=5,help='Number of machines to pass to')
def cli(test,update,workingdir,threads,block,allblocks):
    p = Pool(processes=threads)
    serieslist=[]
    data_dir = '/Volumes/arraymapMirror/arraymap/hg19/'
    workingdir = os.path.expanduser(workingdir)
    with open(os.path.join(workingdir,'affy.tsv'),'r',encoding='utf-8') as f:
        c=0
        for line in f:
            if c%allblocks==block:
                if not line.startswith('G'):
                    serieslist.append('GSE'+line.rstrip())
                else:
                    serieslist.append(line.rstrip())
            c+=1

    if test:
        serieslist = random.sample(len(serieslist),10)
    serieslist = sorted(serieslist,reverse=True)
    print ('Doing %d series: \n%s' %(len(serieslist),serieslist))

#### run CNARA, per series
#    arguments = [series, data_dir, str(update)]
#    subprocess.run('cd {0};/usr/local/bin/r --vanilla <run_CNARA.r --args {1} '.format(workingdir,' '.join(arguments)), shell = True)

## run analysis per array
    summaryfile = 'output/sample_evaluation_summary.tsv'
    arrayPlatform = {}
    arrayTumor = {}
    with open('blockfile_PLATFORM_knownser.tsv', 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        for i in reader:
            arrayPlatform.update({i[1]:i[2]})
            arrayTumor.update({i[1]:i[3]})

    names = ['platform','TumorOrNormal','CNsegments','max1_cn','max2_cn','max3_cn','fracbsegments','max1_fb','max2_fb','max3_fb','skewness_seg','kurtosis_seg','skewness_prob','kurtosis_prob','probeMean','segMedian','positive_segment_ratio_mean','positive_segment_ratio_med']
    if not os.path.isfile(summaryfile):
        with open(summaryfile,'w') as evaluationfile:
            evaluationfile.write('\t'.join(['Series', 'Array'] + names) + '\n')

    ###to check if the array is new in file
    donearrays = []
    with open(summaryfile,'r') as evaluationfile:
        reader = csv.reader(evaluationfile,delimiter='\t')
        for i in reader:
            donearrays.append(i[1])
    donearrays = donearrays[1:]
    serieslist=['GSE11409']
    def runPipeline (series):
        print ('Started series: {}'.format(series))
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
            clen_chr = []
            try:
                for chr in range(23):
                    clen_chr.append(len(cnseg[cnseg.iloc[:,1]==chr+1]))
            except TypeError:
                for chr in range(23):
                    clen_chr.append(len(cnseg[cnseg.iloc[:,1]==str(chr+1)]))
            max3_clen_chr = sorted(clen_chr,reverse=True)[:3]

            blen = len(fbseg)
            blen_chr = []
            try:
                for chr in range(23):
                    blen_chr.append(len(fbseg[fbseg.iloc[:,1]==chr+1]))
            except TypeError:
                for chr in range(23):
                    blen_chr.append(len(fbseg[fbseg.iloc[:,1]==str(chr+1)]))
            max3_blen_chr = sorted(blen_chr,reverse=True)[:3]


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
            values = [platform,tumorOrNormal] + list(map(str,[clen]+max3_clen_chr+[blen]+max3_blen_chr)) + list(map(lambda x: '%.4f'%(x) ,[skseg,ktseg,skprob,ktprob,adjMean,adjMedian,posratio_mean,posratio_med]))
    #        with open(os.path.join(data_dir,series,array,'sample_evaluation.tsv'), 'a') as f:
    #            for i in range(len(names)):
    #                f.write('\t'.join([names[i], values[i]]) + '\n')

            with open(summaryfile,'a') as evaluationfile:
                evaluationfile.write('\t'.join([series, array] + values) + '\n')

            p = cnseg [cnseg['value'] >0]
            np.average(p['value'],weights= p['end']-p['start'])
            n = cnseg [cnseg['value'] <0]
            np.average(n['value'],weights= n['end']-n['start'])

    p.map(runPipeline,serieslist)

def test_example():
    cli('sample_series',1)
    assert()

def weighted_median(pddata,value_colname, weights_StartEndinTuple):
    pddata.sort_values(value_colname,inplace=True)
    cumsum = (pddata[weights_StartEndinTuple[1]] - pddata[weights_StartEndinTuple[0]]).cumsum()
    cutoff = (pddata[weights_StartEndinTuple[1]] - pddata[weights_StartEndinTuple[0]]).sum() / 2.0
    return (pddata[value_colname][cumsum >= cutoff].iloc[0])


##for tests:
#cnseg = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE40546/GSM996083/segments,cn.tsv',sep='\t')
#pddata= cnseg
#value_colname= 'value'
#weights_StartEndinTuple= ('start','end')
# cnprob = pd.read_csv('/Volumes/arraymapMirror/arraymap/hg19/GSE40546/GSM996083/probes,cn.tsv',sep='\t')
if __name__ == '__main__':
    print()
    cli()
