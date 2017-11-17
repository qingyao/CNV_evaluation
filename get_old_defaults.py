#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:26:59 2017

@author: pgweb
"""
#import json
import csv
import os
import sys
import yaml

def write_threshold_to_one_file(lth,gth,file):
    ls_write = [series,arrayID,array[lth],array[gth]]
    file.write('\t'.join(ls_write) + '\n')
    
def write_probemean_to_one_file(probemean,file):
    ls_write = [series,arrayID,array[probemean]]
    file.write('\t'.join(ls_write) + '\n')

def write_ymax_to_one_file(ymax,file):
    ls_write = [series,arrayID,array[ymax]]
    file.write('\t'.join(ls_write) + '\n')


def write_thresholds_to_yaml(lth,gth):
    yaml_write = dict(cna_gain_threshold = array[gth], cna_loss_threshold = array[lth])
    try:
        with open(os.path.join(rootout,series,arrayID,'plotdefaults.yaml'),'w') as yamlf:
            yaml.dump(yaml_write, yamlf, default_flow_style=False)
    except FileNotFoundError:
        ef.write(series + '\t' + arrayID +"\tfile doesn't exist on .22")
        
def check_threshold(lth,gth):
    if lth in array.keys() or gth in array.keys():
        if float(array[lth]) != -0.15 and float(array[gth]) != 0.15:
#            print(float(array['LOSSTHRESH']),float(array['GAINTHRESH']),'float success')
            write_threshold_to_one_file(lth,gth,thf)
            write_thresholds_to_yaml(lth,gth)
        return('found')
    
def check_probemean(probemean):
    if probemean in array.keys():
        if abs(float(array[probemean])) > 0.05:
            write_probemean_to_one_file(probemean,pbf)
        return('found')
    
def check_ymax(ymax):
    if ymax in array.keys():
        if float(array[ymax]) < 3:
            write_ymax_to_one_file(ymax,ymf)
        return('found')
    
done_series=[]
if os.path.isfile('old_defaults_done_ser_68.txt'):
    with open('old_defaults_done_ser_68.txt', 'r') as df:
        content = df.readlines()
    for i in content:
        done_series.append(i.rstrip())

thf = open('old_threshold_68.txt', 'a')
pbf = open('old_probemean_68.txt','a')
ymf = open('old_ymax_68.txt','a')
ef = open('old_defaults_log_68.txt', 'a')
ff = open('old_defaults_file_err_68.txt', 'a')

root = '/Volumes/arrayRAID/arraymap/hg18/'
rootout = '/Volumes/arraymapMirror/arraymap/hg18/'
   
for series in sorted(os.listdir(root)):
    if done_series: 
        if series in done_series: continue
    try:
        arrayIDs = os.listdir(os.path.join(root,series))
        arrayIDs.sort()
        print(series,len(arrayIDs))
    except:
        continue
    for arrayID in arrayIDs:
#            print(arrayID)
        
        if arrayID == '.DS_Store': continue
        try:
            array={}
            with open(os.path.join(root,series,arrayID,'defaults.tab'), 'r') as f:
                reader = csv.reader(f, delimiter="\t")
                for i in reader:
                    if i:
                        array.update({i[0]:i[1]})
                
            ### check threshold
#            for lth,gth in [('LOSSTHRESH','GAINTHRESH'),('segLossThresh','segGainThresh'),('-lth','-gth')]:
#                if lth in array.keys() or gth in array.keys():
#                    try:
#                        if check_threshold(lth,gth) == 'found': continue
#                        
#                    except:
#                        raise
            ### check ymax
            trykey = 'not_found'
            for ymax in ['MAXY']:
                if ymax in array.keys():
                    try:
                        trykey = check_ymax(ymax)
                        if trykey == 'found': continue
                    except:
                        raise
            if trykey == 'not_found': 
                ef.write('\t'.join([series,arrayID,'ymax'])+'\t'+','.join(array.keys())+ '\n')
            
            ### check probemean
            trykey = 'not_found'
            for probemean in ['PROBEMEAN']:
                if probemean in array.keys():
                    try:
                        trykey = check_probemean(probemean)
                        if trykey == 'found': continue
                    except:
                        raise
            if trykey == 'not_found': 
                ef.write('\t'.join([series,arrayID,'probemean'])+'\t'+','.join(array.keys())+ '\n')
            
        except IndexError: ### This is when the files are separated by : and comma endline
            
            try:
                array={}
                with open(os.path.join(root,"GSE9672","GSM244402",'defaults.tab'), 'r') as f:
                    reader = csv.reader(f, delimiter=":")
                    for i in reader:
                        if i:
                            array.update({i[0]:i[1].replace(',','')})
                ### check threshold
#                for lth,gth in [('LOSSTHRESH','GAINTHRESH'),('segLossThresh','segGainThresh'),('-lth','-gth')]:
#                    if lth in array.keys() or gth in array.keys():
#                        try:
#                            if check_threshold(lth,gth) == 'found': continue
#                        except:
#                            raise  
                ### check ymax
                trykey = 'not_found'
                for ymax in ['MAXY']:
                    if ymax in array.keys():
                        try:
                            trykey = check_ymax(ymax)
                            if trykey == 'found': break
                        except:
                            raise
                if trykey == 'not_found': 
                    ef.write('\t'.join([series,arrayID,'ymax'])+'\t'+','.join(array.keys())+ '\n')
                
                ### check probemean
                trykey = 'not_found'
                for probemean in ['PROBEMEAN','BASECORR','segNormal']:
                    if probemean in array.keys():
                        try:
                            trykey = check_probemean(probemean)
                            if trykey == 'found': break
                        except:
                            raise
                if trykey == 'not_found': 
                    ef.write('\t'.join([series,arrayID,'probemean'])+'\t'+','.join(array.keys())+ '\n')
            except:
                raise
    
        except NotADirectoryError:
            ff.write('\t'.join([series,arrayID,'NotADirectoryError'])+'\n')
        except FileNotFoundError:
            ff.write('\t'.join([series,arrayID,'FileNotFoundError'])+'\n')
        except:
            ef.write(series + '\t' + arrayID +'\tUnexpected error:' + str(sys.exc_info()[0]) + '\n')
    with open('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_done_ser_68.txt', 'a') as df:
        df.write(series + '\n')

thf.close()
pbf.close()
ymf.close()

ef.close()

ff.close()
