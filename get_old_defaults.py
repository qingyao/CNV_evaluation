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

def write_threshold_to_one_file(lth,gth):
    ls_write = [series,arrayID,array[lth],array[gth]]
    wf.write('\t'.join(ls_write) + '\n')

def write_thresholds_to_yaml(lth,gth):
    yaml_write = dict(cna_gain_threshold = gth, cna_loss_threshold = lth)
    with open(os.path.join(root2,series,arrayID,'plotdefaults.yaml'),'w') as yamlf:
        yaml.dump(yaml_write, yamlf, default_flow_style=False)

done_series=[]
if os.path.isfile('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_done_ser_68.txt'):
    with open('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_done_ser_68.txt', 'r') as df:
        content = df.readlines()
    for i in content:
        done_series.append(i.rstrip())

wf = open('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_68.txt', 'a')
ef = open('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_log_68.txt', 'a')
ff = open('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_file_err_68.txt', 'a')

root = '/Volumes/arrayRAID/arraymap/hg18/'
root2 = '/Volumes/arraymapMirror/arraymap/hg18/'
   
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
            with open(os.path.join(root,series,arrayID,'defaults.tab'), 'r') as f:
                reader = csv.reader(f, delimiter="\t")
                array = dict(reader)
                
#            with open(os.path.join(root,series,arrayID,'progenetix.json'), 'r') as f:
#                array = json.load(f)
            
            if 'LOSSTHRESH' in array.keys() or 'GAINTHRESH' in array.keys():
                try:
                    if float(array['LOSSTHRESH']) != -0.15 and float(array['GAINTHRESH']) != 0.15:
                        print(float(array['LOSSTHRESH']),float(array['GAINTHRESH']),'float success')
                        write_threshold_to_one_file('LOSSTHRESH','GAINTHRESH')
                except:
                    raise
            elif 'segLossThresh' in array.keys() or 'segGainThresh' in array.keys():
                try:
                    if float(array['LOSSTHRESH']) != -0.15 and float(array['GAINTHRESH']) != 0.15:
                        write_threshold_to_one_file('segLossThresh','segGainThresh')
                except:
                    raise
            else:
                ef.write('\t'.join([series,arrayID])+'\t'+','.join(array.keys())+ '\n')
                
        except NotADirectoryError:
            ff.write('\t'.join([series,arrayID,'NotADirectoryError'])+'\n')
        except FileNotFoundError:
            ff.write('\t'.join([series,arrayID,'FileNotFoundError'])+'\n')

        except:
            ef.write(series + '\t' + arrayID +'\tUnexpected error:' + str(sys.exc_info()[0]) + '\n')
    with open('/Users/pgweb/arraydata/aroma/EvaluationPack/old_defaults_done_ser_68.txt', 'a') as df:
        df.write(series + '\n')

wf.close()

ef.close()

ff.close()



    
    