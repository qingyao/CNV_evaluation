import sys
import subprocess
import os
from time import strftime

serieslist=[]
workingdir = os.getcwd()
with open('affy.tsv','r',encoding='utf-8') as f:
    allseries = f.readlines()
    for i,ser in enumerate(allseries):
    	if i % 5 == 0:
        	serieslist.append(ser.strip())

for series in sorted(serieslist):
  # if series.startswith('GSE'+ str(sys.argv[1])):
  print(strftime("%H:%M:%S"), series)
  subprocess.run('python3 {}/evaluation.py -s {} -c 1 -w ~/aroma/EvaluationPack'.format(workingdir, series), shell=True)
