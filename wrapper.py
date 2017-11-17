import sys
import subprocess
import os
from time import strftime

serieslist=[]
workingdir = os.getcwd()
with open('affy.tsv','r',encoding='utf-8') as f:
    allseries = f.readlines()
    for i in allseries:
        serieslist.append(i.strip())

for series in sorted(serieslist):
  if series.startswith('GSE'+ str(sys.argv[1])):
      print(strftime("%H:%M:%S"), series)
      subprocess.run('python3 {}/evaluation.py -s {}'.format(workingdir, series), shell=True)
