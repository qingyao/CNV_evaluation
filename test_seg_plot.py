import os,subprocess

workdir = "/Users/pgweb/arraydata/aroma/EvaluationPack/test_data/combine_series"
# samples = os.listdir(workdir)
samples = ["GSM1704973","GSM1704975","GSM1704979","GSM1146803","GSM337641", "GSM381297", "GSM433918"]
samples = ['GSM1974257','GSM1974269']
for sample in samples:
# sample = 'GSM1146803'
    for lmd in [0.3,0.5,1,1.5,2,3,6]: #0.5,1
        for gapsize in ['1e4','5e4']: #,'1e5'
            if not os.path.isfile('{0}/segments,cn,5_sdundo_1,lasso{1}_{2}.tsv'.format(os.path.join(workdir,sample),lmd,gapsize)):
                continue
            if int(lmd) != lmd:
                lmdName = int(lmd*10)
            else:
                lmdName = lmd
            # subprocess.run ("perl /Users/pgweb/Dropbox\ \(baudisgroup\)/baudisgroup/dbtools/PGX/Scripts/arrayplotter.pl -genome hg19 -arraypath {0} -value_plot_y_max 3.2 -segfilename segments,cn,{1}_sdundo_{2}.tsv -svgfilename arrayplot,{1},{2}.svg".format(os.path.join(workdir,sample),smrg,undosd),shell=True)
            subprocess.run ("perl /Users/pgweb/Dropbox\ \(baudisgroup\)/baudisgroup/dbtools/PGX/Scripts/arrayplotter.pl -genome hg19 -arraypath {0} -value_plot_y_max 3.2 -segfilename segments,cn,5_sdundo_1,lasso{1}_{2}.tsv -svgfilename arrayplot,{3},{2}.svg".format(os.path.join(workdir,sample),lmd,gapsize,lmdName),shell=True)
