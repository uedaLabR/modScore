from AttentionClassfication import *

checkpoint_path =  "/mnt/share/ueda/RNA004/nanoModFilter/ntmodel.weights.h5"
# m6Apath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.m6A.result.col29.bed"
# m5Cpath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.m5C.result.col29.bed"
# psudepath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.Pseudo.result.col29.bed"
# editingpath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.RNA-editing.result.col29.bed"

fp_ivtpath = "/share/ueda/nanoModiTune/U87IVTpu.bed"

outhistory = "/share/ueda/nanoModiTune/outhistory_mm10.csv"
sourcepath="/share/ueda/nanoModiTune/resource"
# genome="hg38"
genome="mm10"

trainNN(sourcepath, genome, fp_ivtpath, outhistory, eachsize=20000, epoch=100)