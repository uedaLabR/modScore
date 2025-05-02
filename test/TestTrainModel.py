from nnmodel.AttentionClassfication import *
from nnmodel import AttentionClassfication

bases = ["A","C","T"]

for base in bases:

    checkpoint_path =  "/data/source/model_weight/hg38_"+base+".weights.h5"
    fp_ivtpath = "/share/ueda/nanoModiTune/U87IVTpu.bed"
    outhistory = "/share/ueda/nanoModiTune/outhistory_"+base+"hg38.csv"
    sourcepath="/share/ueda/nanoModiTune/resource"
    genome="hg38"
    AttentionClassfication.trainNN(base, sourcepath, genome, fp_ivtpath, outhistory, checkpoint_path)

for base in bases:

    checkpoint_path = "/data/source/model_weight/mm10_" + base + ".weights.h5"
    fp_ivtpath = "/share/ueda/nanoModiTune/U87IVTpu.bed"
    outhistory = "/share/ueda/nanoModiTune/outhistory_"+base+"mm10.csv"
    sourcepath = "/share/ueda/nanoModiTune/resource"
    genome = "mm10"
    AttentionClassfication.trainNN(base, sourcepath, genome, fp_ivtpath, outhistory, checkpoint_path)




