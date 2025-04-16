from bam_manup import BamRecalib

inbam ="/mnt/ssdnas07/nanozero/rna/nanomoditune_v02/HEK293T_DR13/HEK293T_DR13/HEK293T_DR13_sorted.bam"
outbam="/share/ueda/nanoModiTune/Hek293_filter.bam"
filter_bed= "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
BamRecalib.run_recalib(inbam, outbam, filter_bed)