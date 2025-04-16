from bam_manup import BamRecalib
import cProfile
import pstats

inbam ="/mnt/ssdnas07/nanozero/rna/nanomoditune_v02/HEK293T_DR13/HEK293T_DR13/HEK293T_DR13_sorted.bam"
outbam="/share/ueda/nanoModiTune/Hek293_filter.bam"
filter_bed= "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
def main():
    # vt@CoLq
    BamRecalib.run_recalib(inbam, outbam, filter_bed)

main()

# profile_filename = "profile_results.prof"
# cProfile.run("main()", filename=profile_filename)
# stats = pstats.Stats(profile_filename)
# stats.sort_stats("cumulative").print_stats(30)