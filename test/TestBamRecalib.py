from bam_manup import BamRecalib
import cProfile
import pstats

inbam ="/mnt/ssdnas07/nanozero/rna/rna_modscore_v1/Adipocyte_1/Adipocyte_1/Adipocyte_1_sorted.bam"
# inbam="/share/ueda/nanoModiTune/Adipocyte_1_test_filter.bam"
outbam="/share/ueda/nanoModiTune/Adipocyte_1_test_filter2.bam"
filter_bed= "/mnt/ssdnas07/nanozero/rna/rna_modscore_v1/Adipocyte_1/Adipocyte_1/Adipocyte_1_pileup_filter.bed"
ref = "/mnt/ssdnas07/pipeline/rna_modscore_v1/source/genome/mm10.fa"

def main():

    #
    BamRecalib.run_recalib(inbam, outbam, filter_bed,ref)

main()

# profile_filename = "profile_results.prof"
# cProfile.run("main()", filename=profile_filename)
# stats = pstats.Stats(profile_filename)
# stats.sort_stats("cumulative").print_stats(30)