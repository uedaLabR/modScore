from filter.FilterBed import *

# Example usage:
if __name__ == "__main__":

    # Define thresholds and file paths.
    # scorethres = 20
    # ratiothres = 10
    bed = "/share/ueda/nanoModiTune/Hek293pu.bed"
    bed_out = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    source_path = "/share/ueda/nanoModiTune/resource"
    # # bed = "/share/ueda/nanoModiTune/Hek293_recalib_pu.bed"
    # # bed_out = "/share/ueda/nanoModiTune/Hek293_recalib_pu_filter.bed"
    # ref = "/mnt/ssdnas07/pipeline/rna_v08/source/hg38.fa"
    # repeat = "/share/ueda/nanoModiTune/rmsk.txt"
    # # knownDB = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human*.bed"
    # knownDB = "/share/ueda/nanoModiTune/resource"
    # # db_m5C = "/share/ueda/nanoModiTune/db2"
    # genome = "hg38"
    # checkpoint_path = "/mnt/share/ueda/RNA004/nanoModFilter/ntmodel.weights.h5"
    genome = "hg38"

    filterBed(bed, bed_out, source_path, genome)

    # filterBed(bed, bed_out, ref, checkpoint_path, knownDB, genome)

