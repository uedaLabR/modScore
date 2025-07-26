from stats.StatsResult import stats_result

# Example usage:
if __name__ == "__main__":


    bed = "/mnt/ssdnas07/nanozero/rna/rna_modscore_v1/HEK293T_DR13/HEK293T_DR13/HEK293T_DR13_pileup_filter.bed"
    stats_f = "/share/ueda/nanoModiTune/Adipocyte_1_pileup_filter_stats.bed"
    stats_result(bed,stats_f)