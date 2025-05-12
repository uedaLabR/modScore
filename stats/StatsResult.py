

def stats_result(bed_file, output_file):

    stats = {}
    with open(bed_file, "r") as bed:
        print(bed_file)
        for line in bed:
            print(line)
            columns = line.strip().split("\t")
            alt = columns[3]
            passfail = columns[-1]
            flg = columns[-2]

            #
            if alt not in stats:
                stats[alt] = {'Pass': 0, 'Fail': 0, 'Pass_flg': {}, 'Fail_flg': {}}

            # Pass/FailJEg
            stats[alt][passfail] += 1

            # flgJEg
            flg_dict = stats[alt][f'{passfail}_flg']
            if flg not in flg_dict:
                flg_dict[flg] = 0
            flg_dict[flg] += 1

    with open(output_file, "w") as out:
        out.write("alt\tStatus\tCount\tFlg\tFlg_Count\n")
        for alt, alt_data in stats.items():
            for status in ['Pass', 'Fail']:
                count = alt_data[status]
                flg_counts = alt_data[f'{status}_flg']
                if flg_counts:
                    for flg, flg_count in flg_counts.items():
                        print(f"{alt}\t{status}\t{count}\t{flg}\t{flg_count}\n")
                        out.write(f"{alt}\t{status}\t{count}\t{flg}\t{flg_count}\n")
                else:
                    print(f"{alt}\t{status}\t{count}\t-\t0\n")
                    out.write(f"{alt}\t{status}\t{count}\t-\t0\n")


# bed = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
# output = "/share/ueda/nanoModiTune/Hek293pu_stats.txt"
# #
# # # bed = "/share/ueda/nanoModiTune/Hek293_recalib_pu_filter.bed"
# # # output = "/share/ueda/nanoModiTune/Hek293pu_stats_recalib.txt"
# #
# stats_result(bed, output)
