
if __name__ == "__main__":


    bed_file = "/mnt/ssdnas07/nanozero/rna/rna_modscore_v1/HEK293T_DR13/HEK293T_DR13/HEK293T_DR13_pileup_filter_only.bed"
    # bed_file = "/mnt/ssdnas07/nanozero/rna/rna_modscore_v1/HEK293T_DR13/HEK293T_DR13/HEK293T_DR13_pileup.bed"

    chroms = set()
    cnt = 0
    with open(bed_file, "r") as bed:

        for line in bed:

            data = line.split('\t')
            chr = data[0]
            # print(line)
            if "alt" in chr or "random" in chr:
                continue
            cnt+=1
            if cnt%1000==0:
                print(chr,cnt)
            chroms.add(chr)
    print(chroms)