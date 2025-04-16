import glob

import pandas as pd

def loadBed(file1_path,mod):

    cnt = 0
    scorethres = 15
    ratiothres = 8
    poss = set()

    with open(file1_path, "r") as bed_file:
        for line in bed_file:
            columns = line.strip().split("\t")
            alt = columns[3]
            score = int(columns[4])
            ratio = float(columns[10])
            if columns[-1] == "Fail":
                continue


            key0 = str(columns[0]) + ":" + str(columns[2])
            # print(key0)
            if alt == mod:
                poss.add(key0)

    return poss


def loadBed2(file1_path):

    poss = set()
    with open(file1_path, "r") as bed_file:
        for line in bed_file:
            columns = line.strip().split("\t")
            key0 = str(columns[0]) + ":" + str(columns[2])
            # print(key0)
            poss.add(key0)

    return poss

def loadBed2Puls1(file1_path):

    poss = set()
    with open(file1_path, "r") as bed_file:
        for line in bed_file:
            columns = line.strip().split("\t")
            key0 = str(columns[0]) + ":" + str(int(columns[2])+1)
            # print(key0)
            poss.add(key0)

    return poss

def loadBed3(file1_path):

    poss = set()
    with open(file1_path, "r") as bed_file:
        n=0
        for line in bed_file:

            if n==0:
                n+=1
                continue

            columns = line.strip().split("\t")
            key0 = str(columns[0]) + ":" + str(int(str(columns[1])))
            # print(key0)
            poss.add(key0)
            n+=1

    return poss


#run1()

def run1(vcf1):


    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    # vcf1 = "/share/ueda/nanoModiTune/Hek293_recalib_pu_filter.bed"
    poss = loadBed(vcf1,"a")
    file1_path = '/mnt/share/ueda/RNA004/resource/glory.bed'
    poss2 = loadBed2(file1_path)

    print("BED1:", len(poss))
    print("BED2:", len(poss2))

    overlap_count = 0
    for pos in poss:
        chrom, pos_str = pos.split(":")
        pos_num = int(pos_str)
        candidate_keys = {f"{chrom}:{pos_num}", f"{chrom}:{pos_num-1}", f"{chrom}:{pos_num+1}"}
        if candidate_keys & poss2:
            overlap_count += 1

    print("Overlap (allowing }1 base):", overlap_count)

def reverse_complement(seq):
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

def is_drach(seq):
    # 1: A, G, T
    # 2: A, G
    # 3: A
    # 4: C
    # 5: A, C, T
    return (seq[0] in "AGT" and
            seq[1] in "AG" and
            seq[2] == "A" and
            seq[3] == "C" and
            seq[4] in "ACT")

import pysam
def filterDRACH(poss2,ref,adjust=0):

    filtered = set()
    for pos_str in poss2:
        parts = pos_str.split(":")
        if len(parts) != 2:
            continue
        chrom = parts[0]
        try:
            pos = int(parts[1])+adjust
        except ValueError:
            continue

        # pysam  0-indexed gpApos 1-indexedA
        # S pos-1 O2B
        # : [pos-3, pos+2) => 5
        start = pos - 3
        end = pos + 2
        if start < 0:
            continue
        with pysam.FastaFile(ref) as fasta:
            try:
                motif = fasta.fetch(chrom, start, end).upper()
            except ValueError:
                continue

            if len(motif) != 5:
                continue

            # print(motif)
            # T  reverse complement `FbN
            if motif[2] == "A":
                if is_drach(motif):
                    filtered.add(pos_str)
                    # print("DRACH",pos_str,len(filtered))
            elif motif[2] == "T":
                rev_motif = reverse_complement(motif)
                if is_drach(rev_motif):
                    filtered.add(pos_str)
                    # print("DRACH",pos_str,len(filtered))

    return filtered


def run1_1(vcf1):

    ref = "/mnt/ssdnas07/pipeline/rna_v08/source/hg38.fa"
    vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    poss = loadBed(vcf1, "a")
    poss = filterDRACH(poss,ref)

    file1_path = '/mnt/share/ueda/RNA004/resource/glory.bed'

    poss2 = loadBed2(file1_path)
    poss2 = filterDRACH(poss2,ref)



    print("BED1:", len(poss))
    print("BED2:", len(poss2))

    overlap_count = 0
    for pos in poss:
        chrom, pos_str = pos.split(":")
        pos_num = int(pos_str)
        candidate_keys = {f"{chrom}:{pos_num}", f"{chrom}:{pos_num-1}", f"{chrom}:{pos_num+1}"}
        if candidate_keys & poss2:
            overlap_count += 1

    print("Overlap (allowing }1 base):", overlap_count)


def run2(vcf1):

    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    poss = loadBed(vcf1,"m")
    file1_path = '/mnt/share/ueda/RNA004/resource/m5C.bed'
    poss2 = loadBed3(file1_path)


    print("BED1:", len(poss))
    print("BED2:", len(poss2))

    overlap_count = 0
    for pos in poss:
        chrom, pos_str = pos.split(":")
        pos_num = int(pos_str)
        candidate_keys = {f"{chrom}:{pos_num}", f"{chrom}:{pos_num-1}", f"{chrom}:{pos_num+1}"}
        if candidate_keys & poss2:
            overlap_count += 1

    print("Overlap (allowing }1 base):", overlap_count)

# run2()
def addbed(posset,file1_path):

    with open(file1_path, "r") as file1_path:
        n=0
        for line in file1_path:

            if "NONE" in line:
                continue
            if len(line) < 10:
                continue

            # print(line)
            columns = line.strip().split("_")
            # print(columns)
            chr = columns[0]
            poss = columns[1]

            if "-" in poss:
                send = poss.split("-")
                s = int(send[0])
                e = int(send[1])
                key0 = str(chr) + ":" + str(s)
                posset.add(key0)
                key0 = str(chr) + ":" + str(e)
                posset.add(key0)

            else:
                key0 = str(chr) + ":" + str(columns[1])
                posset.add(key0)
            n+=1

    return posset

def run3(vcf1):

    file1_path = '/mnt/share/ueda/RNA004/resource/Y2.bed'
    poss2 = loadBed3(file1_path)
    print("BED2:", len(poss2))
    file1_path = '/mnt/share/ueda/RNA004/resource/Y.bed'
    poss2 = addbed(poss2,file1_path)
    print("BED2:", len(poss2))

    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    poss = loadBed(vcf1,"17802")

    print("BED1:", len(poss))
    print("BED2:", len(poss2))

    overlap_count = 0
    for pos in poss:
        chrom, pos_str = pos.split(":")
        pos_num = int(pos_str)
        candidate_keys = {f"{chrom}:{pos_num}", f"{chrom}:{pos_num-1}", f"{chrom}:{pos_num+1}"}
        if candidate_keys & poss2:
            overlap_count += 1

    print("Overlap (allowing }1 base):", overlap_count)


# run2()

def run4(vcf1):

    file1_path = '/mnt/share/ueda/RNA004/resource/slikseq.bed'
    poss2 = loadBed3(file1_path)

    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    # vcf1 = "/share/ueda/nanoModiTune/Hek293_recalib_pu_filter.bed"
    poss = loadBed(vcf1,"17596")

    print("BED1:", len(poss))
    print("BED2:", len(poss2))

    overlap_count = 0
    for pos in poss:
        chrom, pos_str = pos.split(":")
        pos_num = int(pos_str)
        candidate_keys = {f"{chrom}:{pos_num}", f"{chrom}:{pos_num-1}", f"{chrom}:{pos_num+1}"}
        if candidate_keys & poss2:
            overlap_count += 1

    print("Overlap (allowing }1 base):", overlap_count)

from bisect import bisect_left, bisect_right

def run5(vcf1):

    file1_path = '/mnt/share/ueda/RNA004/resource/slikseq.bed'
    poss2 = loadBed3(file1_path)

    # poss2: dict[chrom] -> sorted list of positions
    poss2_dict = {}
    for pos_str in poss2:
        chrom, pos = pos_str.split(":")
        pos = int(pos)
        if chrom not in poss2_dict:
            poss2_dict[chrom] = []
        poss2_dict[chrom].append(pos)

    for chrom in poss2_dict:
        poss2_dict[chrom].sort()

    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    # vcf1 = "/share/ueda/nanoModiTune/Hek293_recalib_pu_filter.bed"
    poss = loadBed(vcf1, "17596")

    overlap_count = 0

    for pos_str in poss:
        chrom, pos = pos_str.split(":")
        pos = int(pos)
        if chrom not in poss2_dict:
            continue

        # }100bp vfCfbNX
        positions = poss2_dict[chrom]
        left = bisect_left(positions, pos - 100)
        right = bisect_right(positions, pos + 100)

        if left < right:
            overlap_count += 1

    print("BED1:", len(poss))
    print("BED2:", sum(len(v) for v in poss2_dict.values()))
    print("}100bp overlap:", overlap_count)


import pandas as pd
from bisect import bisect_left


# Alu BED[hi`AF start/end \[gXgj
import pandas as pd
from bisect import bisect_left

def load_alu_from_rmsk_bisect(file_path):
    # RepeatMasker`
    alu_dict = {}

    with open(file_path, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue

            chrom = cols[5]
            try:
                start = int(cols[6])
                end = int(cols[7])
            except ValueError:
                continue

            rep_class = cols[11]
            rep_family = cols[10]

            if rep_class == "SINE":
                if chrom not in alu_dict:
                    alu_dict[chrom] = []
                alu_dict[chrom].append((start, end))

    for chrom in alu_dict:
        intervals = sorted(alu_dict[chrom])
        starts = [s for s, e in intervals]
        ends = [e for s, e in intervals]
        alu_dict[chrom] = {'starts': starts, 'ends': ends}

    return alu_dict


def is_alu(chrom, pos, alu_dict):
    if chrom not in alu_dict:
        return False

    starts = alu_dict[chrom]['starts']
    ends = alu_dict[chrom]['ends']

    idx = bisect_left(ends, pos)
    if idx < len(starts) and starts[idx] <= pos < ends[idx]:
        return True
    elif idx > 0 and starts[idx - 1] <= pos < ends[idx - 1]:
        return True

    return False


def is_alu_bisect(chrom, pos, alu_dict):
    if chrom not in alu_dict:
        return False
    starts = alu_dict[chrom]['starts']
    ends = alu_dict[chrom]['ends']

    # pos <= end index
    idx = bisect_left(ends, pos)
    if idx < len(starts) and starts[idx] <= pos < ends[idx]:
        return True
    elif idx > 0 and starts[idx - 1] <= pos < ends[idx - 1]:
        return True
    return False


# gp
alu_dict = load_alu_from_rmsk_bisect("/share/ueda/nanoModiTune/rmsk.txt")

def run6(vcf1,alu_dict):

    file1_path = '/mnt/share/ueda/RNA004/resource/slikseq.bed'
    poss2 = loadBed3(file1_path)
    alucnt1 = 0
    alucnt2 = 0

    # poss2: dict[chrom] -> sorted list of positions
    poss2_dict = {}
    for pos_str in poss2:
        chrom, pos = pos_str.split(":")
        pos = int(pos)
        result = is_alu_bisect(chrom, pos, alu_dict)
        if result:
            alucnt1+=1

    # vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
    poss = loadBed(vcf1, "17596")
    for pos_str in poss:
        chrom, pos = pos_str.split(":")
        pos = int(pos)
        result = is_alu_bisect(chrom, pos, alu_dict)
        if result:
            alucnt2+=1

    print("BED1:", alucnt2,len(poss))
    print("BED2:", alucnt1,len(poss2))



#
vcf1 = "/share/ueda/nanoModiTune/Hek293pu_filter.bed"
run1_1(vcf1)
run1(vcf1)
run2(vcf1)
run3(vcf1)
run4(vcf1)
run5(vcf1)
run6(vcf1,alu_dict)
