import glob

import pandas as pd

def loadBed(file1_path):

    cnt = 0
    scorethres = 15
    ratiothres = 8
    poss = set()

    with open(file1_path, "r") as bed_file:
        for line in bed_file:
            columns = line.strip().split("\t")
            key0 = str(columns[0]) + ":" + str(columns[2])
            poss.add(key0)
            # key0 = str(columns[0]) + ":" + str(columns[1])
            # poss.add(key0)

    return poss


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


def loadBed3(file1_path):

    poss = set()
    with open(file1_path, "r") as bed_file:
        n=0
        for line in bed_file:

            if n==0:
                n+=1
                continue

            columns = line.strip().split("\t")
            key0 = str(columns[0]) + ":" + str(columns[1])
            # print(key0)
            poss.add(key0)
            n+=1

    return poss

def run1():

    file1_path = '/mnt/share/ueda/RNA004/resource/Y2.bed'
    poss2 = loadBed3(file1_path)
    print("BED2:", len(poss2))
    file1_path = '/mnt/share/ueda/RNA004/resource/Y.bed'
    poss2 = addbed(poss2,file1_path)
    print("BED2:", len(poss2))

    vcf1 = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.Pseudo.result.col29.bed"
    poss = loadBed(vcf1)

    overlap_poss = poss.intersection(poss2)

    print("BED1:", len(poss))
    print("BED2:", len(poss2))
    print(":", len(overlap_poss))

# run2()

import pandas as pd
from bisect import bisect_left


run1()
