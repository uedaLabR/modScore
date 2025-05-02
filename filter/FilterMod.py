#Filter mod
from MSUtils import *

# Flg_Error = 0
# Flg_other = 1
# Flg_m6A = 2
# Flg_I = 3
# Flg_m5C = 4
# Flg_Y = 5
def strFlg(prediction):

    if prediction == Flg_Error:
        return "Flg_Error"
    if prediction == Flg_other:
        return "Flg_other"
    if prediction == Flg_m6A:
        return "Flg_m6A"
    if prediction == Flg_I:
        return "Flg_I"
    if prediction == Flg_m5C:
        return "Flg_m5C"
    if prediction == Flg_Y:
        return "Flg_Y"
    return "None"


def _posInKnownPos(knownPos,mod_flg,poskey):

    flg_indb = knownPos.get(poskey,-1)
    return flg_indb == mod_flg

def posInKnownPos(knownPos,mod_flg,chrom,pos):

    key0 = str(chrom) + ":" + str(pos-1)
    key1 = str(chrom) + ":" + str(pos)
    key2 = str(chrom) + ":" + str(pos+1)
    return _posInKnownPos(knownPos,mod_flg,key0) or _posInKnownPos(knownPos,mod_flg,key1) or _posInKnownPos(knownPos,mod_flg,key2)

from Configs import *
def filter_m5C(chrom_group,  knownPos):

    retlist = []
    depth_thres,ratio_thres = SecondndThresFor_m5C
    for chrom in chrom_group:

        chrlist = chrom_group[chrom]
        for columns,sequence,prediction in chrlist:

            pos = int(columns[1])
            score = int(columns[4])
            ratio = float(columns[10])

            flgOK = (prediction == Flg_m5C)
            secondthresOK = score > depth_thres and ratio > ratio_thres

            if posInKnownPos(knownPos,Flg_m5C,chrom,pos):
                columns.append("knownPos")
                columns.append("Pass")
            elif flgOK and secondthresOK:
                columns.append("flgOK")
                columns.append("Pass")
            else:
                if flgOK:
                    columns.append(strFlg(prediction)+",low thres")
                else:
                    columns.append(strFlg(prediction))
                columns.append("Fail")

            retlist.append(columns)

    return retlist

import re
def isDrach(qseq):

    # Define the regular expression pattern for DRACH motif
    drach_pattern = re.compile(r"[AGUT][AG]AC[AUCUT]", re.IGNORECASE)

    # Search for the pattern in the input sequence
    return bool(drach_pattern.search(qseq))


def getFlg_m6A(drach,inknownPos,predictionOK):

    flg=""
    if drach:
        flg = flg+"DRACH"
    if inknownPos:
        if len(flg) > 0:
            flg = flg+","
        flg = flg+"knownPos"
    if predictionOK:
        if len(flg) > 0:
            flg = flg+","
        flg = flg+"prediction"
    return flg

def filter_m6A(chrom_group,  knownPos):

    retlist = []
    for chrom in chrom_group:

        chrlist = chrom_group[chrom]
        for columns, sequence, prediction in chrlist:

            chrom = columns[0]
            pos = int(columns[1])
            key = str(chrom) + ":" + str(pos)
            fivemer = sequence[18:23]
            drach = isDrach(fivemer)
            inknownPos = posInKnownPos(knownPos, Flg_m6A, chrom, pos)
            predictionOK =  (prediction == Flg_m6A)


            if drach or inknownPos or predictionOK:

                flg = getFlg_m6A(drach,inknownPos,predictionOK)
                columns.append(flg)
                columns.append("Pass")

            else:

                columns.append(strFlg(prediction))
                columns.append("Fail")

            retlist.append(columns)

    return retlist


from bisect import bisect_left, bisect_right
def nearKnown(positions,pos,window=100):

    if not positions:
        return False
    left = bisect_left(positions, pos - window)
    right = bisect_right(positions, pos + window)
    return left < right



def load_repeat_from_rmsk_bisect(repeat_file_path,genome):

    # RepeatMasker`
    alu_dict = {}

    with open(repeat_file_path, 'r') as f:
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

            mouse = "mm" in genome
            if mouse:
                if rep_class == "LINE" or rep_class == "SINE":
                    if chrom not in alu_dict:
                        alu_dict[chrom] = []
                    alu_dict[chrom].append((start, end))
            else:
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

def inRepeat(repeatdef,pos):

    if repeatdef is None:
        return False

    starts = repeatdef['starts']
    ends = repeatdef['ends']

    # pos <= end index
    idx = bisect_left(ends, pos)
    if idx < len(starts) and starts[idx] <= pos < ends[idx]:
        return True
    elif idx > 0 and starts[idx - 1] <= pos < ends[idx - 1]:
        return True
    return False

import glob
import pandas as pd
def loadKnown(genome,knownPosDir):

    known_dict = {}
    files = glob.glob(knownPosDir+"/db1/*")
    # print(files)
    tgfile =""
    for file in files:
        if genome in file and "editing" in file:
              tgfile = file
    if len(tgfile)==0:
        print("No def file for RNA Editing")
    else:
        bed_df = pd.read_csv(tgfile, header=None, sep='\t',comment='#')
        for index, row in bed_df.iterrows():

            chrom = row[0]
            pos = int(row[1])
            if chrom not in known_dict:
                known_dict[chrom] = []
            known_dict[chrom].append(pos)

    files2 = glob.glob(knownPosDir + "/db2/*")
    tgfile = None
    for file in files2:
        if genome in file and "editing" in file:
              tgfile = file
    if tgfile:
        bed_df = pd.read_csv(tgfile, header=None, sep='\t',comment='#')
        for index, row in bed_df.iterrows():

            chrom = row[0]
            pos = int(row[1])
            if chrom not in known_dict:
                known_dict[chrom] = []
            known_dict[chrom].append(pos)


    #sort by pos
    for key in known_dict:
        known_dict[key] = sorted(known_dict[key])

    return known_dict


def getFlg_Inosine(inknownPos,inrepeat,nearKnownSites,flgOk):
    flg = ""
    if inknownPos:
        if len(flg) > 0:
            flg = flg + ","
        flg = flg + "knownPos"
    if inrepeat:
        if len(flg) > 0:
            flg = flg + ","
        flg = flg + "inrepeat"
    if nearKnownSites:
        if len(flg) > 0:
            flg = flg + ","
        flg = flg + "nearKnownSites"
    if flgOk:
        if len(flg) > 0:
            flg = flg + ","
        flg = flg + "flgOk"

    return flg


def filter_Inosine(chrom_group, repeat,genome,knownPos,source_path):

    repeat_dict = load_repeat_from_rmsk_bisect(repeat,genome)
    known_dict =  loadKnown(genome,source_path)
    retlist = []
    for chrom in chrom_group:

        repeatdef = repeat_dict.get(chrom, {'starts': [], 'ends': []})
        knowndef = known_dict.get(chrom, [])

        chrlist = chrom_group[chrom]
        for columns, sequence, prediction in chrlist:

            chrom = columns[0]
            pos = int(columns[1])
            key = str(chrom) + ":" + str(pos)
            inrepeat = inRepeat(repeatdef,pos)
            nearKnownSites = nearKnown(knowndef,pos)
            inknownPos = posInKnownPos(knownPos, Flg_I, chrom, pos)
            flgOk = (prediction == Flg_I)
            if inknownPos or (inrepeat and flgOk) or (nearKnownSites and flgOk):

                flg = getFlg_Inosine(inknownPos,inrepeat,nearKnownSites,flgOk)
                columns.append(flg)
                columns.append("Pass")

            else:
                if flgOk:
                    columns.append(strFlg(prediction)+",non near sites")
                else:
                    columns.append(strFlg(prediction))

            retlist.append(columns)

    return retlist


def getFlg_PsudeUridine(inknownPos, predictionOK):
    flg = ""
    if inknownPos:
        if len(flg) > 0:
            flg = flg + ","
        flg = flg + "knownPos"
    if predictionOK:
        if len(flg) > 0:
            flg = flg + ","
        flg = flg + "flgOk"

    return flg
def filter_PsudeUridine(chrom_group,  knownPos):

    retlist = []
    for chrom in chrom_group:

        chrlist = chrom_group[chrom]
        for columns, sequence, prediction in chrlist:

            chrom = columns[0]
            pos = int(columns[1])
            inknownPos = posInKnownPos(knownPos, Flg_Y, chrom, pos)
            predictionOK = (prediction == Flg_Y)
            if inknownPos or predictionOK:

                flg = getFlg_PsudeUridine(inknownPos, predictionOK)
                columns.append(flg)
                columns.append("Pass")

            else:

                columns.append(strFlg(prediction))
                columns.append("Fail")

            retlist.append(columns)



    return retlist
