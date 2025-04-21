import glob
import pandas as pd

Flg_Error = 0
Flg_other = 1
Flg_m6A = 2
Flg_I = 3
Flg_m5C = 4
Flg_Y = 5

def loadKnownPos(knownPosDir,genome):


    knownPos = {}
    files = glob.glob(knownPosDir+"/db1/*")
    # print(files)
    for file in files:

        if genome in file:

            _loadKnownPos(knownPos, file)

    return knownPos

def _loadKnownPos(knownPos,path):

    flg = getFlg(path)
    bed_df = pd.read_csv(path, header=None, sep='\t')
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[1]
        key = str(chr) + ":" + str(pos)
        knownPos[key] = flg

def getFlg(path):

    if "m5C" in path:
        return Flg_m5C
    if "m6A" in path:
        return Flg_m6A
    if "Pseudo" in path:
        return Flg_Y
    if "editing" in path:
        return Flg_I
    return -1



def addData2(knownPos, path, flg,adjust=0):

    if not path or len(path) < 2:
        return

    with open(path, 'r') as f:

        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            if "\t" not in line:
                print(line)
                chrom_startend = line.split("_")
                if len(chrom_startend) > 2 or "NONE" in line:
                    continue
                chrom, startend = chrom_startend
                #
                if "-" in startend:
                    print(startend)
                    startend = startend.split("-")
                    key1 = chrom + ":" + str(int(startend[0])+adjust)
                    knownPos[key1] = flg
                    key2 = chrom + ":" + str(int(startend[1])+adjust)
                    knownPos[key2] = flg

                else:
                    key1 = chrom + ":" + str(int(startend)+adjust)
                    knownPos[key1] = flg

            else:

                cols = line.split("\t")
                key1 = cols[0] + ":" + str(int(cols[1])+adjust)
                knownPos[key1] = flg


def getFiles2(sourcepath,genome):

    m6Apath, m5Cpath, psudepath, editingpath = "","","",""
    files = glob.glob(sourcepath+"/db2/*")
    print("files",files,sourcepath+"/db2/*")
    for file in files:
        if genome in file:
            flg = getFlg(file)
            if flg == Flg_m5C:
                m5Cpath = file
            if flg == Flg_m6A:
                m6Apath = file
            if flg == Flg_Y:
                psudepath = file
            if flg == Flg_I:
                editingpath = file

    return m6Apath, m5Cpath, psudepath, editingpath
def addOtherDB(knownPos, sourcepath,genome):

    m6Apath, m5Cpath, psudepath, editingpath = getFiles2(sourcepath, genome)
    print((m6Apath, m5Cpath, psudepath, editingpath))
    if m6Apath:
        addData2(knownPos, m6Apath, Flg_m6A)
    if m5Cpath:
        addData2(knownPos, m5Cpath, Flg_m5C)
    if psudepath:
        addData2(knownPos, psudepath, Flg_Y)
    if editingpath:
        addData2(knownPos, editingpath, Flg_I)

    return knownPos