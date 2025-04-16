from filter.FilterMod import *
from filter.LoadPos import *

def filterEachMod(datalist_filter,knownPos,knownDB,repeat,genome):

    alt_dict = {}
    for columns,sequence,pre in datalist_filter:
        alt = columns[3]
        score = int(columns[4])
        ratio = float(columns[10])
        # Create a new inner dictionary for this alt if it doesn't exist.
        if alt not in alt_dict:
            alt_dict[alt] = {}
        chrom = columns[0]
        # Within this alt group, create a list for the chromosome if needed.
        if chrom not in alt_dict[alt]:
            alt_dict[alt][chrom] = []
        # Append the current entry (as a list of columns) into the corresponding chromosome group.
        alt_dict[alt][chrom].append((columns,sequence,pre))

    filterlist = []
    for alt, chrom_group in alt_dict.items():

        if alt == "m":
            # Process alt "m" using filter_m5C for all chromosomes.
            result = filter_m5C(chrom_group,  knownPos)
        elif alt == "a":
            # Process alt "a" using filter_m6A.
            result = filter_m6A(chrom_group,  knownPos)
        elif alt == "17596":
            # Process alt "17596" using filter_Inosine.
            result = filter_Inosine(chrom_group, repeat,genome,knownPos,knownDB)
        else:
            # Process any other alt using filter_PsudeUridine.
            result = filter_PsudeUridine(chrom_group, knownPos)
        # Merge the filtered results into one list.
        if result is not None:
            filterlist.extend(result)

    # Return the final filtered list for further processing if needed.
    return filterlist

from Configs import *
def bedToList(bed):

    data = []
    cnt=0
    with open(bed, "r") as bed_file:
        for line in bed_file:

            columns = line.strip().split("\t")
            alt = columns[3]
            score = int(columns[4])
            ratio = float(columns[10])
            supportread = 0.01 * ratio * score
            # Check if the entry meets the filtering criteria.


            if alt == "17596":
                depth_thres,ratio_thres,minsupport_read_thres =  ThresFor_I
            elif alt == "17802":
                depth_thres,ratio_thres,minsupport_read_thres =  ThresFor_Y
            elif alt == "m":
                depth_thres, ratio_thres, minsupport_read_thres = ThresFor_m5C
            elif alt == "a":
                depth_thres, ratio_thres, minsupport_read_thres = ThresFor_m6A

            if score >= depth_thres and ratio > ratio_thres and  supportread >= minsupport_read_thres :
                data.append(columns)
                if len(data) %10000 == 0:
                    print("load",len(data),"out of ",cnt,"scorethres",depth_thres,"ratiothres",ratio_thres)
            cnt+=1

    return data


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



import glob
def getFiles(source_path,genome):

    repeat, ref, checkpoint_path = "","",""
    dirs = ["genome","repeat","model_weight"]
    for s in dirs:

        p = source_path+"/"+s+ "/*"
        if s == "genome":
            p = source_path + "/" + s + "/*.fa"
        files = glob.glob(p)
        for file in files:

            if genome in file:

                if s == "genome":
                    ref = file
                if s == "repeat":
                    repeat = file
                if s == "model_weight":
                    checkpoint_path = file

    return repeat, ref, checkpoint_path


from filter.NNFilter import *
def filterBed(bed, bed_out, source_path, genome):

    print("start filterbe")

    print("loading known pos")
    knownPos = loadKnownPos(source_path, genome)
    print("loading known pos2")
    addOtherDB(knownPos, source_path,genome)

    repeat, ref, checkpoint_path = getFiles(source_path,genome)

    print("loading bed")
    datalist = bedToList(bed)
    print("apply sequeunce filter")
    datalist_filter = applyNNFilter(datalist,ref,checkpoint_path)


    print("filter each modification")
    filterList = filterEachMod(datalist_filter,knownPos,knownDB,repeat,genome)
    output(bed_out, filterList)

def output(bed_out, filterList):

    with open(bed_out, "w") as out_file:
        for columns in filterList:

            line = "\t".join(columns)
            if not line.endswith("\n"):
                line += "\n"
            out_file.write(line)



