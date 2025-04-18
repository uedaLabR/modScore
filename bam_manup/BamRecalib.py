from numba import njit
from numba import njit, int64


def toInt(s):
    try:
        return int(s)
    except ValueError:
        return s



def sortbyMMTagKeyInfo(modkeys,mm_tag):

    order_mapping ={}
    modifications = mm_tag.split(";")
    n=0
    for mod in modifications:
        mod = mod.split(",")[0].replace(".","")
        if mod:
            if '+' in mod:
                base_mod_strand, modkey = mod.split("+", 1)
            elif '-' in mod:
                base_mod_strand, modkey = mod.split("-", 1)
            else:
                continue
            order_mapping[toInt(modkey)] = n
            n+=1

    if len(order_mapping)==0:
        #MM tag parse err, shold not happen
        order_mapping = {17596: 0, 'a': 1, 'm': 2, 17802: 3}

    sorted_tuples = sorted(modkeys, key=lambda x: order_mapping.get(x[2], float('inf')))
    return sorted_tuples

def convertToGenomepos(x,refposs):

    if x < 0:
        return 0
    if x >= len(refposs):
        return 0
    conv =  refposs[x]
    if conv is not None:
        conv = conv+1
    return conv

@njit
def convert_to_genomepos_nb(x, refposs):
    """
    refposs: int64[:], None  -1 \
    return: 1-based u or 0
    """
    if x < 0 or x >= refposs.shape[0]:
        return 0
    rp = refposs[x]
    if rp >= 0:
        return rp + 1
    else:
        return 0


# Load BED file to create a filter dictionary
def load_filter_bed(filter_bed):
    alt_dict = {}
    with open(filter_bed) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0]
            pos = int(cols[1])
            alt = cols[3]
            passfail = cols[-1]
            if passfail == "Pass":
                # K typed.Set
                alt_dict.setdefault(alt, {}) \
                        .setdefault(chrom, set()) \
                        .add(pos)
    # A alt/chrom  Set
    return alt_dict

# Check if position is filtered
def matchedPoss(chrom, pos, modkey,all_dict):

    debug = False
    # if modkey == "17596":
    #     print(modkey,all_dict.keys())
    #     debug = True

    altdict = all_dict.get(modkey,None)
    if debug:
        print(altdict.keys())
    if altdict:
        chrom_posdict = altdict.get(chrom,None)
        if debug:
            print(chrom,len(chrom_posdict),chrom_posdict)
        if chrom_posdict:
            # print(pos)
            if (pos-1 in chrom_posdict) or (pos in chrom_posdict) or (pos+1 in chrom_posdict):

                return True

    return False

@njit
def updateML(ML, mlindex, positions, filter_set):

    clear = 0
    kept  = 0
    n_ml  = ML.shape[0]
    for i in range(positions.shape[0]):
        if mlindex >= n_ml:
            break
        p = positions[i]
        # }1
        if (p   in filter_set) or \
           (p-1 in filter_set) or \
           (p+1 in filter_set):
            kept += 1
        else:
            ML[mlindex] = 0
            clear += 1
        mlindex += 1
    return mlindex, clear, kept




import pysam
import copy
import numpy as np
def run_recalib(inbam, outbam, filter_bed):

    print("start bam recalib")
    empty = set()

    all_dict = load_filter_bed(filter_bed)
    with pysam.AlignmentFile(inbam, "rb") as bam_in, pysam.AlignmentFile(outbam, "wb", template=bam_in) as bam_out:

        readcnt = 0
        cntdict = {}
        for read in bam_in:

            if not read.is_mapped:
                bam_out.write(read)
                continue

            chrom = bam_in.get_reference_name(read.reference_id)

            if read.has_tag("MM") and read.has_tag("ML"):

                refposs_np = np.array([p if p is not None else -1 for p in read.get_reference_positions()],
                                      dtype=np.int64)


                modbase = read.modified_bases
                ML = read.get_tag("ML")
                ML_arr = np.array(ML, dtype=np.int32)

                if modbase is not None:

                    modkeys = list(modbase.keys())
                    mm_tag = read.get_tag("MM")
                    modkeys = sortbyMMTagKeyInfo(modkeys, mm_tag)
                    mlindex = 0

                    for modkey in modkeys:

                        key_str = str(modkey[2])
                        clear, kept = cntdict.get(key_str, (0, 0))
                        filter_set = all_dict.get(key_str, {}).get(chrom, empty)

                        modlist = modbase[modkey]
                        positions = np.array([convert_to_genomepos_nb(x, refposs_np) for x, y in modlist],
                                             dtype=np.int64)
                        mlindex,_clear,_kept = updateML(ML_arr,mlindex,positions,filter_set)
                        cntdict[key_str] = (clear+_clear, kept+_kept)

                ML = ML_arr.tolist()
                read.set_tag("ML", ML)


            readcnt += 1
            if readcnt % 50000 == 0:
                print(readcnt,cntdict)
                # break

            # Write the modified read to the output BAM file
            bam_out.write(read)


