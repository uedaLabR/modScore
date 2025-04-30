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
def convert_to_genomepos_nb(x, qstart, refposs):
    """
    refposs: int64[:], None  -1 \
    return: 1-based u or 0
    """
    x = x-qstart
    if x < 0 or x >= refposs.shape[0]:
        return 0
    rp = refposs[x]
    if rp >= 0:
        return rp
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

# @njit
def updateML(ML, mlindex, positions, filter_set,drachSet):

    clear = 0
    kept  = 0
    n_ml  = ML.shape[0]
    for i in range(positions.shape[0]):

        if mlindex >= n_ml:
            break
        if ML[mlindex] == 0:
            mlindex += 1
            continue

        p = positions[i]
        inset = (p in filter_set) or (p - 1 in filter_set) or (p + 1 in filter_set)
        drach = False
        if drachSet:
            drach = (p in drachSet)
        # print("pos in read",p,inset)
        if inset or drach:
            kept += 1
            # print("pos in read", p,mlindex,ML[mlindex],i, inset)
        else:
            ML[mlindex] = 0
            clear += 1
        mlindex += 1

    return mlindex, clear, kept



import pysam
from typing import Optional, List


def build_read_to_ref_map(read: pysam.AlignedSegment) -> List[Optional[int]]:

    qlen = read.query_length
    mapping: List[Optional[int]] = [None] * qlen
    query_consumed = 0
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            for i in range(length):
                mapping[query_consumed] = ref_pos
                query_consumed += 1
                ref_pos       += 1

        elif op in (2, 3):
            # deletion / ref-skipFref i
            ref_pos += length

        elif op == 1:
            # insertionFquery ii} None j
            query_consumed += length

        elif op == 4:
            # soft clipFquery i
            # query_consumed += length
            pass

        # 5=H,6=P:

    return mapping


def rev_comp(seq: str) -> str:
    """Return the reverse complement of the given DNA sequence."""
    comp_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp_table)[::-1]

import re
from typing import List, Tuple, Union, Set
def checkDRACH(seq_plus,nonstrand,positions) -> Set[int]:

    DRACH_PATTERN = re.compile(r'^[AGT][AG]AC[ACT]$')
    hits: Set[int] = set()

    for pos in positions:

        start = pos - 2
        end = pos + 3  # fetch end is exclusive

        if start < 0:
            continue

        try:
            # seq = fasta.fetch(chrom, start, end).upper()
            seq = seq_plus[start:end]
        except ValueError:
            # out-of-range fetch
            continue

        if nonstrand:
            seq = rev_comp(seq)

        # print(seq)
        if len(seq) == 5 and DRACH_PATTERN.match(seq):
            hits.add(pos)

    return hits


def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

import pysam
import copy
import numpy as np
def run_recalib(inbam, outbam, filter_bed,ref):

    print("start bam recalib")
    fasta = pysam.FastaFile(ref)

    all_dict = load_filter_bed(filter_bed)
    with pysam.AlignmentFile(inbam, "rb") as bam_in, pysam.AlignmentFile(outbam, "wb", template=bam_in) as bam_out:


        all_dict = load_filter_bed(filter_bed)
        #chr6:71,587,038-71,635,009
        for chrom in bam_in.references:

            print("doing",chrom)
            contig_len = bam_in.get_reference_length(chrom)
            seq_plus = fasta.fetch(chrom, 0, contig_len).upper()
            doChrom(all_dict,chrom,seq_plus,contig_len,bam_in,bam_out)

def doChrom(all_dict,chrom,seq_plus,contig_len,bam_in,bam_out):

    empty = set()
    cntdict = {}
    readcnt = 0

    for read in bam_in.fetch(chrom, 0, contig_len):

        if not read.is_mapped:
            continue
        if not read.is_mapped:
            bam_out.write(read)
            continue

        if read.has_tag("MM") and read.has_tag("ML"):

            mp = build_read_to_ref_map(read)
            modbase = read.modified_bases
            ML = read.get_tag("ML")
            ML_arr = np.array(ML, dtype=np.int32)

            if modbase is not None:

                modkeys = list(modbase.keys())
                mm_tag = read.get_tag("MM")
                # print(mm_tag)
                # print(ML)
                modkeys = sortbyMMTagKeyInfo(modkeys, mm_tag)
                # print(modkeys)
                mlindex = 0

                for modkey in modkeys:

                    key_str = str(modkey[2])
                    nonstrand = (modkey[1] == 1)
                    clear, kept = cntdict.get(key_str, (0, 0))
                    filter_set = all_dict.get(key_str, {}).get(chrom, empty)

                    modlist = modbase[modkey]
                    positions = []
                    for x, y in modlist:
                        refp = x - read.qstart
                        if refp < 0:
                            positions.append(0)
                            continue
                        genomepos = mp[refp]
                        if genomepos:
                            positions.append(genomepos)
                        else:
                            positions.append(0)

                    if nonstrand:
                        positions.reverse()

                    #keep DRACH
                    drachSet = None
                    if key_str == "a":
                        drachSet = checkDRACH(seq_plus,nonstrand,positions)

                    # filter_array = np.fromiter(filter_set, dtype=np.int32)
                    # drach_array = np.fromiter(drachSet, dtype=np.int32) if drachSet else None
                    positions = np.array(positions)
                    mlindex,c,k = updateML(ML_arr,mlindex,positions,filter_set,drachSet)
                    cntdict[key_str] = (clear+c, kept+k)



            ML = ML_arr.tolist()
            read.set_tag("ML", ML)


        readcnt += 1
        if readcnt % 1000 == 0:
            print(chrom,read.reference_start,readcnt,cntdict)
            # break

        # Write the modified read to the output BAM file
        bam_out.write(read)



