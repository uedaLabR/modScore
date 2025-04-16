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

# Load BED file to create a filter dictionary
def load_filter_bed(filter_bed):
    alt_dict = {}
    with open(filter_bed, "r") as bed:
        for line in bed:
            columns = line.strip().split("\t")
            chrom, pos, alt, passfail = columns[0], int(columns[1]), str(columns[3]), columns[-1]
            if passfail == "Pass":
                alt_dict.setdefault(alt, {}).setdefault(chrom, set()).add(pos)
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


import pysam
import copy
import numpy as np
def run_recalib(inbam, outbam, filter_bed):

    print("start bem recalib")

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

                refposs = read.get_reference_positions()
                modbase = read.modified_bases
                ML = read.get_tag("ML")
                orginal_array = copy.copy(ML)

                if modbase is not None:
                    modkeys = list(modbase.keys())
                    mm_tag = read.get_tag("MM")
                    modkeys = sortbyMMTagKeyInfo(modkeys, mm_tag)
                    mlindex = 0
                    for modkey in modkeys:

                        filter_set = all_dict.get(str(modkey[2]), {}).get(chrom, None)
                        modlist = modbase[modkey]
                        processed_tuples = [(convertToGenomepos(x, refposs), x, y) for x, y in modlist]
                        for tp in processed_tuples:

                            if mlindex >= len(ML):
                                break

                            pos, localpos, originalscore = tp
                            #
                            counter = cntdict.get(str(modkey[2]), (0,0))
                            clear,kept = counter

                            if filter_set is not None and (
                                    (pos - 1 in filter_set) or (pos in filter_set) or (pos + 1 in filter_set)):
                                kept += 1
                            else:
                                ML[mlindex] = 0
                                clear += 1
                            #
                            cntdict[str(modkey[2])] = (clear,kept)

                            mlindex += 1

                read.set_tag("ML", ML)
                read.set_tag("XM", orginal_array)

            readcnt += 1
            if readcnt % 50000 == 0:
                print(readcnt,cntdict)
                # break

            # Write the modified read to the output BAM file
            bam_out.write(read)


