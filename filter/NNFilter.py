from nnmodel.NNModel import *
from MSUtils import *
import pysam
def applyNNFilter(datalist,ref,checkpoint_path_A , checkpoint_path_C, checkpoint_path_T):

    fasta = pysam.FastaFile(ref)
    chroms = set()
    sequence_dic = {}
    for columns in datalist:

        chrom = columns[0]
        strand = columns[5]
        pos = int(columns[1])
        alt = columns[3]


        start = pos - 20
        end = pos + 21
        if "random" in chrom or "alt" in chrom or  "M" in chrom:
            continue

        try:

            sequence = fasta.fetch(chrom, start, end).upper()

        except ValueError as e:

            print(f"[Warning] Failed to fetch {chrom}:{start}-{end} ({e})")
            sequence = "A"*41

        if strand == "-":
            sequence = reverse_complement(sequence)

        alt2key = {"a": "A", "17596": "A", "m": "C"}
        key = alt2key.get(alt, "T")

        sequence_list = sequence_dic.get(key,[])
        sequence_dic[key] = sequence_list

        if len(sequence) != 41:
            dummyseq = "A" * 41
            sequence = dummyseq

        sequence_list.append((columns,sequence))

    retlist = []
    for key in sequence_dic:

        if key == "A":
            checkpoint_path = checkpoint_path_A
            numclass = 4
        elif key == "C":
            checkpoint_path = checkpoint_path_C
            numclass = 3
        else:
            checkpoint_path = checkpoint_path_T
            numclass = 3

        model = getModel(numclass)
        # model.summary()
        # model.compile()
        # print(sequence_dic.keys())
        # print("checkpoint_path",checkpoint_path_A)
        # print("checkpoint_path", checkpoint_path_C)
        # print("checkpoint_path", checkpoint_path_T)
        model.load_weights(checkpoint_path)

        sequence_list = sequence_dic.get(key, [])
        numlist = toNumberList2(sequence_list)

        X_test = []
        for item in numlist:
            if (item is not None) and (len(item[0])==39):
                X_test.append(item[0])
            else:
                dummyseq = "A" * 41
                s = encode_dna(dummyseq)
                X_test.append(s)


        X_test = np.array(X_test)
        print("X_len",len(X_test))

        # model = tf.keras.models.load_model(checkpoint_path)
        predict = model.predict(X_test)
        idx = 0
        for columns,sequence in sequence_list:

            pre = np.argmax(predict[idx])
            if key == "C" and pre == 2:
                pre = 4
            if key == "T" and pre == 2:
                pre = 5

            retlist.append((columns,sequence,pre))
            chroms.add(columns[0])
            idx+=1

    return retlist


# Flg_Error = 0
# Flg_other = 1
# Flg_m6A = 2
# Flg_I = 3
# Flg_m5C = 4
# Flg_Y = 5