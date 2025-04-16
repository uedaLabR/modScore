from nnmodel.NNModel import *
from MSUtils import *
import pysam
def applyNNFilter(datalist,ref,checkpoint_path):

    fasta = pysam.FastaFile(ref)
    model = getModel()
    model.summary()
    model.compile()
    model.load_weights(checkpoint_path)

    tuple_list = []
    sequence_list = []
    for columns in datalist:

        chrom = columns[0]
        pos = int(columns[1])
        strand = columns[5]
        pos = int(columns[1])
        alt = columns[3]
        strand = columns[5]

        start = pos - 20
        end = pos + 21
        if "M" in chrom:
            break

        sequence = fasta.fetch(chrom, start, end).upper()
        if strand == "-":
            sequence = reverse_complement(sequence)


        if len(sequence) == 41:
            tuple_list.append((sequence, 0))
        else:
            sequence = "A"*41
            tuple_list.append((sequence, 0))

        sequence_list.append((columns,sequence))

    numlist = toNumberList(tuple_list)

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
    retlist = []
    for columns,sequence in sequence_list:

        pre = np.argmax(predict[idx])
        retlist.append((columns,sequence,pre))
        idx+=1

    return retlist