
import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import to_categorical
import pandas as pd
from NNModel import *
from MSUtils import *
import tensorflow as tf

def train(data,weightpath,epoch,outhistory):

    # structlist = calc2ndaryStruct(data,calcstruct=False)
    # X = [item[0] for item in structlist]
    # y = [item[1] for item in structlist]

    numlist = toNumberList(data)
    X = []
    y =[]
    for item in numlist:
        if item is not None:
            s = item[0]
            f = item[1]
            if len(s) ==39:
                X.append(s)
                y.append(f)

    X = np.array(X)
    y = np.array(y)
    print(X.shape)
    print(y.shape)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.15, random_state=42)


    y_train = to_categorical(y_train, num_classes=6)
    y_test = to_categorical(y_test, num_classes=6)
    X_train = np.array(X_train)
    y_train = np.array(y_train)


    print("X_train shape:", X_train.shape)
    print("y_train shape:", y_train.shape)
    print("X_test shape:", X_test.shape)
    print("y_test shape:", y_test.shape)

    print("X_train shape:", X_train.shape)
    print("y_train shape:", y_train.shape)
    print("X_train sample:", X_train[:5])
    print("y_train sample:", y_train[:5])

    model = getModel()
    model.summary()

    checkpoint = ModelCheckpoint(filepath=weightpath,
                                 save_weights_only=True,
                                 save_best_only=True,
                                 monitor='val_acc',
                                 verbose=1)

    # Preparing callbacks.
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=0.001,  # Starting learning rate
        decay_steps=10000,  # After how many steps to apply decay
        decay_rate=0.9,  # The decay rate to apply
        staircase=True)  # If True, decay the learning rate at discrete intervals

    # Initialize the Adam optimizer with the learning rate schedule
    adam = tf.keras.optimizers.Adam(
        learning_rate=lr_schedule,  # Use the learning rate schedule
        beta_1=0.9,
        beta_2=0.999,
        epsilon=1e-07,  # For default value, you can simply omit this or set it to 1e-7
        amsgrad=False)

    model.compile(optimizer=adam,
                  # loss='sparse_categorical_crossentropy',
                  loss='categorical_crossentropy',
                  metrics=['acc'])


    callbacks = [
        checkpoint
    ]

    # seed_everything(42)

    print(X_train)
    # Train the model.
    history = model.fit(x=X_train,
                        y= y_train,
                        batch_size=64,
                        epochs=epoch,
                        callbacks=callbacks,
                        shuffle=True,
                        validation_data=(X_test,y_test))

    history_df = pd.DataFrame(history.history)
    history_df.to_csv(outhistory, index=False)

    output_csv_path = "/share/ueda/nanoModiTune/classval.csv"
    evaluate_validation_set(model,X_test, y_test, output_csv_path)

# data = []
# data.append(("GCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATT",1))
# data.append(("TTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCT",1))

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
import random




import pysam
import random
import re


def fetch_random_sequences(fasta, sequence_length=41, number_of_sequences=5):


    autosomes = [ref for ref in fasta.references if re.match(r'chr\d+$', ref)]
    selected_chromosome = random.choice(autosomes)

    chromosome_length = fasta.get_reference_length(selected_chromosome)

    sequences = []
    for _ in range(number_of_sequences*10):
        start = random.randint(0, chromosome_length - sequence_length)
        end = start + sequence_length
        sequence = fasta.fetch(selected_chromosome, start, end).upper()
        if "N" in sequence:
            continue
        if len(sequences) ==  number_of_sequences:
            break
        if sequence[20] == "G":
            continue
        # print(sequence,len(sequence))
        sequences.append((sequence,Flg_other))

    fasta.close()

    return sequences



def is_drach(sequence):

    # Check if the sequence length is exactly 5
    if len(sequence) != 5:
        return False

    # Define the allowed bases for each position in the DRACH motif
    d_bases = ['A', 'G', 'T']
    r_bases = ['A', 'G']
    a_base = ['A']
    c_base = ['C']
    h_bases = ['A', 'C', 'T']

    # Check each position against the allowed bases
    if sequence[0] in d_bases and sequence[1] in r_bases and sequence[2] in a_base and sequence[3] in c_base and \
            sequence[4] in h_bases:
        return True
    else:
        return False




def addData(data,file,flg,nuc,maxcnt):

    bed_df = pd.read_csv(file, header=None, sep='\t')
    column_19 = bed_df.iloc[:, 18]
    cnt = 0
    print(len(column_19))
    for x in column_19:
        # pseq = x[1:]
        pseq = x
        print(pseq,len(pseq),pseq[20])
        if pseq[20] == nuc:
            pseq = reverse_complement(pseq)
        # pseq = pseq[0:40]
        if len(pseq) == 41:
            data.append((pseq, flg))
            cnt += 1
            # if cnt > maxcnt:
            #     break
    print("append",cnt)

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}  # Ns
    return ''.join(complement.get(base, base) for base in reversed(seq))


import pandas as pd
import pysam
def getData(data,m6Apath,m5Cpath, psudepath,editingpath, fp_ivtpath,fasta_path,fasta_path_hg38):

    maxcnt = 40000

    scorethres = 15
    ratiothres = 10
    sequences = []
    fasta = pysam.FastaFile(fasta_path_hg38)

    with open(fp_ivtpath, "r") as bed_file:
        for line in bed_file:
            columns = line.strip().split("\t")
            alt = columns[3]
            score = int(columns[4])
            ratio = float(columns[10])
            # Check if the entry meets the filtering criteria.
            if score >= scorethres and ratio > ratiothres:

                chrom = columns[0]
                pos = int(columns[1])
                alt = columns[3]
                strand = columns[5]
                # if strand == "-":
                #     pos = pos
                # else:
                #     pos = pos +1

                start = pos - 20
                end = pos + 21
                if "M" in chrom:
                    break

                sequence = fasta.fetch(chrom, start, end).upper()
                if strand == "-":
                    sequence = reverse_complement(sequence)

                #
                # print(line)
                print(sequence,len(sequence),strand,sequence[20],alt)
                sequences.append((sequence, Flg_Error))


    data.extend(sequences)
    print("fp data size=",len(sequences))
    fasta = pysam.FastaFile(fasta_path)
    tuple_list3 = fetch_random_sequences(fasta, sequence_length=41, number_of_sequences=(maxcnt*3))
    data.extend(tuple_list3)
    print("randomseq=", len(tuple_list3))
    print("maxcnt",maxcnt)
    addData(data, editingpath, Flg_I, 'T',maxcnt-1)
    addData(data, m6Apath, Flg_m6A, 'T', maxcnt-1)
    addData(data, m5Cpath, Flg_m5C, 'G', maxcnt-1)
    addData(data, psudepath, Flg_Y, 'A', maxcnt-1)

    #
    return data


Flg_Error = 0
Flg_other = 1
Flg_m6A = 2
Flg_I = 3
Flg_m5C = 4
Flg_Y = 5

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

import glob
def getFiles(sourcepath,genome):

    m6Apath, m5Cpath, psudepath, editingpath = "","","",""
    files = glob.glob(sourcepath+"/db1/*")
    # print("files",files,sourcepath+"/*")
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

def addData2(data,ref,path,flg,adjust=-1):

    posset = set()
    print(path)
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
                chrom,startend = chrom_startend
                #
                if "-" in startend:
                    print(startend)
                    startend = startend.split("-")
                    key1 = chrom+":"+ str(int(startend[0])+adjust)
                    posset.add(key1)
                    key2 = chrom+":"+ str(int(startend[1])+adjust)
                    posset.add(key2)
                else:
                    key1 = chrom + ":" + str(int(startend)+adjust)
                    posset.add(key1)

            else:

                cols= line.split("\t")
                key1 = cols[0] + ":" + str(int(cols[1])+adjust)
                posset.add(key1)

        fasta = pysam.FastaFile(ref)
        cnt = 0
        for poskey in posset:

            chrom,pos = poskey.split(":")
            pos = int(pos)

            start = pos - 20
            end = pos + 21
            if "M" in chrom:
                break
            if start < 0:
                continue

            sequence = fasta.fetch(chrom, start, end).upper()

            if len(sequence) < 35:
                continue

            nuc = sequence[20]
            print(sequence,len(sequence),nuc)
            strand = "+"
            if flg == Flg_m5C:
                if nuc == "G":
                    strand = "-"
            if flg == Flg_m6A or flg == Flg_I:
                if nuc == "T":
                    strand = "-"
            if flg == Flg_Y:
                if nuc == "A":
                    strand = "-"

            if strand == "-":
                sequence = reverse_complement(sequence)


            if len(sequence) != 41:
                continue
            print(sequence,len(sequence),sequence[20])
            #
            data.append((sequence, flg))
            cnt+=1
            # if cnt > 20:
            #     break



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

import random
from collections import defaultdict
def sample_by_flag(data, max_per_flag=50000):

    groups = defaultdict(list)
    for item in data:

        groups[item[1]].append(item)

    sampled_data = []

    for flag, items in groups.items():
        if len(items) > max_per_flag and (flag >=2):
            # 50,000_o
            sampled = random.sample(items, max_per_flag)
        else:
            # 50,000p
            sampled = items
        sampled_data.extend(sampled)

    return sampled_data

import numpy as np
from sklearn.metrics import confusion_matrix
import pandas as pd

def evaluate_validation_set(model,X_test, y_test, output_csv_path):
    """
    Evaluates the model on the validation set and outputs the confusion matrix.
    :param model: Trained model
    :param X_test: Validation dataset features
    :param y_test: Validation dataset true labels (one-hot encoded)
    :param class_names: List of class names
    :param output_csv_path: Path to save the confusion matrix as CSV
    :return: None
    """
    # Perform inference
    predictions = model.predict(X_test)

    # Convert one-hot encoded y_test to class indices
    true_labels = np.argmax(y_test, axis=1)
    #

    # Get predicted class indices
    predicted_labels = np.argmax(predictions, axis=1)

    # Calculate confusion matrix
    conf_matrix = confusion_matrix(true_labels, predicted_labels)

    # Display confusion matrix
    print("Confusion Matrix:")
    print(conf_matrix)

    # Optionally save the confusion matrix to a CSV
    if output_csv_path:
        conf_matrix_df = pd.DataFrame(conf_matrix)
        conf_matrix_df.to_csv(output_csv_path)
        print(f"Confusion matrix saved to {output_csv_path}")

    # Return confusion matrix for further analysis if needed
    return conf_matrix

import glob
def getRef(sourcepath, genome):

    p = sourcepath+"/genome/*"
    files = glob.glob(p)
    for file in files:
        if genome in file:
            return file
    return ""

def trainNN(sourcepath, genome, fp_ivtpath, outhistory, weightpath,eachsize=50000, epoch=100):

    ref = getRef(sourcepath, genome)
    refhg38 = getRef(sourcepath, "hg38")
    m6Apath, m5Cpath, psudepath, editingpath = getFiles2(sourcepath,genome)
    print((m6Apath, m5Cpath, psudepath, editingpath))
    data = []
    if m6Apath:
        addData2(data,ref,m6Apath,Flg_m6A)
    if m5Cpath:
        addData2(data,ref,m5Cpath,Flg_m5C)
    if psudepath:
        addData2(data,ref,psudepath,Flg_Y)
    if editingpath:
        addData2(data,ref,editingpath,Flg_I)
    # #
    m6Apath, m5Cpath, psudepath, editingpath = getFiles(sourcepath,genome)
    print(m6Apath, m5Cpath, psudepath, editingpath)
    print("----")
    data = getData(data,m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,refhg38)
    print("finish get data")

    print("random sampling")
    data = sample_by_flag(data)
    print("end random sampling")
    #
    flg_labels = {
        0: "Error",
        1: "Other",
        2: "m6A",
        3: "Inosine",
        4: "m5C",
        5: "Y"
    }


    from collections import Counter

    flg_counts = Counter()

    for sequence, flg in data:
        flg_counts[flg] += 1

    for flg, count in flg_counts.items():
        label = flg_labels.get(flg, "Unknown")
        print(f"{label} (Flag {flg}): {count} sequences")

    # weightpath = sourcepath+"/model_weight/"+genome+".weights.h5"
    train(data,weightpath,epoch,outhistory)




