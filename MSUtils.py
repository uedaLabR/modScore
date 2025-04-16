Flg_Error = 0
Flg_other = 1
Flg_m6A = 2
Flg_I = 3
Flg_m5C = 4
Flg_Y = 5

def toNumberList(data):

    ret = []
    cnt = 0
    for seq,label in data:
        tokens = tokonizeN(seq)
        ret.append((tokens, label))
        cnt += 1

    return ret

def tokonizeN(seq):

    ret = []
    for n in range(len(seq)-2):
        ps = seq[n:n+3]
        token = encode_dna(ps)
        ret.append(token)

    return ret

def encode_dna(dna):

    base_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3,'N':0}
    # Convert DNA sequence and secondary structure to numerical values
    dna_value = sum(base_mapping[base] * (4 ** i) for i, base in enumerate(reversed(dna)))
    return  dna_value

def reverse_complement(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N': 'N'}
    reverse_comp = ''.join([complement[base] for base in reversed(dna_seq)])
    return reverse_comp
