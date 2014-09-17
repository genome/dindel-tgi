


def HomopolymerLength(seq = [], pos = 0):
    hp_len = 1
    for i in range(pos+1, len(seq)):
        if seq[i]==seq[i-1]:
            hp_len += 1
        else:
            break

    for i in range(pos-1, 0,-1):
        if seq[i]==seq[i+1]:
            hp_len += 1
        else:
            break


    return hp_len
