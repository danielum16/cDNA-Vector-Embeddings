from library import reader
import numpy as np
import warnings
import gzip

import math
from encoder import Encoder


class Mapper:
    __mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    __rmapping = {0: "A", 1: "C", 2: "G", 3: "T"}

    def __init__(self, R: list[list], k, m):
        self.fs = Mapper.feature_set(R, k)
        self.hfs = Mapper.select_high_variance(self.fs, m)

    @classmethod
    def encode(cls, s: list | str) -> int:
        res = 0
        for i, c in enumerate(s):
            res += cls.__mapping[c] * (4**i)
        return res

    @classmethod
    def decode(cls, c: int, k: int) -> str:
        res = ""
        for i in range(k):
            res += cls.__rmapping[c % 4]
            c //= 4
        return res

    @classmethod
    def vec(cls, i: int, j: int, k: int, R: list[list | str]) -> float:
        return (R[i].count(cls.decode(j, k)) / len(R[i]), i)

    @classmethod
    def feature_set(cls, R: list[list | str], k: int):
        return [[cls.vec(i, j, k, R) for j in range(4**k)] for i in range(len(R))]

    @classmethod
    def select_high_variance(cls, feature_set: list[list], m: int) -> list[list]:
        fm = np.array(feature_set)
        # [:,:,0] = [...,0]
        var = np.var(fm[:, :, 0], 0)
        return fm[:, var >= np.sort(var)[-m]]


def kmeans(V, k):
    from sklearn.cluster import KMeans

    res = KMeans(n_clusters=k, random_state=0).fit(V[..., 0])
    return V[res.labels_ == 0], V[res.labels_ == 1]


def dist(v1: np.ndarray, v2: np.ndarray) -> float:
    return np.linalg.norm(v1 - v2)


def mindist(ca0: np.ndarray, cb: np.ndarray) -> np.ndarray:
    d = np.array([dist(cbi, ca0) for cbi in cb])
    return cb[d == np.min(d)][0]


def cluster(V, l=0.3):
    if len(V) < 3:
        return [V]
    else:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c1, c2 = kmeans(V, 2)
        if len(c1) == 0:
            return [c2]
        if len(c2) == 0:
            return [c1]

        if len(c1) + len(c2) >= 4 and (len(c1) == 1 or len(c2) == 1):
            if len(c1) != 1:
                c2, c1 = c1, c2
            ca, cb = c1, c2

            cbmin = mindist(ca[0], cb)
            if dist(ca[0], cbmin) < l * dist(mindist(cbmin, cb[cb != cbmin])):
                ca = np.array([ca[0], cbmin])
                cb = cb[cb != cbmin]
                return cluster(cb, l) + [ca]
        else:
            return cluster(c1, l) + cluster(c2, l)


def transform(V):
    S = []
    for v in V:
        S.append((v[..., 0], int(v[0, 1])))
    return S


def encode_base(base):
    if base == "A":
        # 1000
        return 8
    elif base == "C":
        # 0100
        return 4
    elif base == "G":
        # 0010
        return 2
    elif base == "T":
        # 0001
        return 1
    else:
        # 1111
        return 15


def encode_sequence(sequence):
    val = 0
    for i, base in enumerate(sequence):
        val += encode_base(base) * (16**i)
    return val


def groupings(S, sequences):
    clust = []
    # S = cluster(mapper.hfs)
    for x, V in enumerate(S):
        sub_list = []
        for _, i in transform(V):
            sub_list.append((sequences[i], i))
        clust.append(sub_list)
    return clust


def compress_clusters(int_binary_cluster):
    """
    The problem we are facing is that a byte can only represent numbers from 0-255.
    What we can do is split an integer N into multiple smaller integers {a0, a1, a2,..., a_k}
    where N = Σ (a_j * (256^k), then convert that list of smaller integers into a byte array.
    """

    comp = []
    for int_sequences in int_binary_cluster:
        comp.append(gzip.compress(int_sequences))
    return comp


def get_size_compressed(compressed_cluster):
    return sum([len(cluster) for cluster in compressed_cluster])


def decode_base(int_base):
    if int_base == 8:
        # 1000
        return "A"
    elif int_base == 4:
        # 0100
        return "C"
    elif int_base == 2:
        # 0010
        return "G"
    elif int_base == 1:
        # 0001
        return "T"
    else:
        # 1111
        return "N"


def decode_sequence(val):
    sequence = ""
    n = math.floor(math.log(val) / math.log(16))

    while val > 0:
        next_layer = val % 16**n
        sequence = str(decode_base(int((val - next_layer) / 16**n))) + sequence
        n -= 1
        val = next_layer

    return sequence


def decode_cluster(output_cluster, output_indexes):
    cluster_of_sequences = []

    for i, encoded_sequences in enumerate(output_cluster):
        inner_list = []
        for j, encoded_sequence in enumerate(encoded_sequences):
            inner_list.append((decode_sequence(encoded_sequence), output_indexes[i][j]))
        cluster_of_sequences.append(inner_list)

    return cluster_of_sequences


def convert(num: int) -> list:
    numlist = []
    while num > 0:
        base = num % 16
        match base:
            case 1:
                numlist.append([0, 0, 0, 1])
            case 2:
                numlist.append([0, 0, 1, 0])
            case 4:
                numlist.append([0, 1, 0, 0])
            case 8:
                numlist.append([1, 0, 0, 0])
            case _:
                numlist.append([1, 1, 1, 1])
        num //= 16
    return numlist


if __name__ == "__main__":
    R = reader.Reader()
    mp, count, total_len = R.read_fasta("sample/sandmouse.fa")
    sequences_dict_items = mp.values()
    sequences = list(sequences_dict_items)

    mapper = Mapper(sequences, 2, 3)
    groups_of_similar_kmers = cluster(mapper.hfs)

    cluster_of_sequences = groupings(groups_of_similar_kmers, sequences)

    c = Encoder.encode_clusters(cluster_of_sequences)
    compressed_cluster = compress_clusters(c)
    print(compressed_cluster)
    size = get_size_compressed(compressed_cluster)
    print(size)
