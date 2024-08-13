from heapq import heappop, heappush
from itertools import combinations, permutations
import sys
import numpy as np

sys.path.append("/Users/danielum/Documents/MSCS/Spring_2022/COMS_4762_COMS_E6901/VG/VG")
from library.reader import Reader
from PineConeExperiments.kmeans import Cluster

# pip install -r requirements.txt


class Mapper:
    def __init__(self):
        pass

    # function phi() maps nucleotides to integers A:0, C:1, G:2,T:3
    def phi(self, base):
        if base == "A":
            return 0
        if base == "C":
            return 1
        if base == "G":
            return 2
        if base == "T":
            return 3

    # each k-mer (subsequence s with lenght k) in r_{ij} can be mapped to a 2xk bit integer
    def encoding_function(self, s):
        ret = 0
        k = len(s)
        for a in range(k):
            ret += self.phi(s[a]) * 4 ** (a - 1)
        return ret

    def minimizer(self, s):
        ret = 0
        k = len(s)
        for a in range(1, k + 1):
            ret += s[a]
        ret *= 4 ** (a - 1)
        return self.phi(ret)

    def get_kmer_frequency(self, kmer, read_list):
        kmer = "".join(kmer)
        return read_list.count(kmer)

    # for each reads set R_i, the corresponding feature vector V_i = {v_i1,v_i2,...v_i4k}
    def feature_vector(self, reads_matrix, minimizer_list_of_lists):
        feature_matrix = [[] for _ in range(len(reads_matrix))]
        for i, read_list in enumerate(reads_matrix):
            n = len(read_list)
            for _, kmer in enumerate(minimizer_list_of_lists[i]):

                feature_matrix[i].append(self.get_kmer_frequency(kmer, read_list) / n)
        # return value is a i-len(reads_matrix) rows by j-len(minimizer_list) columns
        # i.e. list of vectors, one for each R_i read set
        return feature_matrix

    def select_m_highest_variance(self, feature_matrix, m):
        # m features are selected from the original 4^k features (where k is length of subsequence s)
        # a set of m-dimensional n vectors are generated, denoted by V = {V_1,...V_n}
        # what are the features?
        # transpose? YES! You need to remove rows with the lowest variance

        # transpose
        feature_matrix_transposed = np.array(feature_matrix).T.tolist()

        # list(map(list, zip(*feature_matrix)))
        # list(map(list, itertools.zip_longest(*l, fillvalue=None)))

        m_heap = []
        variance_list = []
        for feature_list in feature_matrix_transposed:
            avg = sum(feature_list) / len(feature_list)
            var = sum((x - avg) ** 2 for x in feature_list) / len(feature_list)
            variance_list.append(var)

            heappush(m_heap, (-1 * var, feature_list))
        # m-dimensional n vectors
        ret = []
        for _ in range(m):
            tmp = heappop(m_heap)
            ret.append(tmp[1])

        ret = np.array(ret).T.tolist()
        return ret

    def make_kmers(self, sequences, k):

        minimizer_list_of_lists = []

        for sequence in sequences:
            minimizer_list = combinations(sequence, k)
            minimizer_list = list(set(minimizer_list))
            minimizer_list_of_lists.append(minimizer_list)

        return minimizer_list_of_lists

def main():
    R = Reader()
    mp, count, total_len = R.read_fasta("sample/pathogen.fa")
    sequences_dict_items = mp.values()
    sequences = list(sequences_dict_items)

    M = Mapper()
    minimizer_list_of_lists = M.make_kmers(sequences, 2)
    return
    # ret = []
    # for kmer in minimizer_list:
    #     # print("kmer" + str(kmer))
    #     ret.append(M.encoding_function(kmer))

    # print(ret)
    # print("len(sequences): " + str(len(sequences)))
    # print("len(minimizer): " + str(len(minimizer_list)))

    feature_matrix = M.feature_vector(sequences, minimizer_list_of_lists)
    # print(np.shape(feature_matrix))
    # print(feature_matrix[0])

    ret = M.select_m_highest_variance(feature_matrix, 3)
    # print("number of sequences, number of kmer minimizers:")
    # print(np.shape(ret))
    # print("number of sequences, kmers columns with m most variance:")
    # print(ret[0])

    C = Cluster()
    cluster_set = C.k_means(ret, 0.5)
    print(cluster_set)

if __name__ == '__main__':
    main()