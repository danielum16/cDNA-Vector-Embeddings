from heapq import heappop, heappush
from itertools import combinations
import numpy as np

from library.reader import Reader


class Mapper:
    def __init__(self):
        pass

    def phi(self, base):
        if base == "A":
            return 0
        if base == "C":
            return 1
        if base == "G":
            return 2
        if base == "T":
            return 3

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

    def feature_vector(self, reads_matrix, minimizer_list_of_lists):
        feature_matrix = [[] for _ in range(len(reads_matrix))]
        for i, read_list in enumerate(reads_matrix):
            n = len(read_list)
            for _, kmer in enumerate(minimizer_list_of_lists[i]):

                feature_matrix[i].append(self.get_kmer_frequency(kmer, read_list) / n)
        return feature_matrix

    def select_m_highest_variance(self, feature_matrix, m):
        feature_matrix_transposed = np.array(feature_matrix).T.tolist()

        m_heap = []
        variance_list = []
        for feature_list in feature_matrix_transposed:
            avg = sum(feature_list) / len(feature_list)
            var = sum((x - avg) ** 2 for x in feature_list) / len(feature_list)
            variance_list.append(var)

            heappush(m_heap, (-1 * var, feature_list))
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


if __name__ == '__main__':
    main()