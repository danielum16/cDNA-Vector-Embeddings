class Encoder:
    # one-hot encoding
    __mapping = {"A": 8, "C": 4, "G": 2, "T": 1, "N": 15, "E": 0}
    __rmapping = {8: "A", 4: "C", 2: "G", 1: "T", 15: "N", 0: "E"}

    def __init__(self):
        pass

    @classmethod
    def encode_sequence(cls, sequence):
        val = 0
        for i, base in enumerate(sequence):
            val += cls.__mapping[base] * (16**i)
        return val

    @classmethod
    def encode_cluster(cls, cluster):
        ec = []
        N = len(cluster)
        indices = []
        lengths = []

        byte_array = N.to_bytes(4, "little")

        for sequence, index in cluster:
            length = (len(sequence) + 1) // 2

            es = cls.encode_sequence(sequence).to_bytes(length, "little")
            ec.append(es)
            indices.append(index)
            lengths.append(length)

        for index in indices:
            byte_array += index.to_bytes(4, "little")

        for length in lengths:
            byte_array += length.to_bytes(1, "little")

        for es in ec:
            byte_array += es

        return byte_array

    @classmethod
    def encode_clusters(cls, clusters):
        out = []
        for cluster in clusters:
            out.append(cls.encode_cluster(cluster))
        return out

    @classmethod
    def decode_sequence(cls, val):
        import math
        sequence = ""
        n = math.floor(math.log(val) / math.log(16))

        while val > 0:
            next_layer = val % 16**n
            sequence = str(cls.__rmapping[int((val - next_layer) / 16**n)]) + sequence
            n -= 1
            val = next_layer

        return sequence

    @classmethod
    def decode_cluster(cls, byte_array):
        N = int.from_bytes(byte_array[:4], "little", signed=False)
        indices = []
        lengths = []
        cluster = []

        for i in range(1, N + 1):
            indices.append(
                int.from_bytes(byte_array[4 * i : 4 * i + 4], "little", signed=False)
            )

        for i in range(4 * N + 4, 5 * N + 4):
            lengths.append(
                int.from_bytes(byte_array[i : i + 1], "little", signed=False)
            )

        indent = 5 * N + 4
        for i in range(N):
            esb = byte_array[indent : indent + lengths[i]]
            es = int.from_bytes(esb, "little", signed=False)
            cluster.append([cls.decode_sequence(es), indices[i]])
            indent += lengths[i]

        return cluster

    @classmethod
    def decode_clusters(cls, clusters):
        out = []
        for cluster in clusters:
            out.append(cls.decode_cluster(cluster))
        return out


if __name__ == "__main__":
    cluster = [
        [["ACGT", 1], ["TCGA", 1], ["ACGTGTCGAGTGT", 2]],
        [["ACGATGCGCGCTAGGT", 1], ["ACGTGTCGCGCAATCGCTAGAC", 2]],
    ]
    encoded = Encoder.encode_clusters(cluster)
    print(encoded)
    decoded = Encoder.decode_clusters(encoded)
    print(decoded)
