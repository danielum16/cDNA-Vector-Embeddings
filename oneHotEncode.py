from fasta_one_hot_encoder import FastaOneHotEncoder
from library import reader
import tensorflow_io as tfio

__mapping = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "\n": [1, 1, 1, 1],
}

if __name__ == "__main__":
    R = reader.Reader()
    mp, count, total_len = R.read_fasta("sample/sandmouse_copy.fa")
    sequences_dict_items = mp.values()
    sequences = list(sequences_dict_items)

    with open("sample/sandmouse_copy.fa", "r") as file:
        data = file.readlines()

    even_lines, odd_lines = data[::2], data[1::2]

    with open("SandMouseOneHotNewest2.fa", "w") as f1:
        for odd_line in odd_lines:
            tmp = ""
            for i, char in enumerate(odd_line):
                tmp += str(__mapping[char])
            f1.write(tmp)
