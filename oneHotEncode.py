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
    # print(Mapper.encode('GT'))
    R = reader.Reader()
    mp, count, total_len = R.read_fasta("sample/sandmouse_copy.fa")
    sequences_dict_items = mp.values()
    sequences = list(sequences_dict_items)

    # with is like your try .. finally block in this case
    with open("sample/sandmouse_copy.fa", "r") as file:
        # read a list of lines into data
        data = file.readlines()

    even_lines, odd_lines = data[::2], data[1::2]

    with open("SandMouseOneHotNewest2.fa", "w") as f1:
        for odd_line in odd_lines:
            # print(type(odd_line))
            tmp = ""
            for i, char in enumerate(odd_line):
                tmp += str(__mapping[char])
            f1.write(tmp)
        # for even_line in even_lines:
        #     f1.write(even_line)


# encoder = FastaOneHotEncoder(
#     nucleotides="acgt", lower=True, sparse=False, handle_unknown="ignore"
# )
# path = "sample/pathogen.fa"
# encoder.transform_to_df(path, verbose=True).to_csv("oneHot.csv")
