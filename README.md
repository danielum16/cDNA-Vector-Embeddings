# cDNA-Vector-Embeddings
Repository for the paper "Vector Embeddings by Sequence Similarity and Context for Improved Compression, Similarity Search, Clustering, Organization, and Manipulation of cDNA Libraries."

## Repository Structure
cDNA-Vector-Embeddings/
├── data/
│ ├── alternative_splicing_human_10541.fasta
│ ├── 3D_Structure_1_200_seq_length _17552.fasta
│ ├── binding_site_28052.fasta
├── Context_Embeddings_3D_Plot_Generation_Script.ipynb
├── README.md
├── requirements.txt
└── .gitignore


## Contents

- `data/`: Directory containing the FASTA data files.
- `notebooks/`: Directory containing the Jupyter notebook for analysis.
- `src/`: Directory containing additional scripts or modules used in the analysis.
- `README.md`: This file, providing an overview of the repository.
- `LICENSE`: License for the code in this repository.
- `requirements.txt`: List of dependencies required to run the notebook.
- `.gitignore`: Specifies files and directories to be ignored by Git.

## Getting Started

### Prerequisites

To run the code in this repository, you need to have Python installed. The recommended way to manage dependencies is via a virtual environment. You can use `virtualenv` or `conda` for this purpose.

### Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/danielum16/cDNA-Vector-Embeddings.git
    cd cDNA-Vector-Embeddings
    ```

2. Create and activate a virtual environment:
    ```sh
    python -m venv venv
    source venv/bin/activate   # On Windows use `venv\Scripts\activate`
    ```

3. Install the required dependencies:
    ```sh
    pip install -r requirements.txt
    ```

### Running the Notebook

1. Open Notebook

2. Start Jupyter Notebook:
    ```sh
    jupyter notebook
    ```

3. Open and run the `analysis_notebook.ipynb` notebook to see the analysis and results.

## Description of Code

### FASTA to DataFrame Conversion

The `fasta2df` function in the notebook converts a FASTA file to a pandas DataFrame. Each sequence in the FASTA file is converted to uppercase and stored along with its description.

### Encoding and Decoding Sequences

Functions are provided to encode DNA sequences into integer values and decode them back. These are used to transform the sequence data for analysis.

### Data Preprocessing

The `remove_null_AA` function cleans the sequences by removing null amino acids. The `amin` function extracts unique amino acids from the sequences.

### One-Hot Encoding and Model Training

One-hot encoding is used to transform amino acid data into a suitable format for machine learning models. TensorFlow is used to define and train a simple neural network model for embedding the amino acids into a vector space.

### Visualization

The embeddings are visualized using a 3D scatter plot, where each point represents an amino acid and is color-coded based on its properties.

## Contributing

If you wish to contribute to this project, please fork the repository and use a feature branch. Pull requests are warmly welcome.

## Acknowledgements



