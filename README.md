# Determination of the composition of a metagenic sample

Program is used for classification of metagenic samples based on k-mer frequencies.
It uses provided dataset given in `fasta` and `fastq` format to classify the samples based on cosine similarity of k-mer frequencies.

## Instructions on setting up the environment
```bash
git clone https://github.com/TXort/bioinformatics-project
cd bioinformatics-project
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
## Running the program
Run `python main.py -h` for instructions on how to run the script.

## Additional information
Scripts `confusion_calc.py` and `csv_analyzer.py` are not required by main script and are used to provide additional information about the results of the main script.

`confusion_calc.py` is used to calculate confusion matrix of files located in `references` directory.

`csv_analyzer.py` is placed in a folder where results are generated and is used to analyze the results of the main script.
