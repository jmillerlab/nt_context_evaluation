# Evaluating how context impacts nucleotide-resolution predictions from the Nucleotide Transformer/SegmentNT in the top 5 most cited genes 

**Step 0: Set Up Environment to Run the Nucleotide Transformer and SegmentNT**

More information about the Nucleotide Transformer and SegmentNT can be found on their GitHub repository: https://github.com/instadeepai/nucleotide-transformer 

**Step 1: Download the Reference Genome**

The same reference genome used for training the Nucleotide Transformer was used for testing the effects of context on the predictions: GRCh38 (GCF_000001405.26).

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz
```

_md5sum_
```
3d03cc56a6a65413b8abc5a3239d78cd  GCF_000001405.26_GRCh38_genomic.fna.gz
17a5bba01d6af11a17822dd661ce4d7f  GCF_000001405.26_GRCh38_genomic.gff.gz
```

**Step 1a: bgzip and index the fna file so samtools can query it**
```
zcat GCF_000001405.26_GRCh38_genomic.fna.gz | bgzip -c >GCF_000001405.26_GRCh38_genomic.fna.bgz
samtools faidx GCF_000001405.26_GRCh38_genomic.fna.bgz
```

_md5sum_
```
73e1c720bf86dec9e9c5a2f5dc6677a6  GCF_000001405.26_GRCh38_genomic.fna.bgz
6a0bdc7a62ed22744b5da347f8835c68  GCF_000001405.26_GRCh38_genomic.fna.bgz.fai
c7ef1d39c5d8402afac445e2e62f42f7  GCF_000001405.26_GRCh38_genomic.fna.bgz.gzi
```

**Step 2: Extract gene annotations from the gff file**

As an example, here is how the APOE gene was extracted:
```
zcat GCF_000001405.26_GRCh38_genomic.gff.gz |grep "gene=APOE;" >apoe_grch38.txt
```

<img width="1039" alt="Screenshot 2025-03-07 at 9 04 49 AM" src="https://github.com/user-attachments/assets/c674ad62-6172-4cc7-ba13-3707e72d17e9" />

**Step 2a: Select the interval that you want to query from apoe_grch38.txt**

In this case, we wanted to look at the entire APOE gene sequence, which apoe_grch38.txt shows is located at NC_000019.10:44905782-44909393. However, since we want to evaluate how context impacts predictions, we padded that sequence with 50k bp upstream and downstream of the APOE gene. This padding allowed us to potentially test contexts up to 50k bp (e.g., when the first nucleotide of the APOE gene is the last nucleotide of a 50k input sequence). We decided to extract the following interval, which includes the 50k padding: NC_000019.10:44855782-44959393.

**Step 3: Create a smaller fasta file with just the APOE sequence**
```
samtools faidx GCF_000001405.26_GRCh38_genomic.fna.bgz "NC_000019.10:44855782-44959393" >apoe.fasta
```
<img width="833" alt="Screenshot 2025-03-07 at 9 06 02 AM" src="https://github.com/user-attachments/assets/aad033f8-3941-4036-9639-25ffe00c45e9" />

Note: Genomic coordinates with +/- 50k bases for the other genes are listed as follows:

EGFR   NC_000007.14:54969032-55257338

TNF    NC_000006.12:31525567-31628336

TP53   NC_000017.11:7618402-7737550

VEGFA  NC_000006.12:43720209-43836487

Negative control (NC_000003.12:153050001-153053598)   NC_000003.12:153000001-153103598

**Step 4: Run run_nt.py**

```
usage: run_nt.py [-h] -i INPUT -o OUTPUT [-n NUMBER_OF_TOKENS_PER_SEQ] [-t TOKEN_SIZE] [-m MAX_SEQS_PER_RUN] [-p PADDING]

Run SegmentNT multiple times using a sliding window of 1 nt. Specify a gene sequence and context size to iterate over all possible contexts.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        (input) Required input fasta file with the gene sequence. Only 1 header and sequence are allowed, although the sequence can be split on multiple lines.
  -o OUTPUT, --output OUTPUT
                        (output) Required output text file to write the predictions.
  -n NUMBER_OF_TOKENS_PER_SEQ, --number_of_tokens_per_seq NUMBER_OF_TOKENS_PER_SEQ
                        Number of tokens per sequence
  -t TOKEN_SIZE, --token_size TOKEN_SIZE
                        Token size. Most models currently use a token size of 6. Do not change unless a model uses a different token length
  -m MAX_SEQS_PER_RUN, --max_seqs_per_run MAX_SEQS_PER_RUN
                        The number of sequences to supply to SegmentNT for each run. Dependant on GPU capabilities
  -p PADDING, --padding PADDING
                        The number of nucleotides before and after the gene sequence that are included in the fasta file
```

_Example Command using 4096 tokens_

```
python run_nt.py -i apoe.fasta -n 4096 -o apoe_nt_output_4096.txt
```

_Example of the first few lines in an input fasta file. The fasta file for a single gene is usually pretty small <0.5MB_
<img width="433" alt="Screenshot 2025-03-07 at 8 04 15 AM" src="https://github.com/user-attachments/assets/562062dd-3191-49a9-ad35-82304670123f" />

_Example of first few lines in output text file. These output files range from several MB to >100G depending on the length of the gene and the context size_
<img width="1482" alt="Screenshot 2025-03-07 at 8 07 40 AM" src="https://github.com/user-attachments/assets/30fbeeab-6593-49b1-9a07-e5376076e50a" />


**Step 5: Run parseWindow.py**

```
usage: parseWindow.py [-h] -i INPUT -o OUTPUT

Parse output from run_nt.py to make it easier to find predictions at any position in the input sequence for exons and introns

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        (input) Required input file from run_nt.py output.
  -o OUTPUT, --output OUTPUT
                        (output) Required output text file
```

_Example Command_

```
python parseWindow.py -i apoe_nt_output_4096.txt -o apoe_ordered_predictions_4096.txt
```

_Screenshot of the first few characters in each line of the output file. These output files are about the same size as the input file_ 
<img width="808" alt="Screenshot 2025-03-07 at 8 11 23 AM" src="https://github.com/user-attachments/assets/b68ddcff-68bf-4070-a076-4fc5a11f2b82" />

**Step 6: Run getPositions.py**

```
usage: getPositions.py [-h] -i INPUT -o OUTPUT [-n NUMBER_OF_TOKENS_PER_SEQ] [-t TOKEN_SIZE]

Get the prediction values across the entire gene sequence for the first, middle, and last nucleotide. Requires output from parseWindow.py

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        (input) Required input file from parseWindow.py output.
  -o OUTPUT, --output OUTPUT
                        (output) Required output text file
  -n NUMBER_OF_TOKENS_PER_SEQ, --number_of_tokens_per_seq NUMBER_OF_TOKENS_PER_SEQ
                        Number of tokens per sequence
  -t TOKEN_SIZE, --token_size TOKEN_SIZE
                        Token size. Most models currently use a token size of 6. Do not change unless a model uses a different token length
 ```
                       
_Example Command_

```
python getPositions.py -i apoe_ordered_predictions_4096.txt -n 4096 -o apoe_first_middle_last_predictions_4096.txt 
```

_Screenshot of first few characters of apoe_first_middle_last_predictions_4096.txt. Expected output file size is <15MB_
<img width="989" alt="Screenshot 2025-03-07 at 8 15 08 AM" src="https://github.com/user-attachments/assets/63c8e558-fc8e-4bbd-a113-842609ade2a7" />

**Step 7: Run getOnePosition.py**

```
usage: getOnePosition.py [-h] -i INPUT -o OUTPUT [-n NUMBER_OF_TOKENS_PER_SEQ] [-t TOKEN_SIZE] [-p POSITION]

Get the prediction values at a single position within the sequence. Output values represent the prediction within the context window (i.e., is the nucleotide at the beginning of the prediction sequence or the
end?). Requires output from parseWindow.py

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        (input) Required input file from parseWindow.py output.
  -o OUTPUT, --output OUTPUT
                        (output) Required output text file
  -n NUMBER_OF_TOKENS_PER_SEQ, --number_of_tokens_per_seq NUMBER_OF_TOKENS_PER_SEQ
                        Number of tokens per sequence
  -t TOKEN_SIZE, --token_size TOKEN_SIZE
                        Token size. Most models currently use a token size of 6. Do not change unless a model uses a different token length
  -p POSITION, --position POSITION
                        Position within the sequence to query
```

_Example Command_

```
python getOnePosition.py -i apoe_ordered_predictions_4096.txt -n 4096 -p 850 -o apoe_at_pos_850_4096.txt
```

Positions used for different genes in the manuscript:

APOE: 850 (exon 2)

VEGFA: 1500 (exon 9; this is the last exon)

TNF: 180 (exon 1)

EGFR: 183535 (exon 29)

TP53: 6800 (exon 4)

Control: 850

_The output file is expected to look similar to the screenshot below. Expected file size <1MB_
<img width="852" alt="Screenshot 2025-03-07 at 8 17 23 AM" src="https://github.com/user-attachments/assets/7a2f09bd-9bf0-46d4-bf63-1de90fda840a" />

**Step 8: Graph the Results**

The scripts to graph results are located in ```make_graphs```. Scripts to calculate the accuracy of the predictions are located in ```calculate_stats```. Please see the accompanying ```README.md``` in each directory for additional information on how to run each of those scripts.


