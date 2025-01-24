# nt_context_evaluation
Evaluating how context impacts predictions from the Nucleotide Transformer/SegmentNT

**Step 1: Download the Reference Genome**

The same reference genome used for training the Nucleotide Transformer was used for testing the effects of context on the predictions: GRCh38 (GCF_000001405.26).

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz
```
**Step 1a: bgzip and index the fna file so samtools can query it**
```
zcat GCF_000001405.26_GRCh38_genomic.fna.gz | bgzip -c >GCF_000001405.26_GRCh38_genomic.fna.bgz
samtools faidx GCF_000001405.26_GRCh38_genomic.fna.bgz
```
**Step 2: Extract APOE annotations from the gff file**
```
zcat GCF_000001405.26_GRCh38_genomic.gff.gz |grep "APOE" >apoe_grch38.txt
```
**Step 2a: Select the interval that you want to query from apoe_grch38.txt**

In this case, we wanted to look at the entire APOE gene sequence, which apoe_grch38.txt shows is located at NC_000019.10:44905782-44909393. However, since we want to evaluate how context impacts predictions, we padded that sequence with 50k bp upstream and downstream of the APOE gene. This padding allowed us to potentially test contexts up to 50k bp (e.g., when the first nucleotide of the APOE gene is the last nucleotide of a 50k input sequence). We decided to extract the following interval, which includes the 50k padding: NC_000019.10:44855782-44959393.

**Step 3: Create a smaller fasta file with just the APOE sequence**
```
samtools faidx GCF_000001405.26_GRCh38_genomic.fna.bgz "NC_000019.10:44855782-44959393" >apoe.fasta
```
**Step 4: Run run_nt_apoe.py**
```
python run_nt_apoe.py >apoe_sliding_window_output.txt
```
**Step 5: Run parseWindow.py**
```
python parseWindow.py > apoe_ordered_predictions.txt
```
**Step 6: Run getPositions.py**
```
python getPositions.py >apoe_first_middle_last_predictions.txt
```
**Step 7: Run getOnePosition.py**
```
python getOnePosition.py >apoe_850.txt
```
**Step 8: Graph the Results**

**Step 8a: graph_pred_at_position.py** 

This script shows how the position of the nucleotide within the input sequence impacts the predicted probability of being an exon or an intron. Similarly, different intervals can be specified, which show how the start site of the input sequence impacts the predictions as well (i.e., moving the input sequence n nucleotides to the left or right impacts predictions for any nucleotide within the sequence)

_Usage_

usage: graph_pred_at_position.py [-h] -i INPUT -o OUTPUT [-n INTERVAL] [-v] [-z] [-t GRAPH_TITLE] [-oe] [-oi]

Graph the exon/intron probabilities for a single nucleotide using different contexts

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     (input) Required input file with two lines: exon=[list of probabilities] and intron=[list of probabilities]
  -o, --output OUTPUT   (output) Required output file path. Example: APOE_24.html
  -n, --interval INTERVAL
                        Interval size
  -v, --verbose         flag to print the split lists for each interval
  -z, --zscore          Optional normalization with z score
  -t, --graph_title GRAPH_TITLE
                        Change Graph Title
  -oe, --only_exon      Graph only the exon values. Note: Input file must still include the intron list
  -oi, --only_intron    Graph only the intron values. Note: Input file must still include the exon list

_Example Command_
```
python graph_pred_at_position.py -i apoe_at_pos_850.txt -o example_output_apoe_850.html 
```
**Step 8b: graph_different_contexts.py**

This script is intended to show the difference in exon and intron prediction probabilities when predicting the first, middle, or last nucleotide in the input sequence across the entire APOE gene.

usage: graph_different_contexts.py [-h] -i INPUT -o OUTPUT [-v] [-z] [-t GRAPH_TITLE]

Graph the exon/intron probabilities across APOE when predicting the first, middle, or last nucleotide

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     (input) Required input json file containing a dictionary with key:value pairs of probabilities.
  -o, --output OUTPUT   (output) Required output file path. Example: APOE_first_middle_last_graph.html. This will be appended to the name of the feature. For example, if the first feature is "Exon" the file path to the graph with that first feature will be "Exon_APOE_first_middle_last_graph.html"
  -v, --verbose         flag to print the split lists for each interval
  -z, --zscore          Optional normalization with z score
  -t, --graph_title GRAPH_TITLE
                        Change Graph Title

_Example Command_
```
python graph_different_contexts.py -i apoe_first_middle_last_predictions.txt -o example_first_middle_last.html 
```
