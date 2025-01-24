# nt_context_evaluation
Evaluating how context impacts predictions from the Nucleotide Transformer/SegmentNT

**graph_pred_at_position.py** 

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
python graph_pred_at_position.py -i apoe_at_pos_850.txt -o example_output.html 
