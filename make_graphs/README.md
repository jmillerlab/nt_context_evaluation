# Make Graphs
**graph_pred_at_position.py** 

This script shows how the position of the nucleotide within the input sequence impacts the predicted probability of being an exon or an intron. Similarly, different intervals can be specified, which show how the start site of the input sequence impacts the predictions as well (i.e., moving the input sequence n nucleotides to the left or right impacts predictions for any nucleotide within the sequence)

```
usage: graph_pred_at_position.py [-h] -i INPUT -o OUTPUT [-n INTERVAL] [-p POSITION] [-s START] [-e END] [-v] [-z] [-g {apoe,egfr,tnf,tp53,vegfa,control}] [-t GRAPH_TITLE] [-oe] [-oi] [-m]

Graph the exon/intron probabilities for a single nucleotide using different contexts

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     (input) Required input file with two lines: exon=[list of probabilities] and intron=[list of probabilities]
  -o, --output OUTPUT   (output) Required output file path. Example: APOE_24.html
  -n, --interval INTERVAL
                        Interval size
  -p, --position POSITION
                        Position of nucleotide in gene
  -s, --start START     Start of visualization window
  -e, --end END         End of visualization window
  -v, --verbose         flag to print the split lists for each interval
  -z, --zscore          Optional normalization with z score
  -g, --gene {apoe,egfr,tnf,tp53,vegfa,control}
                        gene used in graph
  -t, --graph_title GRAPH_TITLE
                        Change Graph Title
  -oe, --only_exon      Graph only the exon values. Note: Input file must still include the intron list
  -oi, --only_intron    Graph only the intron values. Note: Input file must still include the exon list
  -m, --markers         Draw the plot with the point markers visible
```

_Example Command_

```
python graph_pred_at_position.py -i apoe_at_pos_850_4096.txt -p 850 -g apoe -o apoe_850.html 
```

_The expected output will be an html file looking like the screeshot below. An svg will also be downloaded of the same image_

![image](https://github.com/user-attachments/assets/cb0a0db4-a8cb-47e0-a9a0-085ec2a59fc5)


**graph_different_contexts.py**

```
usage: graph_different_contexts.py [-h] -i INPUT -o OUTPUT -g {apoe,egfr,tnf,tp53,vegfa,control} [-v] [-z] [-s] [-t GRAPH_TITLE] [-m]

Graph the exon/intron probabilities across APOE when predicting the first, middle, or last nucleotide

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     (input) Required input json file containing a dictionary with key:value pairs of probabilities.
  -o, --output OUTPUT   (output) Required output file path. Example: APOE_first_middle_last_graph.html. This will be appended to the name of the feature. For example, if the first feature is "Exon" the file
                        path to the graph with that first feature will be "Exon_APOE_first_middle_last_graph.html"
  -g, --gene {apoe,egfr,tnf,tp53,vegfa,control}
                        gene used to draw exons on plot.
  -v, --verbose         flag to print the split lists for each interval
  -z, --zscore          Optional normalization with z score
  -s, --svg             Optional Save SVG file with same name as html. Will save to Downloads directory
  -t, --graph_title GRAPH_TITLE
                        Change Graph Title
  -m, --markers         Draw the plot with the point markers visible
```

_Example Command_

```
python graph_different_contexts.py -i apoe_first_middle_last_predictions.txt -o example_first_middle_last -g apoe
```

_Two html files will be created looking like the screeshots below._
_File 1: example_first_middle_last_Exon.html_
<img width="1487" alt="Screenshot 2025-03-07 at 8 42 26 AM" src="https://github.com/user-attachments/assets/7d221c1e-7bb9-41b2-bae4-b2552be3dbe5" />

_File 2: example_first_middle_last_Intron.html_
<img width="1487" alt="Screenshot 2025-03-07 at 8 43 08 AM" src="https://github.com/user-attachments/assets/cd04e89a-e7ac-4f38-9963-6e43f8ad4f6c" />

**make_violin_plots.py**
```
usage: makeViolinPlots.py [-h] -i INPUT -o OUTPUT -g {apoe,egfr,tnf,tp53,vegfa,control} -f {exon,intron} [-r] [-v] [-z] [-t GRAPH_TITLE]

Graph violin plots comparing predictions made when the prediction is for the first, middle, or last nucleotide in the context

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     (input) Required input json file containing a dictionary with key:value pairs of probabilities.
  -o, --output OUTPUT   (output) Required output file path. If you want the image saved as an svg, include the extension .svg in the file name
  -g, --gene {apoe,egfr,tnf,tp53,vegfa,control}
                        gene used to draw exons on plot.
  -f, --feature {exon,intron}
                        Graph "exon" or "intron" probabilities
  -r, --reverse         Get probabilities of Introns in Exons or Exons in Introns
  -v, --verbose         flag to print the split lists for each interval
  -z, --zscore          Optional normalization with z score
  -t, --graph_title GRAPH_TITLE
                        Change Graph Title
```

_Example Command_

```
python makeViolinPlots.py -i apoe_first_middle_last_predictions.txt -o violin_plots.svg -g apoe -f exon
```

_One output file is expected: An svg like the one shown below_
<img width="793" alt="Screenshot 2025-03-07 at 8 55 58 AM" src="https://github.com/user-attachments/assets/d43c790c-8c6a-47db-9a37-393ce397d80c" />
