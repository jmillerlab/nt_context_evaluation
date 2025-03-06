# Calculate statistics and create graphics for the main manuscript and supplement

**Calculate AUC and Average Exon/Intron Lengths**

Example:
```
python graph_diff_Exon_Intron.py -i example_tnf_sliding_window_4096.txt -o example.html -g tnf
```

That command will create an output file called ```example.html``` that graphs the first, middle, and last exon predictions across the entire TNF gene.

Four other output files will also be created or appended to:
```TNF_auc_plot.svg```: The AUC plot with AUC values calculated for the first, middle, and last positions
```stats_last.txt```: Statistics for the last prediction including average±standard deviation (median) exon predictions in exons and intron, true positive, false positive, false negative, true negative, sensitivity, specificy, and accuracy. 
```stats_middle.txt``` The same as above but for the middle prediction
```stats_first.txt``` The same as above but for the predictions of the first nucleotide in a sequence

**Calculate AUC and Average Exon/Intron Lengths**

Example 1:
```
python graph_change_in_predictions_based_on_length.py -i example_exon_probabilities_at_different_lengths.txt -o example -g exon
```

That command will create ```example.html``` which plots the average exon probabilities for each input sequence length. The mean and standard deviations can be calculated by running ```graph_diff_Exon_Intron.py``` on all different contexts and formatting the table as shown in ```example_exon_probabilities_at_different_lengths.txt```

Example 2:
```
python graph_change_in_predictions_based_on_length.py -i example_auc_at_all_lengths.txt -o example -g auc
```

Similarly, ```example_auc_at_all_lengths.txt``` can be created by running ```graph_diff_Exon_Intron.py``` and putting the AUC values in a table of the same format as ```example_auc_at_all_lengths.txt```

That command will produce 9 output files:
```
example_subtraction_Last.html
example_subtraction_Middle.html
example_subtraction_First.html
example_zscore_Last.html
example_zscore_Middle.html
example_zscore_First.html
example_raw_Last.html
example_raw_Middle.html
example_raw_First.html
```

Each file plots the AUC values from different context lengths given the normalization method (raw, z score, or subtraction) and for the context of the predicted nucleotide (First, Middle, or Last).

