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

_The expected output of example.html is shown below_
<img width="1155" alt="Screenshot 2025-03-07 at 9 10 31 AM" src="https://github.com/user-attachments/assets/5e1d83d8-df06-45a4-8365-fe476a6a0a5a" />


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

_**Each plot is shown below**_

_example_subtraction_Last.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 14 48 AM" src="https://github.com/user-attachments/assets/2ef2d13f-c115-405b-928c-ac88450df2bc" />


_example_subtraction_Middle.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 15 06 AM" src="https://github.com/user-attachments/assets/4741dcf3-1764-466d-8206-3fae3d11f1ea" />


_example_subtraction_First.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 14 15 AM" src="https://github.com/user-attachments/assets/c287dd5a-f667-41d0-9bb3-692652f91eac" />


_example_zscore_Last.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 15 36 AM" src="https://github.com/user-attachments/assets/5f316c13-a9f4-4911-ac72-faef0d1501a0" />


_example_zscore_Middle.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 15 50 AM" src="https://github.com/user-attachments/assets/caa9b184-0ad7-4a0f-9a99-9fb2b0b0eb66" />


_example_zscore_First.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 15 20 AM" src="https://github.com/user-attachments/assets/946e5f68-a974-49b4-8e6b-424fe8e733ee" />


_example_raw_Last.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 13 28 AM" src="https://github.com/user-attachments/assets/36609159-08a0-46e5-b419-da4d617fec1f" />


_example_raw_Middle.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 13 58 AM" src="https://github.com/user-attachments/assets/eef65a2f-ce9b-4462-bc14-dbab211c5701" />


_example_raw_First.html_
<img width="1470" alt="Screenshot 2025-03-07 at 9 13 42 AM" src="https://github.com/user-attachments/assets/7fd9776a-ccb3-490b-933b-e1c6f5815af2" />


