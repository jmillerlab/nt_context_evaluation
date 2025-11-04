import numpy as np
import json
import argparse
import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt

def parseArgs():
    parser = argparse.ArgumentParser(description='Calculate differences between oscillations using context file that has two lines with Exon=[...] and Intron=[...]')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input json file containing a dictionary with key:value pairs of probabilities.')
    parser.add_argument('-n', '--number_remove', type=int, default=4800, required=False, help='Number of predictions to remove from both ends to account for variability in tails.')
    parser.add_argument('-w', '--window', type=int, default=24,required=False,help='Window size to smooth predictions. Default: 24 nt')
    args = parser.parse_args()
    assert args.number_remove % args.window == 0,'args.number_remove must be divisible by args.window'
    return args

def makeGroups(preds,window):
    # 2. Divide the list into 24 groups
    # Elements at index i, i+24, i+48 etc. go into the same group (group i).
    groups = []
    for i in range(window):
        # Use list slicing with a step of 24: data_list[start:stop:step]
        group_data = preds[i::window]
        groups.append(group_data)
    # 'groups' is now a list of 24 lists, where each inner list contains the data for one group.
    return groups

def anova(groups):
    # You can unpack the groups for f_oneway using the * operator
    f_statistic, p_value = f_oneway(*groups)
    print(f"ANOVA F-statistic: {f_statistic:.3f}")
    print(f"ANOVA p-value: {p_value:.3f}")


def tukey(groups):
    # 3. Perform Tukey's HSD test (post-hoc test)
    # This requires a slightly different data format: a single column of all data and a corresponding column of group labels.
    
    # Create a single array of all data
    all_data = np.concatenate(groups)
    
    # Create a corresponding array of group labels
    group_labels = np.repeat(range(1, 25), len(groups[0])) # Labels 1 to 24
    
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(endog=all_data, groups=group_labels, alpha=0.05)
    
    print("\nTukey HSD Results:")
    print(tukey_result)
    
    # Visualize the results using the plot_simultaneous() method
    #tukey_result.plot_simultaneous()
    #plt.show()

def main():
    args=parseArgs()
    with open(args.input) as inputF:
        line1=inputF.readline().strip()
        assert line1.startswith("exon=["),'First line must start with exon=['
        preds= json.loads(line1[5:])
        exons=makeGroups(preds[args.number_remove:(-1 * args.number_remove)],args.window) #remove variable tails (~20% on either end)
        print("Exon Comparison:")
        anova(exons)
        tukey(exons)
        #outputF.write("exon=" +json.dumps(exons) +"\n")
        line2=inputF.readline().strip()
        assert line2.startswith("intron=["),'First line must start with exon=['
        preds= json.loads(line2[7:])
        introns=makeGroups(preds[args.number_remove:(-1*args.number_remove)],args.window)
        print("\n\n\n\nIntron Comparison:")
        anova(introns)
        tukey(introns)
        #outputF.write("intron=" +json.dumps(introns) +"\n")

        

if __name__ == '__main__':
    main()

