import json
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy import stats

def parseArgs():
    parser = argparse.ArgumentParser(description='Graph violin plots comparing predictions made when the prediction is for the first, middle, or last nucleotide in the context')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input json file containing a dictionary with key:value pairs of probabilities.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output file path. If you want the image saved as an svg, include the extension .svg in the file name')
    parser.add_argument('-g', '--gene', type=str.lower, required=True,choices=['apoe','egfr','tnf','tp53','vegfa','control'],help='gene used to draw exons on plot.')
    parser.add_argument('-f', '--feature', type=str.lower, required=True,choices=['exon','intron'],help='Graph "exon" or "intron" probabilities')
    parser.add_argument('-r', '--reverse', action="store_true", help='Get probabilities of Introns in Exons or Exons in Introns')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print the split lists for each interval')
    parser.add_argument('-z', '--zscore', action="store_true", help='Optional normalization with z score')
    parser.add_argument('-t', '--graph_title', default='How Context Impacts Exon Predictions in ',type=str,help='Change Graph Title')
    args = parser.parse_args()
    args.graph_title= args.graph_title +args.gene.upper()
    return args

def apoe(pred_list):
    exons='''NC_000019.10    BestRefSeq      exon    44905782        44905841
    NC_000019.10    BestRefSeq      exon    44906602        44906667
    NC_000019.10    BestRefSeq      exon    44907760        44907952
    NC_000019.10    BestRefSeq      exon    44908533        44909393'''
    gene_start=44905782
    gene_end=44909393
    strand=1 #positive
    new_pred_list=[]
    exons=exons.split("\n")
    for exon in exons:
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start) - gene_start
        end=int(end) - gene_start +1
        new_pred_list += pred_list[start:end]
    return new_pred_list


def vegfa(pred_list):
    exons='''NC_000006.12    BestRefSeq      exon    43770209        43771312
    NC_000006.12    BestRefSeq      exon    43771985        43772101
    NC_000006.12    BestRefSeq      exon    43774341        43774392
    NC_000006.12    BestRefSeq      exon    43777469        43777665
    NC_000006.12    BestRefSeq      exon    43778460        43778536
    NC_000006.12    BestRefSeq      exon    43778889        43778918
    NC_000006.12    BestRefSeq      exon    43780732        43780854
    NC_000006.12    BestRefSeq      exon    43781956        43782087
    NC_000006.12    BestRefSeq      exon    43784541        43786487'''
    gene_start=43770209
    gene_end=43786487
    strand=1 #positive
    exons=exons.split("\n")
    new_pred_list=[]
    for exon in exons:
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start) - gene_start
        end=int(end) - gene_start +1
        new_pred_list += pred_list[start:end]
    return new_pred_list

def tnf(pred_list):
    exons='''NC_000006.12    BestRefSeq      exon    31575567        31575927
    NC_000006.12    BestRefSeq      exon    31576534        31576579
    NC_000006.12    BestRefSeq      exon    31576767        31576814
    NC_000006.12    BestRefSeq      exon    31577116        31578336'''
    gene_start=31575567
    gene_end=31578336
    strand=1 #positive
    exons=exons.split("\n")
    new_pred_list=[]
    for exon in exons:
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start) - gene_start
        end=int(end) - gene_start +1
        new_pred_list += pred_list[start:end]
    return new_pred_list


def egfr(pred_list):
    exons='''NC_000007.14    BestRefSeq      exon    55019032        55019365
    NC_000007.14    BestRefSeq      exon    55142286        55142437
    NC_000007.14    BestRefSeq      exon    55143305        55143488
    NC_000007.14    BestRefSeq      exon    55146606        55146740
    NC_000007.14    BestRefSeq      exon    55151294        55151362
    NC_000007.14    BestRefSeq      exon    55152546        55152664
    NC_000007.14    BestRefSeq      exon    55154011        55154152
    NC_000007.14    BestRefSeq      exon    55155830        55155946
    NC_000007.14    BestRefSeq      exon    55156533        55156659
    NC_000007.14    BestRefSeq      exon    55156759        55156951
    NC_000007.14    BestRefSeq      exon    55157663        55157753
    NC_000007.14    BestRefSeq      exon    55160139        55160338
    NC_000007.14    BestRefSeq      exon    55161499        55161631
    NC_000007.14    BestRefSeq      exon    55163733        55163823
    NC_000007.14    BestRefSeq      exon    55165280        55165437
    NC_000007.14    BestRefSeq      exon    55168523        55168635
    NC_000007.14    BestRefSeq      exon    55170307        55171045
    NC_000007.14    BestRefSeq      exon    55171175        55171213
    NC_000007.14    BestRefSeq      exon    55172983        55173124
    NC_000007.14    BestRefSeq      exon    55173921        55174043
    NC_000007.14    BestRefSeq      exon    55174722        55174820
    NC_000007.14    BestRefSeq      exon    55181293        55181478
    NC_000007.14    BestRefSeq      exon    55191719        55191874
    NC_000007.14    BestRefSeq      exon    55192766        55192841
    NC_000007.14    BestRefSeq      exon    55198717        55198863
    NC_000007.14    BestRefSeq      exon    55200316        55200413
    NC_000007.14    BestRefSeq      exon    55201188        55201355
    NC_000007.14    BestRefSeq      exon    55201735        55201782
    NC_000007.14    BestRefSeq      exon    55202517        55202625
    NC_000007.14    BestRefSeq      exon    55205256        55207338'''
    gene_start=55019032
    gene_end=55207338
    strand=1 #positive
    exons=exons.split("\n")
    new_pred_list=[]
    for exon in exons:
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start) - gene_start
        end=int(end) - gene_start +1
        new_pred_list += pred_list[start:end]
    return new_pred_list


def tp53(pred_list):
    '''
    Some overlapping exons were removed. Longest exon in overlapping region was chosen.
    '''
    exons='''NC_000017.11    BestRefSeq      exon    7668402 7669690
    NC_000017.11    BestRefSeq      exon    7670609 7670715
    NC_000017.11    BestRefSeq      exon    7673207 7673339
    NC_000017.11    BestRefSeq      exon    7673535 7673608
    NC_000017.11    BestRefSeq      exon    7673701 7673837
    NC_000017.11    BestRefSeq      exon    7674181 7674290
    NC_000017.11    BestRefSeq      exon    7674859 7674971
    NC_000017.11    BestRefSeq      exon    7675053 7675493
    NC_000017.11    BestRefSeq      exon    7675994 7676272
    NC_000017.11    BestRefSeq      exon    7676382 7676622
    NC_000017.11    BestRefSeq      exon    7687377 7687550'''
    gene_start=7668402
    gene_end=7687550
    strand=-1 #on the negative strand, so exon order is reversed
    exons=exons.split("\n")
    new_pred_list=[]
    for exon in exons:
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start) - gene_start
        end=int(end) - gene_start +1
        new_pred_list += pred_list[start:end]
    return new_pred_list



def control(pred_list):
    #No exons need to be drawn
    return pred_list

def makeViolinPlots(args):
    
    geneFunc={'apoe':apoe,'egfr':egfr,'tnf':tnf,'vegfa':vegfa,'control':control,'tp53':tp53}
    with open(args.input) as inputF:
        allValues = json.loads(inputF.read())
    if args.feature == 'exon':
        fml_dict = allValues['Exon']
        if args.reverse:
            args.graph_title = args.graph_title + " Introns"
    elif args.feature == 'intron':
        fml_dict = allValues['Intron']
        args.graph_title = args.graph_title.replace("Exon","Intron")
        if args.reverse:
            args.graph_title = args.graph_title + " Exons"
    else:
        import sys
        print("Invalid args.feature option:", args.feature)
        sys.exit()
        
    assert 'first' in fml_dict, 'input dictionary requires key word "first"'
    assert 'middle' in fml_dict, 'input dictionary requires key word "middle"'
    assert 'last' in fml_dict, 'input dictionary requires key word "last"'
    assert (len(fml_dict['first']) == len(fml_dict['middle'])) and (len(fml_dict['middle'])== len(fml_dict['last'])),"All lists must be the same length"
    
    if args.zscore:
        for key,value in fml_dict.items():
            fml_dict[key] = list(stats.zscore(value))
    for position in ['first','middle','last']:
        exon_list = geneFunc[args.gene] (fml_dict[position]) #gets just the exonic regions
        if (args.feature =="intron" and not args.reverse) or (args.feature == 'exon' and args.reverse):
            for x in exon_list:
                fml_dict[position].remove(x)
        else:
            fml_dict[position] = exon_list
    df= pd.DataFrame(fml_dict)
    
    # Perform statistical tests (e.g., t-tests)
    p_values = {}
    for col1 in df.columns:
        for col2 in df.columns:
            if col1 != col2 and (col2, col1) not in p_values:
                t_stat, p_val = stats.ttest_ind(df[col1], df[col2])
                p_values[(col1, col2)] = p_val
    
    # Create the violin plot
    plt.figure(figsize=(8, 6))
    colors={"last":'forestgreen',"middle":'Red',"first":'Blue'}
    sns.violinplot(data=df,palette=colors)
    
    # Add significance markers
    max_val = df.max().max()
    y_pos = max_val * 1.2 
    for (col1, col2), p_val in p_values.items():
        x1, x2 = df.columns.get_loc(col1), df.columns.get_loc(col2)
        plt.plot([x1, x2], [y_pos, y_pos], color='black', linewidth=1.5)
    
        if p_val == 0:
            symbol="P < 2.2251e-308"

        elif p_val<0.001:
            symbol=f"P = {p_val:.4e}"
        else:
            symbol=f"P = {p_val:.4g}"
        
        plt.text((x1 + x2) / 2, y_pos + max_val*0.02, symbol, ha='center', va='bottom', color='black')
        y_pos += max_val*0.2 
    # Customize the plot
    plt.title(args.graph_title)
    plt.xlabel('Context of Nucleotide')
    if args.feature == 'exon':
        plt.ylabel('Predicted Probability of Exon')
    elif args.feature == 'intron':
        plt.ylabel('Predicted Probability of Intron')
    plt.tight_layout()
    plt.savefig(args.output)
    plt.show()
    

def main():
    args=parseArgs()
    makeViolinPlots(args)

if __name__ == '__main__':
    main()

