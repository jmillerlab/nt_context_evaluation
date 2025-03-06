import plotly.express as px
import pandas as pd
import argparse
import json
import plotly.offline as pio
from statistics import mean,stdev,median
import numpy as np
import sys


def parseArgs():
    parser = argparse.ArgumentParser(description='Graph the difference between exon/intron probabilities across a gene when predicting the first, middle, or last nucleotide')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input json file containing a dictionary with key:value pairs of probabilities.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output file path. Example: APOE_diff_exon_intron.html.')
    parser.add_argument('-g', '--gene', type=str.lower, required=True,choices=['apoe','egfr','tnf','tp53','vegfa','control'],help='gene used to draw exons on plot.')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print the split lists for each interval')
    parser.add_argument('-z', '--zscore', action="store_true", help='Optional normalization with z score')
    parser.add_argument('-n', '--subtract_intron_from_exon', action="store_true", help='Optional normalize by subtracting intron prediction from exon prediction')
    parser.add_argument('-s', '--svg', action="store_true", help='Optional Save SVG file with same name as html. Will save to Downloads directory')
    parser.add_argument('-t', '--graph_title', default=' Difference Between Exon and Intron Predictions Using Difference Contexts',type=str,help='Change Graph Title')
    parser.add_argument('-m', '--markers', action="store_true", help='Draw the plot with the point markers visible')
    args = parser.parse_args()
    args.graph_title= args.gene.upper() + args.graph_title
    return args

def addExon(fig,num,x_start,x_end):
        # Add Exon 1
        fig.add_shape(type='line',
                      x0=x_start,
                      y0=1.0,
                      x1=x_end,
                      y1=1.0,
                      line=dict(color='Black',),
                      name=str(x_start) +"-" +str(x_end),
                      xref='x',
                      yref='y')
       
        middle_annotate=x_start + int((x_end-x_start)/2)-3
        fig.add_annotation(
            x=middle_annotate, y=1.02,
            text="Exon " +str(num),
            showarrow=False,
            font=dict(size=14)
        )
        return fig

def apoe(fig):
    '''
    NC_000019.10    Gnomon  exon    44905807        44905923
    NC_000019.10    Gnomon  exon    44906602        44906667
    NC_000019.10    Gnomon  exon    44907760        44907952
    NC_000019.10    Gnomon  exon    44908533        44909393
    '''
    gene_start=44905782
    gene_end=44909393
    strand=1 #positive

    fig=addExon(fig,1,1,60)
    fig=addExon(fig,2,821,886)
    fig=addExon(fig,3,1979,2171)
    fig=addExon(fig,4,2752,3612)
    fig.update_layout(xaxis=dict(range=[0,gene_end-gene_start+1]))
    return fig

def vegfa(fig):
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
    num=0
    for exon in exons:
        num+=1
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start)
        end=int(end)
        fig=addExon(fig,num,start-gene_start,end-gene_start)
    fig.update_layout(xaxis=dict(range=[0,gene_end-gene_start+1]))
    return fig

def tnf(fig):
    exons='''NC_000006.12    BestRefSeq      exon    31575567        31575927
    NC_000006.12    BestRefSeq      exon    31576534        31576579
    NC_000006.12    BestRefSeq      exon    31576767        31576814
    NC_000006.12    BestRefSeq      exon    31577116        31578336'''
    gene_start=31575567
    gene_end=31578336
    strand=1 #positive
    exons=exons.split("\n")
    num=0
    for exon in exons:
        num+=1
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start)
        end=int(end)
        fig=addExon(fig,num,start-gene_start,end-gene_start)
    fig.update_layout(xaxis=dict(range=[0,gene_end-gene_start+1]))
    return fig


def egfr(fig):
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
    num=0
    for exon in exons:
        num+=1
        exon=exon.strip().split()
        start,end = exon[3:]
        start =int(start)
        end=int(end)
        fig=addExon(fig,num,start-gene_start,end-gene_start)
    fig.update_layout(xaxis=dict(range=[0,gene_end-gene_start+1]))
    return fig


def tp53(fig):
    '''
    Some overlapping exons were removed. Longest exon in overlapping region was chosen.
    NC_000017.11    BestRefSeq      exon    7668402 7669690
    NC_000017.11    BestRefSeq      exon    7670609 7670715
    NC_000017.11    BestRefSeq      exon    7673207 7673339
    NC_000017.11    BestRefSeq      exon    7673535 7673608
    NC_000017.11    BestRefSeq      exon    7673701 7673837
    NC_000017.11    BestRefSeq      exon    7674181 7674290
    NC_000017.11    BestRefSeq      exon    7674859 7674971
    NC_000017.11    BestRefSeq      exon    7675053 7675493
    NC_000017.11    BestRefSeq      exon    7675994 7676272
    NC_000017.11    BestRefSeq      exon    7676382 7676622
    NC_000017.11    BestRefSeq      exon    7687377 7687550
    '''
    gene_start=7668402
    gene_end=7687550
    strand=-1 #on the negative strand, so exon order is reversed
    fig=addExon(fig,1,7687377-gene_start,7687550-gene_start)
    fig=addExon(fig,2,7676382-gene_start,7676622-gene_start)
    fig=addExon(fig,3,7675994-gene_start,7676272-gene_start)
    fig=addExon(fig,4,7675053-gene_start,7675493-gene_start)
    fig=addExon(fig,5,7674859-gene_start,7674971-gene_start)
    fig=addExon(fig,6,7674181-gene_start,7674290-gene_start)
    fig=addExon(fig,7,7673701-gene_start,7673837-gene_start)
    fig=addExon(fig,8,7673535-gene_start,7673608-gene_start)
    fig=addExon(fig,9,7673207-gene_start,7673339-gene_start)
    fig=addExon(fig,10,7670609-gene_start,7670715-gene_start)
    fig=addExon(fig,11,7668402-gene_start,7669690-gene_start)
    fig.update_layout(xaxis=dict(range=[0,gene_end-gene_start+1]))
    return fig



def control(fig):
    #No exons need to be drawn
    fig.update_layout(xaxis=dict(range=[0,3613]))
    return fig

def formatFig(fig,graph_title,x_axis_title):
    '''
    Change fig layout
    '''
    # Update layout to make it look prettier
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',  # Set the paper background color to transparent
        plot_bgcolor='rgba(0,0,0,0)',
        title=graph_title,
        title_font_color='black',
        title_x=0.5,
        xaxis_title=x_axis_title,
        yaxis_title=f'Difference in Exon and Intron Predictions',
        legend_title='Position of Predicted Nucleotide',
        xaxis=dict(
            tickformat=",",  # Use comma as thousands separator
            title_font_color="black",
            tickfont=dict(color="black"),
            color='black',  # Set x-axis color to black
            visible=True,  # Ensure the x-axis is visible
            ticks="outside",  # Add tick marks
            showline=True,
            linecolor='black',
            linewidth=2
        ),
        yaxis=dict(
            title="Probability",
            title_font_color="black",
            tickfont=dict(color="black"),
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks="outside"  # Add tick marks
        ),
        legend=dict(font=dict(color='black'))
    )
    fig.add_hline(y=0, line_dash="dash", line_color="black")
    #fig.update_layout(yaxis=dict(range=[0, 1.025])) #Makes y-axis start at 0. Reduces white space at bottom of graph
    ###fig.update_traces(line=dict(width=0.5)) #Makes line width thinner.
    return fig


def drawAUC(predictions,predictions2,predictions3,labels,gene):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score
    
    # Sample true labels (0 or 1)
    y_true = np.array(labels)
    
    # Corresponding predicted probabilities for the positive class
    y_prob = np.array(predictions)
    y_prob2 = np.array(predictions2)
    y_prob3 = np.array(predictions3)
    
    # Calculate the ROC curve
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    fpr2, tpr2, thresholds2 = roc_curve(y_true, y_prob2)
    fpr3, tpr3, thresholds3 = roc_curve(y_true, y_prob3)
    
    # Calculate the AUC score
    roc_auc = roc_auc_score(y_true, y_prob)
    roc_auc2 = roc_auc_score(y_true, y_prob2)
    roc_auc3 = roc_auc_score(y_true, y_prob3)
    '''
    Optionally print all of the AUC values
    print(f"AUC f/m/l {gene} {roc_auc:.4f} {roc_auc2:.4f} {roc_auc3:.4f}")
    '''
    
    # Plot the ROC curve
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'First ROC curve (AUC = {roc_auc:.4f})')
    plt.plot(fpr2, tpr2, color='red', lw=2, label=f'Middle ROC curve (AUC = {roc_auc2:.4f})')
    plt.plot(fpr3, tpr3, color='green', lw=2, label=f'Last ROC curve (AUC = {roc_auc3:.4f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'{gene} Receiver Operating Characteristic (ROC)')
    plt.legend(loc="lower right")
    plt.savefig(gene + "_auc_plot.svg",transparent=True)
    ###plt.show()


def getGeneExons():
    ##Format {gene_name : [exon_locations, first exon start site]}
    gene_exons={'tnf':['''NC_000006.12    BestRefSeq      exon    31575567        31575927
    NC_000006.12    BestRefSeq      exon    31576534        31576579
    NC_000006.12    BestRefSeq      exon    31576767        31576814
    NC_000006.12    BestRefSeq      exon    31577116        31578336''',31575567,31578336],
                'tp53': ['''NC_000017.11    BestRefSeq      exon    7668402 7669690
    NC_000017.11    BestRefSeq      exon    7670609 7670715
    NC_000017.11    BestRefSeq      exon    7673207 7673339
    NC_000017.11    BestRefSeq      exon    7673535 7673608
    NC_000017.11    BestRefSeq      exon    7673701 7673837
    NC_000017.11    BestRefSeq      exon    7674181 7674290
    NC_000017.11    BestRefSeq      exon    7674859 7674971
    NC_000017.11    BestRefSeq      exon    7675053 7675493
    NC_000017.11    BestRefSeq      exon    7675994 7676272
    NC_000017.11    BestRefSeq      exon    7676382 7676622
    NC_000017.11    BestRefSeq      exon    7687377 7687550''', 7668402, 7687550],
                'egfr': ['''NC_000007.14    BestRefSeq      exon    55019032        55019365
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
    NC_000007.14    BestRefSeq      exon    55205256        55207338''',55019032,55207338],
                'vegfa':['''NC_000006.12    BestRefSeq      exon    43770209        43771312
    NC_000006.12    BestRefSeq      exon    43771985        43772101
    NC_000006.12    BestRefSeq      exon    43774341        43774392
    NC_000006.12    BestRefSeq      exon    43777469        43777665
    NC_000006.12    BestRefSeq      exon    43778460        43778536
    NC_000006.12    BestRefSeq      exon    43778889        43778918
    NC_000006.12    BestRefSeq      exon    43780732        43780854
    NC_000006.12    BestRefSeq      exon    43781956        43782087
    NC_000006.12    BestRefSeq      exon    43784541        43786487''',43770209,43786487],
    'apoe': ['''NC_000019.10    BestRefSeq  exon    44905782        44905841
    NC_000019.10    BestRefSeq  exon    44906602        44906667
    NC_000019.10    BestRefSeq  exon    44907760        44907952
    NC_000019.10    BestRefSeq  exon    44908533        44909393''',44905782,44909393]
                }
   
    return gene_exons


def write_stats(gene_name,position,in_exon_exon,in_intron_exon_preds,predictions,threshold):
     gene_name=gene_name.upper()
     output=open(f"stats_{position}.txt",'a')
     stats = calculate_stats(in_exon_exon)
     output.write(f"Exon IQR {position} {gene_name} {stats['iqr']}\n")
     output.write(f"Average exon score in exons for {gene_name}: {stats['mean']}±{stats['sd']} ({stats['median']})\n")
     stats = calculate_stats(in_intron_exon_preds)
     output.write(f"Average exon score in introns for {gene_name}: {stats['mean']}±{stats['sd']} ({stats['median']})\n")
     totalAboveHalfPreds =sum(i >= threshold for i in predictions) 
     tp=round(sum(i >= threshold for i in in_exon_exon),4)
     fp=round(totalAboveHalfPreds-tp,4)
     fn=round(len(in_exon_exon)-tp,4)
     tn=round(len(in_exon_exon)+len(in_intron_exon_preds)-len(in_exon_exon)-fp,4)
     output.write(f"TP, FP,FN, TN, {gene_name}, {tp}, {fp}, {fn}, {tn}\n")
     sensitivity=round((tp / (tp +fn)),4)
     specificity=round((tn/(tn+fp)),4)
     accuracy=round(((tp+tn)/(tp+tn+fp+fn)),4)
     output.write(f"Sensitivity, Specificity, Accuracy, {gene_name}, {sensitivity}, {specificity}, {accuracy}\n")

def calcDiff(gene_name,allVals,threshold,subtract):
    gene_exons=getGeneExons()
    exon_minus_intron=dict()
    middle_exon_preds=[]
    for position in ['first','middle','last']:
        exon_preds= allVals['Exon'][position]
        in_intron_exon_preds= exon_preds[:]
        intron_preds= allVals['Intron'][position]
        in_exon_exon=[]
        exons=gene_exons[gene_name][0].split("\n")
        gene_start=int(gene_exons[gene_name][1])
        num=0
        diffVals=[a - b for a, b in zip(exon_preds, intron_preds)]  #calculates difference between introns and exons for all positions
        if subtract:
            exon_preds=diffVals[:]
            in_intron_exon_preds=diffVals[:]
        for exon in exons:
            num+=1
            exon=exon.strip().split()
            start,end = exon[3:]
            start =int(start)
            end=int(end) +1 # +1 needed to make sure it's inclusive, not exclusive
            in_exon_exon_preds=exon_preds[(start-gene_start) : (end-gene_start)]
            in_exon_exon += in_exon_exon_preds
            for x in in_exon_exon_preds:
                in_intron_exon_preds.remove(x)
        if position=='first':
            exon_minus_intron['first']=diffVals
            predictions= in_exon_exon + in_intron_exon_preds
            labels =[1]*len(in_exon_exon) + [0] *len(in_intron_exon_preds)
            write_stats(gene_name,"first",in_exon_exon,in_intron_exon_preds,predictions,threshold) 
        if position=='middle':
            middle_exon_preds=in_exon_exon[:]
            middle_intron_exon_preds=in_intron_exon_preds[:]
            exon_minus_intron['middle']=diffVals
            predictions2= in_exon_exon+ in_intron_exon_preds
            write_stats(gene_name,"middle",in_exon_exon,in_intron_exon_preds,predictions2,threshold) 

        if position=='last':
            exon_minus_intron['last']=diffVals
            predictions3= in_exon_exon+ in_intron_exon_preds
            write_stats(gene_name,"last",in_exon_exon,in_intron_exon_preds,predictions3, threshold) 
    drawAUC(predictions,predictions2,predictions3,labels,gene_name.upper())
    return exon_minus_intron,middle_exon_preds,middle_intron_exon_preds

def calculate_stats(data):
  """Calculates the mean, standard deviation, and median of a list.

  Args:
    data: A list of numerical values.

  Returns:
    A dictionary containing the mean, standard deviation, and median.
  """
  if not data:
      return {"mean": None, "sd": None, "median": None, "iqr": None}
  q75, q25 = np.percentile(data, [75, 25])
  iqr = q75 - q25

  return {
      "mean":   round(mean(data),4),
      "sd":     round(stdev(data),4),
      "median": round(median(data),4),
      "iqr": round(iqr,4)
  }

def plotVals(args):
    threshold=0.5
    with open(args.input) as inputF:
        allValues = json.loads(inputF.read())
        if args.zscore:
            threshold=0.0
            from scipy.stats import zscore
            newDict=dict()
            for key,values in allValues.items():
                df=pd.DataFrame(values)
                df = df.apply(zscore)
                new_values=df.to_dict('list')
                newDict[key] = new_values
            allValues=newDict


    drawExon={'apoe':apoe,'egfr':egfr,'tnf':tnf,'vegfa':vegfa,'control':control,'tp53':tp53}
    if args.subtract_intron_from_exon:
        threshold=0.0
    difference_preds,middle_exon_preds,middle_intron_exon_preds = calcDiff(args.gene,allValues,threshold,args.subtract_intron_from_exon)
    df= pd.DataFrame(difference_preds)
    #Make the figure
    fig = px.line(df,x=range(1,len(difference_preds['first'])+1), y=['first','middle','last'], log_y=False,markers=args.markers)
    gene=args.gene.upper() 
    graph_title=f'{gene} Exon Minus Intron Predictions from Different Contexts'
    x_axis_title=f'{gene} Gene Sequence'
    
    #Add Exon label to graph
    fig = drawExon[args.gene](fig)
    #Format Fig
    fig =formatFig(fig,graph_title,x_axis_title)
    
    #Write to output files
    fig.write_html(args.output)
    if args.svg:
        #Save svg file too. Will save to Downloads.
        pio.plot(fig,image_filename=args.output,image="svg")
    if args.verbose:
        print(f"Wrote {args.output}") #tells the user where the output files were written

def main():
    args=parseArgs()
    plotVals(args)

if __name__ == '__main__':
    main()

