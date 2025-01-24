import plotly.express as px
import pandas as pd
import argparse
import json

def splitList(interval,myList,exon_intron,header,datadict):
    for x in range(interval):
        new_list = myList[x::interval]
        label= exon_intron + " "+str(x+1)
        datadict[label] =new_list
        header.append(label)
    return datadict,header

def parseArgs():
    parser = argparse.ArgumentParser(description='Extract the individual Ramp sequences from a collection of genes')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input file with two lines: exon=[list of probabilities] and intron=[list of probabilities]')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output file path. Example: APOE_24.html')
    parser.add_argument('-n', '--interval', type=int, default=24, help='Interval size')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print the split lists for each interval')
    parser.add_argument('-z', '--zscore', action="store_true", help='Optional normalization with z score')
    parser.add_argument('-t', '--graph_title', default='Predictions at APOE position 850 Using Different Contexts and Window Start Sites',type=str,help='Change Graph Title')
    parser.add_argument('-oe', '--only_exon', action="store_true", help='Graph only the exon values. Note: Input file must still include the intron list')
    parser.add_argument('-oi', '--only_intron', action="store_true", help='Graph only the intron values. Note: Input file must still include the exon list')

    args = parser.parse_args()
    return args


def makeHTML(inputF,outputF,interval,z_score,graph_title,only_exon,only_intron):
    assert not (only_intron and only_exon), "Must specify exon or introns to be displayed. Cannot use both -oe and -oi"
    header=[]
    datadict = dict()
    with open(inputF) as input_file:
        line1=input_file.readline().strip()
        assert line1.startswith("exon=["),'First line must start with exon=['
        exon= json.loads(line1[5:])
        line2=input_file.readline().strip()
        assert line2.startswith("intron=["),'Second line must start with intron=['
        intron= json.loads(line2[7:])

    if not only_intron:
        exon_intron="Exon"
        datadict,header=splitList(interval,exon,exon_intron,header,datadict)
    if not only_exon:
        exon_intron="Intron"
        datadict,header=splitList(interval,intron,exon_intron,header,datadict)
    df= pd.DataFrame(datadict)
    #Optional Z score normalization
    if z_score:
        from scipy.stats import zscore
        df = df.apply(zscore)
    
    assert len(exon) == len(intron),'Length of the Exon list must be the same as the intron list' #Ensures that the intron and exon lists are the same length and can be plotted on the same graph
    assert len(exon) % interval == 0, 'Length of the exon/intron lists must be divisible by the interval.' #Ensure that the number of predictions is divisible by the interval. That makes sure that each split list has the same number of items
    
    #Create the figure
    fig = px.line(df,x=range(1,len(exon),interval), y=header, log_y=False)
    fig.update_layout(
        title=graph_title,
        title_x=0.5,
        xaxis_title='Context of Predicted Nucleotide',
        yaxis_title='Probability',
        legend_title='Feature'
    )
    
    #fig.write_html("APOE_predictions_token_" +str(interval) + "_pos_850_zscore.html")
    fig.write_html(outputF)
    return datadict


def verbosePrint(datadict):
    print("Verbose: interval dictionaries")
    for key,value in datadict.items():
        print(key, value[int(len(value)/2)])
    print()


def main():
    args=parseArgs()
    datadict = makeHTML(args.input,args.output,args.interval,args.zscore,args.graph_title, args.only_exon, args.only_intron)
    if args.verbose:
        verbosePrint(datadict)

if __name__ == '__main__':
    main()
