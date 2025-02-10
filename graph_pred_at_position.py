import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import argparse
import json
import plotly.offline as pio

def splitList(fig,interval,myList,exon_intron,header,datadict,markers,start,end):
    for x in range(interval):
        new_list = myList[x::interval]
        label= exon_intron + " "+str(x+1)
        datadict[label] =new_list
        header.append(label)
        x_vals=range(start,end,interval)
        if markers:
            fig.add_trace(go.Scatter(x=[x_vals[a]+x for a in range(len(x_vals))],y=new_list,mode='lines+markers',name=label))
        else:
            fig.add_trace(go.Scatter(x=[x_vals[a]+x for a in range(len(x_vals))],y=new_list,mode='lines',name=label))
    return fig,datadict,header

def parseArgs():
    parser = argparse.ArgumentParser(description='Graph the exon/intron probabilities for a single nucleotide using different contexts')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input file with two lines: exon=[list of probabilities] and intron=[list of probabilities]')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output file path. Example: APOE_24.html')
    parser.add_argument('-n', '--interval', type=int, default=24, help='Interval size')
    parser.add_argument('-p', '--position', type=int, default=850, help='Position of nucleotide in gene')
    parser.add_argument('-s', '--start', type=int, default=0, help='Start of visualization window')
    parser.add_argument('-e', '--end', type=int, default=-1, help='End of visualization window')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print the split lists for each interval')
    parser.add_argument('-z', '--zscore', action="store_true", help='Optional normalization with z score')
    parser.add_argument('-g', '--gene', type=str.lower, required=False,choices=['apoe','egfr','tnf','tp53','vegfa','control'],help='gene used in graph')
    parser.add_argument('-t', '--graph_title', default='APOE Predictions at position 850 using Different Contexts and Window Start Sites',type=str,help='Change Graph Title')
    parser.add_argument('-oe', '--only_exon', action="store_true", help='Graph only the exon values. Note: Input file must still include the intron list')
    parser.add_argument('-oi', '--only_intron', action="store_true", help='Graph only the intron values. Note: Input file must still include the exon list')
    parser.add_argument('-m', '--markers', action="store_true", help='Draw the plot with the point markers visible')
    args = parser.parse_args()
    genePos={'apoe':850,'egfr':183535,'tnf':180,'tp53':6800,'vegfa':1500,'control':850}
    if args.gene is not None:
        args.position=genePos[args.gene]
        args.graph_title=args.graph_title.replace('850',str(args.position)) 
        args.graph_title=args.graph_title.replace('APOE',args.gene.upper()) 
    return args

def formatFig(fig,end,graph_title):
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',  # Set the paper background color to transparent
        plot_bgcolor='rgba(0,0,0,0)',
        title=dict(text=graph_title,font=dict(color="black")),
        title_x=0.5,
        legend_title='Feature'
    )
    fig.update_xaxes(title="Context of Predicted Nucleotide",title_font_color="black",tickfont=dict(color="black"),showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(title="Probability",title_font_color="black",tickfont=dict(color="black"),showline=True, linewidth=2, linecolor='black')
    fig.update_layout(legend=dict(font=dict(color='black')))
    #Add tick marks
    fig.update_xaxes(ticks="outside")
    fig.update_yaxes(ticks="outside")
    fig.update_layout(
        xaxis=dict(
            tickformat=","  # Use comma as thousands separator
        )
    )
    fig.update_layout(yaxis=dict(range=[0, 1.025])) #Makes y-axis start at 0. Reduces white space at bottom of graph
    fig.update_traces(line=dict(width=1)) #Makes line width thinner.
    #set x-axis to start at 0 and go to the end of the context size
    fig.update_layout(xaxis=dict(range=[0,end+1]))
    return fig

def makeHTML(inputF,outputF,interval,z_score,graph_title,only_exon,only_intron,markers,start,end):
    assert not (only_intron and only_exon), "Must specify exon or introns to be displayed. Cannot use both -oe and -oi"
    header=[]
    datadict = dict()
    fig=go.Figure()
    with open(inputF) as input_file:
        line1=input_file.readline().strip()
        assert line1.startswith("exon=["),'First line must start with exon=['
        exon= json.loads(line1[5:])
        if end==-1: 
            end=len(exon)
        exon=exon[start:end]
        line2=input_file.readline().strip()
        assert line2.startswith("intron=["),'Second line must start with intron=['
        intron= json.loads(line2[7:])
        intron=intron[start:end]

    if not only_intron:
        exon_intron="Exon"
        fig,datadict,header=splitList(fig,interval,exon,exon_intron,header,datadict,markers,start,end)
    if not only_exon:
        exon_intron="Intron"
        fig,datadict,header=splitList(fig,interval,intron,exon_intron,header,datadict,markers,start,end)
    df= pd.DataFrame(datadict)
    #Optional Z score normalization
    if z_score:
        from scipy.stats import zscore
        df = df.apply(zscore)
    
    assert len(exon) == len(intron),'Length of the Exon list must be the same as the intron list' #Ensures that the intron and exon lists are the same length and can be plotted on the same graph
    assert len(exon) % interval == 0, 'Length of the exon/intron lists must be divisible by the interval.' #Ensure that the number of predictions is divisible by the interval. That makes sure that each split list has the same number of items
    
    if z_score:
        fig.update_layout(
                yaxis_title='Probability Z score'
        )
    fig =formatFig(fig,end,graph_title) #Format Fig
    fig.write_html(outputF)
    pio.plot(fig,image_filename=outputF.replace(".html",""),image="svg")
    return datadict


def verbosePrint(datadict):
    print("Verbose: interval dictionaries")
    for key,value in datadict.items():
        print(key, value[int(len(value)/2)])
    print()


def main():
    args=parseArgs()
    datadict = makeHTML(args.input,args.output,args.interval,args.zscore,args.graph_title, args.only_exon, args.only_intron,args.markers,args.start,args.end)
    if args.verbose:
        verbosePrint(datadict)

if __name__ == '__main__':
    main()
