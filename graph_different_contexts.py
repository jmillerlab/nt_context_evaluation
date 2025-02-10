import plotly.express as px
#import plotly.graph_objects as go
import pandas as pd
import argparse
import json
import plotly.offline as pio


def parseArgs():
    parser = argparse.ArgumentParser(description='Graph the exon/intron probabilities across APOE when predicting the first, middle, or last nucleotide')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input json file containing a dictionary with key:value pairs of probabilities.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output file path. Example: APOE_first_middle_last_graph.html. This will be appended to the name of the feature. For example, if the first feature is "Exon" the file path to the graph with that first feature will be "Exon_APOE_first_middle_last_graph.html"')
    parser.add_argument('-g', '--gene', type=str.lower, required=True,choices=['apoe','egfr','tnf','tp53','vegfa','control'],help='gene used to draw exons on plot.')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print the split lists for each interval')
    parser.add_argument('-z', '--zscore', action="store_true", help='Optional normalization with z score')
    parser.add_argument('-s', '--svg', action="store_true", help='Optional Save SVG file with same name as html. Will save to Downloads directory')
    parser.add_argument('-t', '--graph_title', default=' Exon and Intron Predictions Using Difference Contexts',type=str,help='Change Graph Title')
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

def formatFig(fig,graph_title,x_axis_title,feature):
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
        yaxis_title=f'Probability of {feature}',
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
    
    fig.update_layout(yaxis=dict(range=[0, 1.025])) #Makes y-axis start at 0. Reduces white space at bottom of graph
    fig.update_traces(line=dict(width=0.5)) #Makes line width thinner.
    return fig

def plotVals(args):
    with open(args.input) as inputF:
        allValues = json.loads(inputF.read())

    drawExon={'apoe':apoe,'egfr':egfr,'tnf':tnf,'vegfa':vegfa,'control':control,'tp53':tp53}
    
    for feature,fml_dict in allValues.items(): #fml_dict = first_middle_last_dict
     
        df= pd.DataFrame(fml_dict)
        assert 'first' in fml_dict, 'input dictionary requires key word "first"'
        assert 'middle' in fml_dict, 'input dictionary requires key word "middle"'
        assert 'last' in fml_dict, 'input dictionary requires key word "last"'
        assert (len(fml_dict['first']) == len(fml_dict['middle'])) and (len(fml_dict['middle'])== len(fml_dict['last'])),"All lists must be the same length"
        
        #Optionally convert all probabilities to z scores
        if args.zscore:
            from scipy.stats import zscore
            df = df.apply(zscore)
        
        #Make the figure
        fig = px.line(df,x=range(1,len(fml_dict['first'])+1), y=['first','middle','last'], log_y=False,markers=args.markers)
        gene=args.gene.upper() 
        graph_title=f'{gene} {feature} Predictions from Different Contexts'
        x_axis_title=f'{gene} Gene Sequence'

        #Add Exon label to graph
        fig = drawExon[args.gene](fig)

        #Format Fig
        fig =formatFig(fig,graph_title,x_axis_title,feature)
        #Write to output files
        output_file=args.output + f"_{feature}"
        fig.write_html(output_file +".html")
        if args.svg:
            #Save svg file too. Will save to Downloads.
            pio.plot(fig,image_filename=output_file,image="svg")
        if args.verbose:
            print(f"Wrote {feature}_{args.output}") #tells the user where the output files were written

def main():
    args=parseArgs()
    plotVals(args)

if __name__ == '__main__':
    main()

