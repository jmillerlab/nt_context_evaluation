import plotly.express as px
import pandas as pd
import argparse
import json


def parseArgs():
    parser = argparse.ArgumentParser(description='Graph the exon/intron probabilities across APOE when predicting the first, middle, or last nucleotide')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input json file containing a dictionary with key:value pairs of probabilities.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output file path. Example: APOE_first_middle_last_graph.html. This will be appended to the name of the feature. For example, if the first feature is "Exon" the file path to the graph with that first feature will be "Exon_APOE_first_middle_last_graph.html"')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print the split lists for each interval')
    parser.add_argument('-z', '--zscore', action="store_true", help='Optional normalization with z score')
    parser.add_argument('-t', '--graph_title', default='Predictions at APOE position 850 Using Different Contexts and Window Start Sites',type=str,help='Change Graph Title')

    args = parser.parse_args()
    return args

def plotVals(args):
    with open(args.input) as inputF:
        allValues = json.loads(inputF.read())
    
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
        fig = px.line(df,x=range(1,len(fml_dict['first'])+1), y=['first','middle','last'], log_y=False)
        
        #Update layout to make it look prettier
        fig.update_layout(
            title=f'APOE {feature} Predictions from Different Contexts',
            title_x=0.5,
            xaxis_title='APOE Gene Sequence',
            yaxis_title=f'Probability of {feature}',
            legend_title='Position of Predicted Nucleotide'
        )
        if args.zscore:
            fig.update_layout(
                    yaxis_title='Probability of Exon Z score'
            )
        # Add Exon 1
        fig.add_shape(type='line',
                      x0=1,
                      y0=1.0,
                      x1=60,
                      y1=1.0,
                      line=dict(color='Black',),
                      name="1-60",
                      xref='x',
                      yref='y')
        
        fig.add_annotation(
            x=25, y=1.02,
            text="Exon 1",
            showarrow=False,
            font=dict(size=14)
        )
        
        # Add Exon 2
        fig.add_shape(type='line',
                      x0=821,
                      y0=1.0,
                      x1=886,
                      y1=1.0,
                      line=dict(color='Black',),
                      name="821-886",
                      xref='x',
                      yref='y')
        
        fig.add_annotation(
            x=845, y=1.02,
            text="Exon 2",
            showarrow=False,
            font=dict(size=14)
        )
        # Add Exon 3
        fig.add_shape(type='line',
                      x0=1979,
                      y0=1.0,
                      x1=2171,
                      y1=1.0,
                      line=dict(color='Black',),
                      name="1979-2171",
                      xref='x',
                      yref='y')
        
        fig.add_annotation(
            x=2050, y=1.02,
            text="Exon 3",
            showarrow=False,
            font=dict(size=14)
        )

        # Add Exon 4
        fig.add_shape(type='line',
                      x0=2752,
                      y0=1.0,
                      x1=3612,
                      y1=1.0,
                      line=dict(color='Black',),
                      name="2752-3612",
                      xref='x',
                      yref='y')
        
        fig.add_annotation(
            x=3150, y=1.02,
            text="Exon 4",
            showarrow=False,
            font=dict(size=14)
        )
        fig.write_html(feature + "_" +args.output)
        print(f"Wrote {feature}_{args.output}") #tells the user where the output files were written

def main():
    args=parseArgs()
    plotVals(args)
    #if args.verbose:
    #    verbosePrint(datadict)

if __name__ == '__main__':
    main()

