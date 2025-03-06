import plotly.express as px
import plotly.graph_objs as go
import pandas as pd
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description='Graph how input sequence length impacts exon predictions or auc')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input tsv. Column headers are required as shown in examples')
    parser.add_argument('-o', '--output_prefix', type=str, required=True,help='(output) Output prefix.')
    parser.add_argument('-g', '--graph', type=str.lower, required=True,choices=['exon','auc'],help='Choose either the exon or auc graph')
    args = parser.parse_args()
    return args

def line(error_y_mode=None, **kwargs):
    '''
    Credit to https://stackoverflow.com/a/69594497 for most of this function
    '''
    """Extension of `plotly.express.line` to use error bands."""
    ERROR_MODES = {'bar','band','bars','bands',None}
    if error_y_mode not in ERROR_MODES:
        raise ValueError(f"'error_y_mode' must be one of {ERROR_MODES}, received {repr(error_y_mode)}.")
    if error_y_mode in {'bar','bars',None}:
        fig = px.line(**kwargs)
    elif error_y_mode in {'band','bands'}:
        if 'error_y' not in kwargs:
            raise ValueError(f"If you provide argument 'error_y_mode' you must also provide 'error_y'.")
        figure_with_error_bars = px.line(**kwargs)
        fig = px.line(**{arg: val for arg,val in kwargs.items() if arg != 'error_y'})
        #fig['data'][0]['line']['color']='black'
        #.update_traces(selector={'Gene':'VEGFA'}, line_color='black')
        for data in figure_with_error_bars.data:
            x = list(data['x'])
            y_upper = list(data['y'] + data['error_y']['array'])
            y_lower = list(data['y'] - data['error_y']['array'] if data['error_y']['arrayminus'] is None else data['y'] - data['error_y']['arrayminus'])
            color = f"rgba({tuple(int(data['line']['color'].lstrip('#')[i:i+2], 16) for i in (0, 2, 4))},.3)".replace('((','(').replace('),',',').replace(' ','')
            fig.add_trace(
                go.Scatter(
                    x = x+x[::-1],
                    y = y_upper+y_lower[::-1],
                    fill = 'toself',
                    fillcolor = color,
                    line = dict(
                        color = 'rgba(255,255,255,0)'
                        #color = 'rgba(0,0,0,0)'
                    ),
                    hoverinfo = "skip",
                    showlegend = False,
                    legendgroup = data['legendgroup'],
                    xaxis = data['xaxis'],
                    yaxis = data['yaxis'],
                )
            )
        # Reorder data as said here: https://stackoverflow.com/a/66854398/8849755
        reordered_data = []
        for i in range(int(len(fig.data)/2)):
            reordered_data.append(fig.data[i+int(len(fig.data)/2)])
            reordered_data.append(fig.data[i])
        fig.data = tuple(reordered_data)
    x_pos=[24,48,96,192,384,768,1536,3072,6144,12288,24576]
    #title="Input Sequence Length Impact on AUC",
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',  # Set the paper background color to transparent
        plot_bgcolor='rgba(0,0,0,0)',
        title_font_color='black',
        xaxis_type="log",
        xaxis_exponentformat="power",
        title_x=0.5,
        #xaxis_title="Input Sequence Length",
        legend_title='Gene',
        xaxis=dict(
            tickformat=",",  # Use comma as thousands separator
            tickvals=x_pos,
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
            #title="AUC",
            title_font_color="black",
            tickfont=dict(color="black"),
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks="outside"  # Add tick marks
        ),
        legend=dict(font=dict(color='black'))
    )
    return fig


def main():
    args=parseArgs()
    df=pd.read_csv(args.input, sep="\t")
    if args.graph == 'exon':
        #use this code for exon probability graphs
        fig = line(
            data_frame = df,
            x = 'Input',
            y = 'Exon_Probability',
            error_y = 'SD',
            error_y_mode = 'band', # Here you say `band` or `bar`.
            color = 'Gene',
            title = f'Input Sequence Length Impact on Exon Probabilities',
            markers = '.',
        )
        fig.update_layout(
            xaxis=dict(
                title="Sequence Length",
            ),
            yaxis=dict(
                title="Exon Probability"))
        fig.write_html(f"{args.output_prefix}.html")
    elif args.graph == 'auc':
        #use this code for auc graphs
        for normalization_method in ['raw','zscore','subtraction']:
            df_filtered=df[df['Normalized'] == normalization_method]

            for position in ['First','Middle','Last']:
                ##print(df)
                fig = line(
                    data_frame = df_filtered,
                    x = 'Input',
                    y = position,
                    error_y = 'SD',
                    error_y_mode = 'band', # Here you say `band` or `bar`.
                    color = 'Gene',
                    title = f'Input Sequence Length Impact on AUC',
                    markers = '.',
                )
                fig.update_layout(
                    xaxis=dict(
                        title="Sequence Length",
                    ),
                    yaxis=dict(
                        title="AUC"))
                fig.write_html(f"{args.output_prefix}_{normalization_method}_{position}.html")

if __name__ == '__main__':
    main()


