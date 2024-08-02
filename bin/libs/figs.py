from scipy.ndimage.filters import gaussian_filter1d
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np
import scipy.stats as stats

interpolation_point_count_dict = {
    'read_length_distribution': (None, 10000, 3),
    'yield_plot': (None, 10000, 3),
    'phred_score_density': (None, 1000, 3),
    'over_time_graph': (None, 1000, 3),
    'scatterplot': (10000, 4000, 3),
    'phred_violin': (10000, 4000, 3),
}

toulligqc_colors = {'all': '#fca311',  # Yellow
                    'all_1d2': '#fca311',  # Yellow
                    'pass': '#51a96d',  # Green
                    'fail': '#d90429',  # Red
                    'barcode_pass': '#79bf90',  # Green
                    'barcode_fail': '#fb1941',  # Red
                    'sequence_length_over_time': '#205b47',
                    'phred_score_over_time': '#7aaceb',
                    'speed_over_time': '#AE3F7B',
                    'nseq_over_time': '#edb773',
                    'pie_chart_palette': ["#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#786fa6", "#f8a5c2",
                                          "#63cdda", "#ea8685", "#596275"],
                    'green_zone_color': 'rgba(0,100,0,.1)'
                    }

def make_plots(df_stats, positions, q1, q_scores, q3):
    fig = make_subplots(rows=3, cols=1)
    fig.add_trace(
        over_time_graph(time_series=positions,
                     percentile_25_series=q1,
                     percentile_50_series=q_scores,
                     percentile_75_series=q3,
                     result_directory='your/directory/path',
                     graph_name='Qscore over position',
                     yaxis_title='Q Score',
                     log=False,
                     sigma=1,
                     yaxis_starts_zero=False,
                     green_zone_starts_at=None),
        row=3, col=1)
    
    fig.add_trace(
        ATGC_graph(time_series=positions,
                     percentage_G=df_stats["%G"],
                     percentage_C=df_stats["%C"],
                     percentage_A=df_stats["%A"],
                     percentage_T=df_stats["%T"],
                     result_directory='your/directory/path',
                     graph_name='Base % over sequence',
                     yaxis_title='Q Score',
                     log=False,
                     sigma=1,
                     yaxis_starts_zero=False,
                     green_zone_starts_at=None),
        row=3, col=1)
                  
        
    fig.add_trace(
        GC_content_plot(df_stats['GC_content']),
        row=3, col=1)

    fig.update_layout(height=600, width=800, title_text="Side By Side Subplots")
    fig.write_html("cc.html")


def over_time_graph(time_series,
                     percentile_25_series,
                     percentile_50_series,
                     percentile_75_series,
                     result_directory,
                     graph_name,
                     yaxis_title,
                     color=toulligqc_colors['phred_score_over_time'],
                     log=False,
                     sigma=1,
                     yaxis_starts_zero=False,
                     green_zone_starts_at=None,
                     green_zone_color=toulligqc_colors['phred_score_over_time']):

    # Apply Gaussian filter to the percentile series
    filtered_25 = gaussian_filter1d(percentile_25_series[0], sigma=sigma)
    filtered_50 = gaussian_filter1d(percentile_50_series[0], sigma=sigma)
    filtered_75 = gaussian_filter1d(percentile_75_series[0], sigma=sigma)

    # Create a Plotly figure
    fig = go.Figure()

    # Add the 25th percentile trace
    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_25, 
        name="Lower percentile",
        mode='lines',
        line=dict(color=color, width=2),
        fill=None
    ))

    # Add the 75th percentile trace
    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_75, 
        name="Uper percentile",
        mode='lines',
        line=dict(color=color, width=2),
        fill='tonexty', # Fill the area between this line and the previous line
        fillcolor='rgba(0, 100, 0, 0.2)' # Semi-transparent fill
    ))

    # Add the 50th percentile (median) trace
    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_50, 
        name="Mean",
        mode='lines',
        line=dict(color='black', width=2)
    ))

    # Define the green zone if required
    if green_zone_starts_at is not None:
        min_x = min(time_series[0])
        max_x = max(time_series[0])
        max_y = max(filtered_75) * 1.05
        fig.add_trace(go.Scatter(
            mode="lines",
            x=[min_x, max_x],
            y=[max_y, max_y],
            line=dict(width=0),
            hoverinfo="skip",
            showlegend=False,
        ))
        fig.add_trace(go.Scatter(
            mode="lines",
            name="Green Zone",
            x=[min_x, max_x],
            y=[green_zone_starts_at, green_zone_starts_at],
            fill='tonexty',
            fillcolor=green_zone_color,
            line=dict(width=0),
            hoverinfo="skip",
        ))
    
    ######
    ### add REVERSE
    #####
        
    # Apply Gaussian filter to the percentile series
    filtered_25 = gaussian_filter1d(percentile_25_series[1], sigma=sigma)
    filtered_50 = gaussian_filter1d(percentile_50_series[1], sigma=sigma)
    filtered_75 = gaussian_filter1d(percentile_75_series[1], sigma=sigma)
                         
    # Add the 25th percentile trace
    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_25, 
        name="Lower percentile",
        mode='lines',
        line=dict(color=color, width=2),
        fill=None,
        visible=False
    ))

    # Add the 75th percentile trace
    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_75, 
        name="Uper percentile",
        mode='lines',
        line=dict(color=color, width=2),
        fill='tonexty', # Fill the area between this line and the previous line
        fillcolor='rgba(0, 100, 0, 0.2)',
        visible=False
    ))

    # Add the 50th percentile (median) trace
    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_50, 
        name="Mean",
        mode='lines',
        line=dict(color='black', width=2),
        visible=False
    ))

    # Define the green zone if required
    if green_zone_starts_at is not None:
        min_x = min(time_series[1])
        max_x = max(time_series[1])
        max_y = max(filtered_75) * 1.05
        fig.add_trace(go.Scatter(
            mode="lines",
            x=[min_x, max_x],
            y=[max_y, max_y],
            line=dict(width=0),
            hoverinfo="skip",
            showlegend=False,
            visible=False
        ))
        fig.add_trace(go.Scatter(
            mode="lines",
            name="Green Zone",
            x=[min_x, max_x],
            y=[green_zone_starts_at, green_zone_starts_at],
            fill='tonexty',
            fillcolor=green_zone_color,
            line=dict(width=0),
            hoverinfo="skip",
            visible=False
        ))
    
    ########
    # Add buttons
    ###############
                         
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list([
                    dict(
                        args=[{'visible': [True, True, True, False, False, False, ]},
                              {**_xaxis('Q score', dict(visible=True)),
                               **_yaxis('Position', dict(visible=True)),
                               'plot_bgcolor': '#e5ecf6'}],
                        label="Begining ",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False, True, True, True, ]},
                              {**_xaxis('Q score', dict(visible=True)),
                               **_yaxis('Position', dict(visible=True)),
                               'plot_bgcolor': '#e5ecf6'}],
                        label="End",
                        method="update"
                    )
                ]),
                pad={"r": 20, "t": 160, "l": 40, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            )
        ]
    )
                         
    y_axis_range_mode = 'tozero' if yaxis_starts_zero else 'normal'

    fig.update_layout(
        title=graph_name,
        autosize=True,
        xaxis_title='<b> Position</b> ',
        yaxis_title='<b>' +yaxis_title +'</b>',
        hovermode='x unified',
        yaxis=dict(rangemode=y_axis_range_mode)
    )


    if log:
        fig.update_yaxes(type="log")

    return fig 


def ATGC_graph(time_series,
                     percentage_G,
                     percentage_C,
                     percentage_A,
                     percentage_T,
                     result_directory,
                     graph_name,
                     yaxis_title,
                     color=toulligqc_colors['phred_score_over_time'],
                     log=False,
                     sigma=1,
                     yaxis_starts_zero=False,
                     green_zone_starts_at=None,
                     green_zone_color=toulligqc_colors['phred_score_over_time']):

    filtered_G = gaussian_filter1d(percentage_G[0], sigma=sigma)
    filtered_C = gaussian_filter1d(percentage_C[0], sigma=sigma)
    filtered_A = gaussian_filter1d(percentage_A[0], sigma=sigma)
    filtered_T = gaussian_filter1d(percentage_T[0], sigma=sigma)


    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_G, 
        name="Base G Percentage",
        mode='lines',
        line=dict(color="red", width=2),
        fill=None
    ))

    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_A, 
        name="Base A Percentage",
        mode='lines',
        line=dict(color="green", width=2),
    ))


    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_T, 
        name="Base T Percentage",
        mode='lines',
        line=dict(color="blue", width=2),

    ))

    fig.add_trace(go.Scatter(
        x=time_series[0], 
        y=filtered_C, 
        name="Base C Percentage",
        mode='lines',
        line=dict(color='black', width=2)
    ))

    ######
    ### add REVERSE
    #####
    filtered_G = gaussian_filter1d(percentage_G[1], sigma=sigma)
    filtered_C = gaussian_filter1d(percentage_C[1], sigma=sigma)
    filtered_A = gaussian_filter1d(percentage_A[1], sigma=sigma)
    filtered_T = gaussian_filter1d(percentage_T[1], sigma=sigma)                   

    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_G, 
        name="Base G Percentage",
        mode='lines',
        line=dict(color="red", width=2),
        fill=None,
        visible=False
    ))

    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_A, 
        name="Base A Percentage",
        mode='lines',
        line=dict(color="green", width=2),
        visible=False

    ))

    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_T, 
        name="Base T Percentage",
        mode='lines',
        line=dict(color="blue", width=2),
        visible=False
    ))


    fig.add_trace(go.Scatter(
        x=time_series[1].multiply(other = -1), 
        y=filtered_C, 
        name="Base C Percentage",
        mode='lines',
        line=dict(color='black', width=2),
        visible=False
    ))
                         
    ########
    # Add buttons
    ###############
                         
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list([
                    dict(
                        args=[{'visible': [True, True, True, True, False, False, False, False]},
                              {**_xaxis('Base percent', dict(visible=True)),
                               **_yaxis('Position', dict(visible=True)),
                               'plot_bgcolor': '#e5ecf6'}],
                        label="Begining ",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False, False, True, True, True, True]},
                              {**_xaxis('Base percent', dict(visible=True)),
                               **_yaxis('Position', dict(visible=True)),
                               'plot_bgcolor': '#e5ecf6'}],
                        label="End",
                        method="update"
                    )
                ]),
                pad={"r": 20, "t": 160, "l": 40, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            )
        ]
    )

    y_axis_range_mode = 'tozero' if yaxis_starts_zero else 'normal'


    fig.update_layout(
        title=graph_name,
        autosize=True,
        xaxis_title='<b> Position </b>',
        yaxis_title= '<b>' + yaxis_title + '</b>',
        hovermode='x unified',
        yaxis=dict(rangemode=y_axis_range_mode)
    )

    if log:
        fig.update_yaxes(type="log")
    return fig 


def GC_content_plot(gc_content_percentages):
    max_gc = 100  
    gc_distribution = [0] * (max_gc + 1)
    for gc_content in gc_content_percentages:
        gc_distribution[int(gc_content)] += 1
    
    real_distribution_smooth = np.convolve(gc_distribution, np.ones(5)/5, mode='same')
    
    mean_gc = np.mean(gc_content_percentages)
    std_gc = np.std(gc_content_percentages)
    theoretical_distribution = [stats.norm.pdf(gc, mean_gc, std_gc) * len(gc_content_percentages) for gc in range(max_gc + 1)]
    

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(max_gc + 1)), y=real_distribution_smooth, mode='lines', name='Distribution GC réelle'))
    fig.add_trace(go.Scatter(x=list(range(max_gc + 1)), y=theoretical_distribution, mode='lines', name='Distribution GC théorique'))
    
    fig.update_layout(title='GC content distribution',
                      autosize=True,
                      #width=1400,
                      #height=400,
                      xaxis_title='<b> GC content (%) </b>',
                      yaxis_title='<b> Count </b>',
                      legend_title='<b> Légende </b>')
    

    return fig


def interpolation_points(series, graph_name):
    count = len(series)
    threshold, npoints, sigma = interpolation_point_count_dict[graph_name]

    if threshold is not None:
        if count > threshold:
            result = npoints
        else:
            result = count
    else:
        result = npoints

    return result, sigma


def _xaxis(title, args=None):
    axis_dict = dict(
        title='<b>' + title + '</b>',
        titlefont_size=14,
        tickfont_size=12)

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(xaxis=axis_dict)

def _yaxis(title, args=None):
    axis_dict = dict(
        title='<b>' + title + '</b>',
        titlefont_size=14,
        tickfont_size=12,
        fixedrange=True)

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(yaxis=axis_dict)