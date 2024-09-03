import pandas as pd
import argparse

from QC_tools import qc_generator
from QC_tools.figs import over_time_graph
from QC_tools.figs import ATGC_graph
from QC_tools.figs import GC_content_plot
from QC_tools.reporte_maker import reporter
import plotly.graph_objects as go

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-q', '--fastq', action='append', dest='fastq', required=True,
                          help='FASTQ file (necessary if no sequencing summary file), ' +
                               'can also be in a tar.gz archive')
parser.add_argument('-r', "--conf_params", action='append', dest="conf_params", help="config parameters")
parser.add_argument('-w', "--work_params", action='append', dest="work_params", help="workflow parameters")
parser.add_argument('-e', "--env", action='append', dest="env", help="packages versions")

parser.add_argument('-s', "--stats", action='store', dest="stats", help="template maker stats")
parser.add_argument('-p', "--positions", action='store', dest="positions", help="Number of base position to show", type=int, default=1000)
parser.add_argument('-t', "--thread", action='store', dest="thread", help="Number of threads", type=int, default=2)
parser.add_argument('-b', "--batch-size", action='store', dest="batch_size", help="Batch size", type=int, default=500)
parser.add_argument('-n', "--project", action='store', dest="project", help="project name", type=str, default="test_project")
args = parser.parse_args()

columns = ["pos",
           "Q_score",
           "GC_content", 
           "%A", 
           "%T",
           "%G", 
           "%C", 
           "percentile_25", 
           "percentile_75"]


def compQC(fastq, positions, reverse=False):
    _stats, quality_average, length = qc_generator.process_fastq_file(fastq, positions,  reverse=reverse)
    df_stats = pd.DataFrame(_stats, columns=columns)
    df_stats['GC_percentage'] =  df_stats["%G"]+ df_stats["%C"]
    
    df_qc = pd.DataFrame({"qscore":quality_average, 
                          "length":length})
    df_qc['passes_filtering'] = ['pass' if x > 9 else 'fail' for x in df_qc.qscore]

    q_scores = df_stats.Q_score
    positions = df_stats.pos
    q1 = df_stats.percentile_25
    q3 = df_stats.percentile_75
    
    return df_stats, df_qc, positions, q_scores, q1, q3


def pie_unknown(df):
    labels = ['Simulated feature counts','Unknow feature counts']
    values = [int(df['Simulated UMI counts']), 
              int(df['Unknown transcript counts'])]

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, pull=[0, 0, 0.2, 0])])
    return fig


def create_table(stats_data, table_width=500, table_height=400, margin=10):
    data=[go.Table(
        header=dict(
            values=['', '<b>Value</b>'],
            fill_color='#EBEBEB',
            font=dict(color='black', size=12)
        ),
        cells=dict(
            values=[
                ['<b>Simulated Cell BC</b>', '<b>Simulated Filtered-Out</b>', '<b>Simulated UMI counts</b>', '<b>Mean UMI per cell</b>'],
                [
                    stats_data['Simulated Cell BC'], 
                    stats_data['Simulated Filtered-Out'], 
                    stats_data['Simulated UMI counts'],
                    int(int(stats_data['Simulated UMI counts']) / int(stats_data['Simulated Cell BC']))
                ]
            ],
            fill_color='#f4f4f4',
            font=dict(color='black', size=11)
        )
    )]

    fig = go.Figure(data)

    fig.update_layout(
        #width=table_width,
        #height=table_height,
        margin=dict(l=margin, 
                    r=margin, 
                    t=20, 
                    b=0),
    )
    return fig


def main():
    df_stats, df_qc, pos, q_scores, q1, q3 = compQC(args.fastq, args.positions)
    df_stats_rev, df_qc_rev, pos_rev, q_scores_rev, q1_rev, q3_rev = compQC(args.fastq, args.positions, reverse=True)

    Q_over_time = over_time_graph(time_series=[pos, pos_rev],
                     percentile_25_series=[q1, q1_rev],
                     percentile_50_series=[q_scores, q_scores_rev],
                     percentile_75_series=[q3, q3_rev],
                     result_directory='your/directory/path',
                     graph_name='Qscore over position',
                     yaxis_title='Q Score',
                     log=False,
                     sigma=1,
                     yaxis_starts_zero=False,
                     green_zone_starts_at=None)

    atgc = ATGC_graph(time_series=[pos, pos_rev],
                     percentage_G=[df_stats["%G"], df_stats_rev["%G"]],
                     percentage_C=[df_stats["%C"], df_stats_rev["%C"]],
                     percentage_A=[df_stats["%A"], df_stats_rev["%A"]],
                     percentage_T=[df_stats["%T"], df_stats_rev["%T"]],
                     result_directory='your/directory/path',
                     graph_name='Base % over sequence',
                     yaxis_title='Q Score',
                     log=False,
                     sigma=1,
                     yaxis_starts_zero=False,
                     green_zone_starts_at=None)

    gc = GC_content_plot(df_stats['GC_percentage'])

    df_stats = pd.read_csv(args.stats)
    pie_fig = pie_unknown(df_stats)
    stats_table = create_table(df_stats)

    reporter(stats_table, 
             pie_fig, 
             Q_over_time, 
             atgc, 
             gc, 
             args.conf_params, 
             args.work_params, 
             args.project,
             args.env)


if __name__ == "__main__":
    main()
    
    