import pandas as pd
import argparse

from libs import qc_generator
from libs.figs import over_time_graph
from libs.figs import ATGC_graph
from libs.figs import GC_content_plot
from libs.report import make_report

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-q', '--fastq', action='append', dest='fastq', required=True,
                          help='FASTQ file (necessary if no sequencing summary file), ' +
                               'can also be in a tar.gz archive')
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
    _stats, quality_average, length = toullig.process_fastq_file(fastq, positions,  reverse=reverse)
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

    make_report(Q_over_time, atgc, gc, args.project, args.positions)

if __name__ == "__main__":
    main()
    
    