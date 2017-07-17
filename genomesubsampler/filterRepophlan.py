# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Program to filter repophlan results by the quality scores.

Example usage:

    filter_repophlan.py -r ./repophlan_microbes_wscores.txt -s 0.6

Will produce two files:

    repophlan_filepaths.good and repophlan_filepaths.bad

each will contain the filenames in the following columns:

    faa_lname,ffn_lname,fna_lname,frn_lname

as well as the original score values and the summarized score.

Score can be calculated as in the repophlan paper, using a normalized and
averaged score, or by passing the -a/--avg flag, to just use an average of
each of the four scores per genome.
"""

import pandas as pd
import click


def calc_norm_score(genome_df, cols=['score_faa',
                                     'score_fna', 'score_rrna', 'score_trna']):
    sub_df = genome_df[cols].copy()
    sub_df_norm = ((sub_df - sub_df.mean()) / sub_df.std()).mean(axis=1)
    score_avg = (sub_df_norm - sub_df_norm.min()
                 ) / (sub_df_norm.max() - sub_df_norm.min())

    return(score_avg)


def calc_avg_score(genome_df, cols=['score_faa',
                                    'score_fna', 'score_rrna', 'score_trna']):
    sub_df = genome_df[cols].copy()
    score_avg = sub_df.mean(axis=1)

    return(score_avg)


@click.command()
@click.option('-r', '--repophlan_fp', type=str,
              help='path to repophlan results')
@click.option('-s', '--score_threshold',
              type=float, default=0.8, show_default=True,
              help='score threshold to include in good file \
                   (default: %(default)s)')
@click.option('-o', '--output_fp',
              type=str, default='./repophlan_filepaths', show_default=True,
              help='path for output (will output [output_fp].good and \
                 [output_fp].bad; default: %(default)s)')
@click.option('-c', '--score_cols',
              type=str, default='score_faa,score_fna,score_rrna,score_trna',
              show_default=True, help='comma delimited list of \
              columns to use for score calculation (default: %(default)s)')
@click.option('-f', '--fp_cols',
              type=str, default='faa_lname,ffn_lname,fna_lname,frn_lname',
              show_default=True, help='comma delimited list of columns \
              to include in output (default: %(default)s)')
@click.option('-a', '--avg', is_flag=True,
              help='use simple average rather than normalized \
              average for score')
def _main(repophlan_fp, score_threshold, output_fp,
          score_cols, fp_cols, avg):

    input_file = repophlan_fp
    score_cols = score_cols.strip().split(',')
    fp_cols = fp_cols.strip().split(',')

    genome_df = pd.read_csv(input_file, sep='\t', header=0, index_col=0)

    if avg:
        comb_col = 'score_avg'
        genome_df[comb_col] = calc_avg_score(genome_df, cols=score_cols)
    else:
        comb_col = 'score_norm'
        genome_df[comb_col] = calc_norm_score(genome_df, cols=score_cols)

    good_fp = output_fp + '.good'
    bad_fp = output_fp + '.bad'

    genome_df.loc[genome_df[comb_col] >= score_threshold, fp_cols +
                  score_cols + [comb_col]].to_csv(good_fp, sep='\t')
    genome_df.loc[genome_df[comb_col] < score_threshold, fp_cols +
                  score_cols + [comb_col]].to_csv(bad_fp, sep='\t')


if __name__ == "__main__":
    _main()
