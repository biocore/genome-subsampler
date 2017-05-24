#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

#
# parser for RepoPhlAn-downloaded genomes
#

import click
import pandas as pd


def parse_repophlan(repophlan_wscores_fp):
    """Compute basic statistics of RepoPhlAn-downloaded genomes.

    Parameters
    ----------
    repophlan_wscores_fp : str
        File path to RepoPhlAn summary table with scores.

    Returns
    -------
    list of str
        Human-readable report of basic statistics of genomes.
    """
    df = pd.read_table(repophlan_wscores_fp, index_col=0, header=0)
    out = []
    out.append('Total number of genomes: %s.' % df.shape[0])
    out.append('Number of RefSeq genomes: %s.'
               % df['assembly_accession'].str.contains('GCF_').sum())
    out.append('With genome sequences (fna): %s.' % df['fna_lname'].count())
    out.append('With protein sequences (faa): %s.' % df['faa_lname'].count())
    out.append('With protein-coding DNA sequences (ffn): %s.'
               % df['ffn_lname'].count())
    out.append('With RNA-coding DNA sequences (frn): %s.'
               % df['frn_lname'].count())
    return out


@click.command()
@click.option('--repophlan-wscores-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='RepoPhlAn summary table with scores')
def _main(repophlan_wscores_fp):
    """Parser for RepoPhlAn-downloaded genomes.
    """
    out = parse_repophlan(repophlan_wscores_fp)
    click.echo('\n'.join(out))
    click.echo('Task completed.')


if __name__ == "__main__":
    _main()
