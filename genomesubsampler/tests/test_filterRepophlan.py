# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from unittest import TestCase, main
from click.testing import CliRunner
from tempfile import mkdtemp
from skbio.util import get_data_path
from shutil import rmtree
from genomesubsampler.filterRepophlan import (
    calc_norm_score, calc_avg_score, _main)


class filterRepophlan(TestCase):
    """Tests for filterRepophlan.py."""

    def setUp(self):
        """  Create working directory and test file
        """
        self.wkdir = mkdtemp()
        self.repophlan_fp = get_data_path('repophlan_microbes_wscores.txt')
        self.score_threshold = 0.8
        self.output_fp = './repoplan_filepaths'
        self.score_cols = 'score_faa,score_fna,score_rrna,score_trna'
        self.fp_cols = 'faa_lname,ffn_lname,fna_lname,frn_lname'
        self.avg = 'score_avg'

    def tearDown(self):
        """ Delete working directory and test files.
        """
        rmtree(self.wkdir)

    def test_calc_norm_score(self):
        """ test function calc_norm_score on simulated data
        """
        input_dict = {'score_faa': {'G000014725': 0.10000000000000001,
                                    'G000254175': 0.10000000000000001,
                                    'G000775715': 0.10000000000000001,
                                    'G000881615': 0.10000000000000001,
                                    'G000955195': 0.10000000000000001,
                                    'G001076295': 0.10000000000000001,
                                    'G001380675': 0.0},
                      'score_fna': {'G000014725': 1.0,
                                    'G000254175': 0.48599999999999999,
                                    'G000775715': 0.90200000000000002,
                                    'G000881615': 1.0,
                                    'G000955195': 1.0,
                                    'G001076295': 0.78500000000000003,
                                    'G001380675': 0.0},
                      'score_rrna': {'G000014725': 1.0,
                                     'G000254175': 1.0,
                                     'G000775715': 1.0,
                                     'G000881615': 0.0,
                                     'G000955195': 0.0,
                                     'G001076295': 0.90000000000000002,
                                     'G001380675': 0.0},
                      'score_trna': {'G000014725': 1.0,
                                     'G000254175': 0.80000000000000004,
                                     'G000775715': 0.90000000000000002,
                                     'G000881615': 0.0,
                                     'G000955195': 0.0,
                                     'G001076295': 0.80000000000000004,
                                     'G001380675': 0.0}}

        exp_norm_scores = {'G000881615': 0.568649417027694,
                           'G000014725': 1,
                           'G000775715': 0.94937522318669,
                           'G000254175': 0.807963889308086,
                           'G001380675': 0,
                           'G000955195': 0.568649417027694,
                           'G001076295': 0.872837563705059}

        input_df = pd.DataFrame(input_dict)
        obs_scores = calc_norm_score(input_df)
        for item in input_df.index:
            self.assertAlmostEqual(obs_scores[item],
                                   exp_norm_scores[item])

    def test_calc_avg_score(self):
        """ test function calc_avg_score on simulated data
        """
        input_dict = {'score_faa': {'G000014725': 0.10000000000000001,
                                    'G000254175': 0.10000000000000001,
                                    'G000775715': 0.10000000000000001,
                                    'G000881615': 0.10000000000000001,
                                    'G000955195': 0.10000000000000001,
                                    'G001076295': 0.10000000000000001,
                                    'G001380675': 0.0},
                      'score_fna': {'G000014725': 1.0,
                                    'G000254175': 0.48599999999999999,
                                    'G000775715': 0.90200000000000002,
                                    'G000881615': 1.0,
                                    'G000955195': 1.0,
                                    'G001076295': 0.78500000000000003,
                                    'G001380675': 0.0},
                      'score_rrna': {'G000014725': 1.0,
                                     'G000254175': 1.0,
                                     'G000775715': 1.0,
                                     'G000881615': 0.0,
                                     'G000955195': 0.0,
                                     'G001076295': 0.90000000000000002,
                                     'G001380675': 0.0},
                      'score_trna': {'G000014725': 1.0,
                                     'G000254175': 0.80000000000000004,
                                     'G000775715': 0.90000000000000002,
                                     'G000881615': 0.0,
                                     'G000955195': 0.0,
                                     'G001076295': 0.80000000000000004,
                                     'G001380675': 0.0}}

        exp_avg_scores = {'G000881615': 0.27500,
                          'G000014725': 0.77500,
                          'G000775715': 0.72550,
                          'G000254175': 0.59650,
                          'G001380675': 0.00000,
                          'G000955195': 0.27500,
                          'G001076295': 0.64625}

        input_df = pd.DataFrame(input_dict)
        obs_scores = calc_avg_score(input_df)
        for item in input_df.index:
            self.assertAlmostEqual(obs_scores[item],
                                   exp_avg_scores[item])

    def test__main(self):
        """ test for the main process following Click.
        """
        params = [['--repophlan_fp', self.repophlan_fp],
                  ['--score_threshold', self.score_threshold],
                  ['--output_fp', self.output_fp],
                  ['--score_cols', self.score_cols],
                  ['--fp_cols', self.fp_cols],
                  ['--avg', self.avg]]
        res = CliRunner().invoke(_main, params)
        self.assertEqual(res.exit_code, -1)


if __name__ == '__main__':
    main()
