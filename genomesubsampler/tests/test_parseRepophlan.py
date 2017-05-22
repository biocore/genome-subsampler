# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from click.testing import CliRunner
from shutil import rmtree
from tempfile import mkdtemp
from skbio.util import get_data_path

from genomesubsampler.parseRepophlan import (parse_repophlan,
                                             _main)


class ParseRepophlanTests(TestCase):
    """ Tests for parseRepophlan.py """

    def setUp(self):
        """ Create working directory and test files
        """
        self.wkdir = mkdtemp()
        self.repophlan_fp = get_data_path('repophlan_microbes_wscores.txt')

    def tearDown(self):
        """ Delete working directory and test files
        """
        rmtree(self.wkdir)

    def test_parse_repophlan(self):
        """ Test function parse_repophlan
        """
        obs = parse_repophlan(self.repophlan_fp)
        exp = str_basic_stats.split('\n')
        self.assertListEqual(obs, exp)

    def test__main(self):
        """ Test for the main process following Click
        """
        params = ['--repophlan-wscores-fp', self.repophlan_fp]
        res = CliRunner().invoke(_main, params)
        self.assertEqual(res.exit_code, 0)
        self.assertIn(str_basic_stats, res.output)
        self.assertIn('Task completed.', res.output)


str_basic_stats = ('Total number of genomes: 9.\n'
                   'Number of RefSeq genomes: 7.\n'
                   'With genome sequences (fna): 9.\n'
                   'With protein sequences (faa): 9.\n'
                   'With protein-coding DNA sequences (ffn): 9.\n'
                   'With RNA-coding DNA sequences (frn): 9.')

if __name__ == '__main__':
    main()
