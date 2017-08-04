# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join
from tempfile import mkdtemp
from shutil import rmtree
import gzip

from skbio.util import get_data_path

from genomesubsampler.mashDistance import (
    compute_mash_distance,
    make_distance_matrix)


class MashDistanceTests(TestCase):
    """Tests for mashDistance.py."""

    def setUp(self):
        """Create working directory and test files.
        """
        self.wkdir = mkdtemp()
        self.genome_dir = get_data_path('test_fna')
        self.out_mash_dist_fp = join(self.wkdir, 'out.dist.gz')
        self.mash_dist_fp = get_data_path('genomes.dist.gz')
        self.out_dist_mat_fp = join(self.wkdir, 'out_dist_mat.txt')
        self.dist_mat_fp = get_data_path('dist_matrix.txt')

    def tearDown(self):
        """Delete working directory and test files.
        """
        rmtree(self.wkdir)

    def test_compute_mash_distance(self):
        """Test the output of function "compute_mash_distance.
        """
        compute_mash_distance(self.genome_dir, self.out_mash_dist_fp)
        with gzip.open(self.out_mash_dist_fp, 'r') as fh:
            obs = fh.read()
        with gzip.open(self.mash_dist_fp, 'r') as fh:
            exp = fh.read()
        self.assertEqual(obs, exp)

    def test_make_distance_matrix(self):
        """Test the output of function "make_distance_matrix.
        """
        make_distance_matrix(self.mash_dist_fp, self.out_dist_mat_fp)
        with open(self.out_dist_mat_fp, 'r') as fh:
            obs = fh.read()
        with open(self.dist_mat_fp, 'r') as fh:
            exp = fh.read()
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
