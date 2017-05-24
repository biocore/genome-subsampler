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

from genomesubsampler.genomeSubsampler import _main


class GenomeSubsamplerTests(TestCase):
    """Tests for genomeSubsampler.py."""

    def setUp(self):
        """Create working directory and test files.
        """
        self.wkdir = mkdtemp()

    def tearDown(self):
        """Delete working directory and test files.
        """
        rmtree(self.wkdir)

    def test__main(self):
        """Test for the main process following Click.
        """
        params = []
        res = CliRunner().invoke(_main, params)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.output, 'Task completed.\n')


if __name__ == '__main__':
    main()
