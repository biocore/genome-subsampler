#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

#
# main front-end of the genome-subsampler
#

import click


@click.command()
def _main():
    """Main front-end of the genome-subsampler.
    """
    click.echo('Task completed.')


if __name__ == "__main__":
    _main()
