#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, genome-subsampler development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Calculate MinHash distances for input genomes.
"""

import os
import gzip
from subprocess import call


def compute_mash_distance(genome_dir, genome_ext='fna', cpus=1):
    """Compute pairwise MinHash distances between genomes in a given path.

    Parameters
    ----------
    genome_dir : str
        directory containing input genomes
    genome_ext : str
        extension name of input genome files (default: fna)
    cpus : int
        number of CPU cores to use (default: 1)
    """
    with open('genomes.txt', 'w') as fh:
        for fname in os.listdir(genome_dir):
            g = fname.split('.')[0]
            fpath = os.path.join(genome_dir, fname)
            if fname.endswith('.%s.gz' % genome_ext):
                call('zcat %s > %s' % (fpath, g), shell=True)
            elif fname.endswith('.%s.bz2' % genome_ext):
                call('bzcat %s > %s' % (fpath, g), shell=True)
            elif fname.endswith('.%s.xz' % genome_ext):
                call('xzcat %s > %s' % (fpath, g), shell=True)
            elif fname.endswith('.fna'):
                call('ln -s %s %s' % (fpath, g), shell=True)
            call('mash sketch %s' % g, shell=True)
            call('rm %s' % g, shell=True)
            fh.write('%s\n' % g)

    # paste sketches in one file
    call('sed "s/$/.msh/" genomes.txt  > genome_mshs.txt', shell=True)
    call('mash paste -l genomes.msh genome_mshs.txt', shell=True)
    call('rm genome_mshs.txt', shell=True)

    # compute pairwise distances based on sketches
    call('mash dist -p %s genomes.msh genomes.msh | gzip > genomes.dist.gz'
         % cpus, shell=True)


def make_distance_matrix(dist_dir):
    """Compute distance matrix based on pairwise MinHash distances

    Parameters
    ----------
    dist_dir: str
        path to directory with pairwise MinHash distances between genomes
    """
    ids = []
    vals = []
    reading_ids = True
    fpath = os.path.join(dist_dir, 'genomes.dist.gz')
    with open('dist_matrix.txt', 'w') as fh:
        with gzip.open(fpath, 'r') as fdist:
            for line in fdist:
                line=line.decode("utf-8")
                source, target, distance = line.rstrip('\r\n').split('\t')[:3]
                if reading_ids:
                    # read ids and values of the first id
                    if not ids and source != target:
                        raise ValueError('First source'
                                         'and target do not match')
                    if not ids or source != ids[0]:
                        ids.append(source)
                        vals.append(distance)
                    else:
                        # reading ids is completed (the same id occurs again)
                        reading_ids = False
                        # write column headers (ids) and values of the first id
                        fh.write('\t%s\n' % '\t'.join(ids))
                        fh.write('%s\t%s\n' % (ids[0], '\t'.join(vals)))
                        # read first value of second row (target)
                        isrc, itag, ctag, vals = 1, 1, target, [distance]
                        if target != ids[1]:
                            raise ValueError('Second source'
                                             'and target do not match')
                else:
                    # read and write the rest of values
                    if target == ctag:
                        # add one value to current row (target)
                        if source != ids[isrc]:
                            raise ValueError('Pair %s-%s is not in order.'
                                             % (source, target))
                        vals.append(distance)
                        isrc += 1
                    else:
                        # finish last row by writing it, and move to target row
                        itag += 1
                        if target != ids[itag]:
                            raise ValueError('Target IDs are not in order.')
                        if source != ids[0]:
                            raise ValueError('Source IDs are not in order.')
                        fh.write('%s\t%s\n' % (ctag, '\t'.join(vals)))
                        isrc, ctag, vals = 1, target, [distance]
        fh.write('%s\t%s\n' % (ctag, '\t'.join(vals)))
