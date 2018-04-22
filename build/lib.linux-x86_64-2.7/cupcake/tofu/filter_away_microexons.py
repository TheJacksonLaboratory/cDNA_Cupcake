#!/usr/bin/env python

__author__ = 'adeslat@jax.org'

"""
with big assist from eteng@pacb.com

Filter away exons that are less then default size of 12 bp

Required input: <input_prefix>.gff, .rep.fq, .abundance.txt

Example:
    filter_away_microexon.py test.collapsed --micro_exon_size=12

"""

import os, sys
from collections import defaultdict
from csv import DictReader, DictWriter

from Bio import SeqIO

from cupcake.io import GFF
from cupcake.tofu import compare_junctions

def sanity_check_collapse_input(input_prefix):
    """
    Check that
    1. the count, gff, rep files exist
    2. the number of records agree among the three
    """
    group_filename = input_prefix + '.group.txt'
    count_filename = input_prefix + '.abundance.txt'
    gff_filename = input_prefix + '.gff'
    rep_filename = input_prefix + '.rep.fq'
    if not os.path.exists(count_filename):
        print >> sys.stderr, "File {0} does not exist. Abort!".format(count_filename)
        sys.exit(-1)
    if not os.path.exists(gff_filename):
        print >> sys.stderr, "File {0} does not exist. Abort!".format(gff_filename)
        sys.exit(-1)
    if not os.path.exists(rep_filename):
        print >> sys.stderr, "File {0} does not exist. Abort!".format(rep_filename)
        sys.exit(-1)

    pbids1 = set([r.name.split('|')[0] for r in SeqIO.parse(open(rep_filename),'fastq')])
    pbids2 = set([r.seqid for r in GFF.collapseGFFReader(gff_filename)])
    pbids3 = set(read_count_file(count_filename)[0].keys())

    if len(pbids1)!=len(pbids2) or len(pbids2)!=len(pbids3) or len(pbids1)!=len(pbids3):
        print >> sys.stderr, "The number of PBID records in the files disagree! Sanity check failed."
        print >> sys.stderr, "# of PBIDs in {0}: {1}".format(rep_filename, len(pbids1))
        print >> sys.stderr, "# of PBIDs in {0}: {1}".format(gff_filename, len(pbids2))
        print >> sys.stderr, "# of PBIDs in {0}: {1}".format(count_filename, len(pbids3))
        f = open('pbids.rep.fq.txt', 'w')
        for i in pbids1:
           f.write(i)
           f.write("\n")
        f.close
        f = open('pbids.gff.txt','w')
        for i in pbids2:
           f.write(i)
           f.write("\n")
        f.close
        f = open('pbids.abundance.txt','w')
        for i in pbids3:
           f.write(i)
           f.write("\n")
        f.close
        sys.exit(-1)

    return count_filename, gff_filename, rep_filename


def read_count_file(count_filename):
    f = open(count_filename)
    count_header = ''
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith('#'):
            f.seek(cur_pos)
            break
        else:
            count_header += line
    d = dict((r['pbid'], r) for r in DictReader(f, delimiter='\t'))
    f.close()
    return d, count_header


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix", help="Input prefix (ex: filtered.microexon)")
    parser.add_argument("--micro_exon_size", type=int, default=12, help="Filter away microexons < micro_exon_size (default: 12bp)")

    args = parser.parse_args()
    output_prefix = args.input_prefix + '.filtered'

    args = parser.parse_args()
    output_prefix = args.input_prefix + '.filtered.microexon'

    count_filename, gff_filename, rep_filename = sanity_check_collapse_input(args.input_prefix)

    recs = defaultdict(lambda: [])
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        assert r.seqid.startswith('PB.')
        recs[int(r.seqid.split('.')[1])].append(r)

    good = []
    f = open(output_prefix + '.gff', 'w')
    keys = recs.keys()
    keys.sort()
    for k in recs:
        xxx = recs[k]
        for r in xxx:
          min_exon_size = min(e.end-e.start for e in r.ref_exons)
          if min_exon_size > 12: # minimum exon must be > default 12 bp
            GFF.write_collapseGFF_format(f, r)
            good.append(r.seqid)

    f.close()

    # read abundance first
    d, count_header = read_count_file(count_filename)

    # write output rep.fq
    f = open(output_prefix + '.rep.fq', 'w')
    for r in SeqIO.parse(open(rep_filename), 'fastq'):
        if r.name.split('|')[0] in good:
            SeqIO.write(r, f, 'fastq')
    f.close()

    # write output to .abundance.txt
    f = open(output_prefix + '.abundance.txt', 'w')
    f.write(count_header)
    writer = DictWriter(f, fieldnames=['pbid','count_fl','count_nfl','count_nfl_amb','norm_fl','norm_nfl','norm_nfl_amb'], \
                        delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for k in good:
        r = d[k]
        writer.writerow(r)
    f.close()

    print >> sys.stderr, "Output written to:", output_prefix + '.gff'
    print >> sys.stderr, "Output written to:", output_prefix + '.rep.fq'
    print >> sys.stderr, "Output written to:", output_prefix + '.gff'

if __name__ == "__main__":
    main()
