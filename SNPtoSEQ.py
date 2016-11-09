#!/usr/bin/env python3

"""Get the contextual sequence around SNPs from a VCF file"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")


import os, argparse

try:
    import vcf
    from Bio import SeqIO
except ImportError as error:
    sys.exit("Please install " + error.name)


class Bed(object):
    """This class describes a line of a BED file"""

    def __init__(self, chrom, chrom_start, chrom_end, name=None, score=0, strand='+'):
        try:
            assert isinstance(chrom, (str, int))
            assert isinstance(chrom_start, int)
            assert isinstance(chrom_end, int)
            assert isinstance(name, str) or name is None
            assert isinstance(score, int) or score is None
            assert strand in ('+', '-') or strand is None
        except AssertionError:
            raise TypeError
        self._chrom = str(chrom)
        self._start = chrom_start
        self._end = chrom_end
        self._name = name
        self._score = score
        self._strand = strand

    def format_bed(self):
        """Return this annotation in BED format"""
        bedline = [
            self._chrom,
            str(self._start),
            str(self._end),
        ]
        for option in (self._name, self._score, self._strand):
            if option:
                bedline.append(option)
            else:
                break
        return '\t'.join(bedline)

    def format_gff(self):
        """Return this annotation in GFF format"""
        return NotImplemented


class Seq(object):
    """This class describes a sequence"""

    def __init__(self, seqid, sequence):
        try:
            assert isinstance(seqid, (str, int, float))
            assert isinstance(sequence, str)
        except AssertionError:
            raise TypeError("'seqid' must be of type 'str', 'int' or 'float'; 'sequence' must be of type 'str'")
        self._seqid = str(seqid)
        self._sequence = sequence

    def get_seqid(self):
        """Get the SeqID"""
        return self._seqid

    def get_sequence(self):
        """Get the sequence"""
        return self._sequence

    def format_fasta(self):
        """Return this sequence in FASTA format"""
        seq = (
            '>' + self._seqid,
            self._sequence
        )
        return '\n'.join(seq)
