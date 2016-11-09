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


_ARGUMENTS = (
    'vcf',
    'reference',
    'outname',
    'no_fasta',
    'no_bed',
    'gff'
)

def _outname(value):
    try:
        assert isinstance(value, str) or value is sys.stdout
        return value
    except AssertionError:
        raise argparse.ArgumentTypeError("%s is not a valid output name" % value)


#   A class for a line of a BED file
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


#   A class for a sequence
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


def _arguments():
    #   VCF in, always required
    #   Reference in, required unless --no-fasta
    #   Output name, defaults to stdout
    #   --no-fasta, no FASTA output, mutually exclusive with --no-bed
    #   --no-bed, no BED output, mutually exclusive with --no-bed and --gff
    #   --gff, output GFF instead of BED, mutually exclusive with --no-bed
    (vcffile, reference, outname, no_fasta, no_bed, gff) = _ARGUMENTS
    parser = argparse.ArgumentParser()
    inputs = parser.add_argument_group(
        title='Input Options',
        description='Hot potato'
    )
    inputs.add_argument( # VCF File
        '-v',
        '--vcf',
        dest=str(vcffile),
        type=str,
        default=None,
        required=True,
        metavar='VCF FILE',
        help="Input VCF file"
    )
    inputs.add_argument(
        '-r',
        '--reference',
        dest=str(reference),
        type=str,
        default=None,
        required=False,
        metavar='REFERENCE FASTA FILE',
        help="Reference genome in FASTA format"
    )
    outputs = parser.add_argument_group(
        title='Output Options',
        description='Set output options'
    )
    outputs.add_argument(
        '-o',
        '--outname',
        dest=str(outname),
        type=_outname,
        default=sys.stdout,
        required=False,
        metavar='OUTPUT NAME',
        help="Basename for the output file(s), defaults to stdout"
    )
    outputs.add_argument(
        '--no-fasta',
        dest=str(no_fasta),
        action='store_const',
        const=True,
        default=False,
        required=False,
        metavar='NO FASTA OUTPUT',
        help="Do we suppress FASTA output? Incompatible with '--no-bed'"
    )
    outputs.add_argument(
        '--no-bed',
        dest=str(no_bed),
        action='store_const',
        const=True,
        default=False,
        required=False,
        metavar='NO BED OUTPUT',
        help="Do we suppress BED output? Incompatible with '--no-fasta' and '--gff'"
    )
    outputs.add_argument(
        '--gff',
        dest=str(gff),
        action='store_const',
        const=True,
        default=False,
        required=False,
        metavar='GFF OUTPUT',
        help="Output GFF instead of BED format, incompatible with '--no-bed'"
    )
    return parser


def _validate_args(args):
    try:
        assert isinstance(args, dict)
    except AssertionError:
        raise TypeError("Args must be a dictionary!")
    try:
        assert set(_ARGUMENTS) == set(args.keys())
    except AssertionError:
        raise ValueError("Your arguments do not match the required argument list!")
    (vcffile, reference, outname, no_fasta, no_bed, gff) = _ARGUMENTS
    try:
        assert os.path.isfile(args[vcffile])
    except AssertionError:
        raise FileNotFoundError("%s is not a file!" % args[vcffile])
    if not args[reference] and not args[no_fasta]:
        raise argparse.ArgumentTypeError("%s requires a reference FASTA file to extract sequences from" % sys.argv[0])
    try:
        assert args[outname]
    except AssertionError:
        raise argparse.ArgumentTypeError("Something fucked up")
    if args[no_bed] and (args[gff] or args[no_fasta]):
        raise argparse.ArgumentTypeError("'--no-bed' is incompatible with '--no-fasta' and '--gff'")
    if args[gff] and args[no_bed]:
        raise argparse.ArgumentTypeError("'--gff' is incompatible with '--no-bed'")
    if args[no_fasta] and args[no_bed]:
        raise argparse.ArgumentTypeError("'--no-fasta' is incompatible with '--no-bed'")


#   Main
def main():
    """potato"""
    parser = _arguments()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args())
    _validate_args(args=args)
    sys.exit(args)
    try:
        reference = SeqIO.to_dict(SeqIO.parse(args['reference'], 'fasta'))
    except FileNotFoundError as error:
        sys.exit("Failed to find " + error.filename)


if __name__ == '__main__':
    main()

