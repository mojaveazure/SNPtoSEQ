#!/usr/bin/env python3

"""Get the contextual sequence around SNPs from a VCF file"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")


import os, argparse

try:
    # import vcf
    from Bio import SeqIO
    from overload import overload
except ImportError as error:
    sys.exit("Please install " + error.name)


_ARGUMENTS = (
    'vcf',
    'reference',
    'outname',
    'window',
    'no_fasta',
    'no_bed'
    # 'gff'
)
_WINDOW = 120

#   A class for a line of a BED file
class Bed(object):
    """This class describes a line of a BED file"""

    def __init__(self, chrom, start, end, name=None, score=0, strand=None):
        try:
            assert isinstance(chrom, (str, int))
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert isinstance(name, str) or name is None
            assert isinstance(score, int)
            assert strand in ('+', '-') or strand is None
        except AssertionError:
            raise TypeError("'chrom' must be of type 'str' or 'int'; 'start' and 'end' must be of type 'int'; 'name' must be of type 'str' or None; 'score' must be of type 'int'; 'strand' must be '+', '-', or None")
        self._chrom = str(chrom)
        self._start = start
        self._end = end
        self._name = name
        self._score = score
        self._strand = strand

    def __repr__(self):
        return self._name + '@' + self._chrom + ':' + str(self._start) + '-' + str(self._end)

    def get_chrom(self):
        """Get the chromosome"""
        return self._chrom

    def get_start(self):
        """Get the start position"""
        return self._start

    def get_end(self):
        """Get the end position"""
        return self._end

    def capture_sequence(self, reference):
        """Capture the sequence defined by this region"""
        try:
            assert isinstance(reference, dict)
        except AssertionError:
            raise TypeError("'reference' must be a dictionary")
        try:
            assert self._chrom in reference.keys()
        except AssertionError:
            raise ValueError("Cannot find %s in the reference genome" % self._chrom)
        sequence = str(reference[self._chrom][self._start:self._end].seq)
        if self._name:
            name = self._name
        else:
            name = self._chrom + ':' + str(self._start) + '-' + str(self._end)
        seq = Seq(seqid=name, sequence=sequence)
        return seq

    def format_bed(self):
        """Return this annotation in BED format"""
        bedline = [
            self._chrom,
            str(self._start),
            str(self._end),
        ]
        for option in (self._name, self._score, self._strand):
            if option or option is 0:
                bedline.append(str(option))
            else:
                break
        return '\t'.join(bedline)


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

    def __repr__(self):
        return self._seqid

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


#   Make an argument parser
def _arguments():
    (vcffile, reference, outname, window, no_fasta, no_bed) = _ARGUMENTS
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Capture SNP contextual sequences into BED and/or FASTA format"
    )
    inputs = parser.add_argument_group(
        title='Input Options',
        description='Provide a VCF file of SNPs and a reference genome in FASTA format'
    )
    inputs.add_argument( # VCF file
        '-v',
        '--vcf',
        dest=str(vcffile),
        type=str,
        default=None,
        required=True,
        metavar='VCF FILE',
        help="SNPs in VCF format"
    )
    inputs.add_argument( # Reference genome
        '-r',
        '--reference',
        dest=str(reference),
        type=str,
        default=None,
        required=True,
        metavar='REFERENCE FASTA FILE',
        help="Reference genome in FASTA format"
    )
    windowsize = parser.add_argument_group(
        title='Window options',
        description='Select a one-sided window size to capture in contextual sequence; this size represents the amount of sequence on each side of the SNP to capture'
    )
    windowsize.add_argument( # Window size
        '-w',
        '--window',
        dest=str(window),
        type=int,
        default=_WINDOW,
        required=False,
        metavar='WINDOW SIZE',
        help="Set the window size, defaults to %d" % _WINDOW
    )
    outputs = parser.add_argument_group(
        title='Output Options',
        description='Provide an output name as well as choose if we suppress some output'
    )
    outputs.add_argument( # Outname
        '-o',
        '--outname',
        dest=str(outname),
        type=str,
        default=None,
        required=False,
        metavar='OUTPUT NAME',
        help="Basename for the output file(s), defaults to the basename of '-v | --vcf'"
    )
    switches = outputs.add_mutually_exclusive_group(required=False)
    switches.add_argument( # No FASTA ouptut switch
        '--no-fasta',
        dest=str(no_fasta),
        action='store_const',
        const=True,
        default=False,
        required=False,
        metavar='NO FASTA OUTPUT',
        help="Do we suppress FASTA output? Incompatible with '--no-bed'"
    )
    switches.add_argument( # No BED output switch
        '--no-bed',
        dest=str(no_bed),
        action='store_const',
        const=True,
        default=False,
        required=False,
        metavar='NO BED OUTPUT',
        help="Do we suppress BED output? Incompatible with '--no-fasta'"
    )
    return parser


#   Validate arguments
def _validate_args(args):
    try:
        assert isinstance(args, dict)
    except AssertionError:
        raise TypeError("'args' must be a dictionary!")
    try:
        assert set(_ARGUMENTS) == set(args.keys())
    except AssertionError:
        raise ValueError("Your arguments do not match the required argument list!")
    #   Temporary shit
    gff = 'gff'
    args[gff] = False
    (vcffile, reference, outname, window, no_fasta, no_bed) = _ARGUMENTS
    if not os.path.isfile(args[vcffile]):
        raise FileNotFoundError("%s is not a file!" % args[vcffile])
    if not os.path.isfile(args[reference]):
        raise FileNotFoundError("%s doesn't exist!" % args[reference])
    try:
        assert args[outname]
        assert args[window]
    except AssertionError:
        raise argparse.ArgumentTypeError("Something fucked up")
    if args[no_bed] and (args[gff] or args[no_fasta]):
        raise argparse.ArgumentTypeError("'--no-bed' is incompatible with '--no-fasta' and '--gff'")
    if args[gff] and args[no_bed]:
        raise argparse.ArgumentTypeError("'--gff' is incompatible with '--no-bed'")
    if args[no_fasta] and args[no_bed]:
        raise argparse.ArgumentTypeError("'--no-fasta' is incompatible with '--no-bed'")


#   Funky subtraction
def funky_sub(center, windowsize, maximum=None):
    """Get lower and upper bounds around a center, bounded by 0 and an optional upper bound"""
    try:
        assert isinstance(center, int)
        assert isinstance(windowsize, int)
        assert isinstance(maximum, int) or maximum is None
    except AssertionError:
        raise TypeError("'center' and 'windowsize' must be of type 'int'; 'maximum' must be of type 'int' or None")
    lower = center - windowsize
    if lower < 0:
        lower = 0
    upper = center + windowsize
    if maximum and upper > maximum:
        upper = maximum
    return(lower, upper)


#   Capture region in BED format
@overload
def capture_region(chromosome, position, name, reference, window):
    """Capture a region defined by a window center and window size"""
    try:
        assert isinstance(chromosome, (str, int))
        assert isinstance(position, int)
        assert isinstance(name, str)
        assert isinstance(window, int)
        assert isinstance(reference, dict)
    except AssertionError:
        raise TypeError("'chromosome' must be of type 'str' or 'int'; 'position' must be of type 'int'; 'name' must be of type 'str'; 'window' must be of type 'int', 'reference' must be a dictionary")
    try:
        assert chromosome in reference.keys()
    except AssertionError:
        raise ValueError("Cannot find chromosome %s in the reference" % chromosome)
    lower, upper = funky_sub(center=position, windowsize=window)
    return Bed(chrom=chromosome, start=lower-1, end=upper, name=name, strand='+')


@capture_region.add
def capture_region(snp_tuple, reference, window):
    """Capture a region defined by a window center and window size"""
    try:
        assert isinstance(snp_tuple, tuple)
        (chromosome, position, name) = snp_tuple
        assert isinstance(chromosome, (str, int))
        assert isinstance(position, int)
        assert isinstance(name, str)
        assert isinstance(reference, dict)
        assert isinstance(window, int)
    except AssertionError:
        raise TypeError("'snp_tuple' must be a tuple containing three objects of types 'str/int', 'int', and 'str', in that order; 'reference' must be a dictionary; 'window' must be of type 'int'")
    return capture_region(
        chromosome=chromosome,
        position=position,
        name=name,
        reference=reference,
        window=window
    )


#   Main
def main():
    """potato"""
    #   Handle arguments
    (vcffile, reference, outname, window, no_fasta, no_bed) = _ARGUMENTS # Get the names of our arguments
    parser = _arguments() # Create an argument parser
    if not sys.argv[1:]: # If we don't have our arguments, exit printing the help message
        sys.exit(parser.print_help())
    args = vars(parser.parse_args()) # Parse our arguments into a dictionary
    if not args[outname]: # Ensure we have an output name
        args[outname] = os.path.splitext(args[vcffile])[0]
    # _validate_args(args=args) # Validate our arguments
    #   Read in the reference gneome if we have one
    print("Reading in reference FASTA file:", args[reference], file=sys.stderr)
    reference = SeqIO.to_dict(SeqIO.parse(args[reference], 'fasta'))
    #   Read in the VCF file into BED format
    snps = list() # Holding list
    print("Reading in VCF file:", args[vcffile], file=sys.stderr)
    with open(args[vcffile], 'r') as vcf:
        for line in vcf:
            if line.startswith('#'): # Skip header lines
                continue
            split = line.strip().split() # Split our line into corresponding pieces
            chrom, pos, name = split[0], int(split[1]), split[2] # Take only the parts we need
            snps.append((chrom, pos, name)) # Wrap the pieces into a tuple and add to our list of SNPs
    #   Capture regions sequences
    print("Capturing regions and sequences", file=sys.stderr)
    bedfile = tuple(capture_region(snp_tuple=t, reference=reference, window=args[window]) for t in snps) # Create a tuple of BED objects
    sequences = tuple(bed.capture_sequence(reference=reference) for bed in bedfile) # Create a tuple of sequences
    #   Write BED/GFF output
    if not args[no_bed]: # If we are writing to a BED file
        bedout = args[outname] + '.bed' # Create an output name
        with open(bedout, 'w') as bo:
            print("Writing regions to", bedout, file=sys.stderr)
            #   Sort the BED objects by chromosome, start position, and end position
            for bed in sorted(bedfile, key=lambda b: (b.get_chrom(), b.get_start(), b.get_end())):
                bo.write(bed.format_bed())
                bo.write('\n')
    #   Write FASTA output
    if not args[no_fasta]: # If we are writing to a FASTA file
        fastaout = args[outname] + '.fasta' # Create an output name
        with open(fastaout, 'w') as fo:
            print("Writing FASTA to", fastaout, file=sys.stderr)
            for seq in sorted(sequences, key=lambda seq: seq.get_seqid()): # Sort sequences by sequence identifiers
                fo.write(seq.format_fasta())
                fo.write('\n')


if __name__ == '__main__':
    main()

