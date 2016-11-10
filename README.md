# SNPtoSEQ.py

`SNPtoSEQ.py` is a simple Python program that captures SNP contextual sequence and provides a BED and/or FASTA file with the sequences

```
$ ./SNPtoSEQ.py
usage: SNPtoSEQ.py -v VCF FILE -r REFERENCE FASTA FILE [-w WINDOW SIZE]
                   [-o OUTPUT NAME] [--no-fasta | --no-bed]

Capture SNP contextual sequences into BED and/or FASTA format

Input Options:
  Provide a VCF file of SNPs and a reference genome in FASTA format

  -v VCF FILE, --vcf VCF FILE
                        SNPs in VCF format
  -r REFERENCE FASTA FILE, --reference REFERENCE FASTA FILE
                        Reference genome in FASTA format

Window options:
  Select a one-sided window size to capture in contextual sequence; this
  size represents the amount of sequence on each side of the SNP to capture

  -w WINDOW SIZE, --window WINDOW SIZE
                        Set the window size, defaults to 120

Output Options:
  Provide an output name as well as choose if we suppress some output

  -o OUTPUT NAME, --outname OUTPUT NAME
                        Basename for the output file(s), defaults to the name
                        of '-v | --vcf' with modified extensions
  --no-fasta            Do we suppress FASTA output? Incompatible with '--no-
                        bed'
  --no-bed              Do we suppress BED output? Incompatible with '--no-
                        fasta'
```

## Inputs

`SNPtoSEQ.py` requires a VCF file with SNP information and a reference genome in FASTA format. Please ensure that the chromosome information in the VCF matches the sequence identifiers in the FASTA file. To check, use UNIX `grep` and `cut` to find the chromosome information in both files:

```bash
grep -v '#' ${MY_VCF} | cut -f 1
grep '>' ${MY_REFERENCE}
```

## Outputs

`SNPtoSEQ.py` creates between one and two output files:

| File name | Contents |
| --------- | -------- |
| *output*.bed | BED file describing the SNP contextual sequence; suppressed with `--no-bed` |
| *output*.fasta | FASTA file with containing the SNP contextual sequences; sequence identifiers are the SNP name; suppressed-with `--no-fasta` |

## Dependencies

`SNPtoSEQ.py` depends on the following:
 - [Python 3](https://www.python.org/)
 - [BioPython](http://biopython.org/wiki/Biopython)
 - [overload](https://pypi.python.org/pypi/overload)

BioPython and overload are available through [PyPi](https://pypi.python.org/pypi) and can be downloaded using [pip3](https://pip.pypa.io/en/latest/installing/) (included with Python 3.4 or greater)