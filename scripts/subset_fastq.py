#!/usr/bin/env python3

"""Sample fastq file.
Usage: cat input.fastq | ./subset_fastq.py - proportion > out.fastq
example: cat input.fastq | ./subset_fastq.py - 0.01 > out.fastq 
[Subset 1% reads of input.fastq to new file: out.fastq]
"""
import sys

def iter_fastq(infile):
    data = []
    for line in infile:
        if len(data) == 4:
            yield data
            data = []
        data.append(line)

def conv_proportion_to_skip_number(num):
    ''' Convert the subseting proportion to the number of records
        to skip between two output records.
    '''
    if num > 1:
        sys.exit("The proportion must be smaller than 1. While {0} is given.\n".format(num))
    num_skip = int(1/num - 1)
    return num_skip

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(\
    formatter_class = argparse.RawDescriptionHelpFormatter,\
    description = __doc__)

    parser.add_argument('infile',nargs='?',help="Input fastq file, \"-\" for stdin ",type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-p','--proportion',nargs='?',help="Proportion of reads to be output",type=float)
    args = parser.parse_args(argv[1:])

    # ns: number of records to skip
    ns = conv_proportion_to_skip_number(args.proportion)
    fq_iterator = iter_fastq(args.infile)
    # cycle: 0,1,2,3,4,....1+ns
    # The 0-th records will be output
    # And the other records is omitted.
    count = 0
    for data in fq_iterator:
        count += 1
        if count == 1:
            args.outfile.write("".join(data))
        if count == ns+1:
            count = 0

if __name__ == '__main__':

    main(sys.argv)
