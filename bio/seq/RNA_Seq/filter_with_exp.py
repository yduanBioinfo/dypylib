#!/usr/bin/env python3

"""Filter transcripts with its expression value.
Expression table and threshold is needed for filtering.

Two mode is provide,
1. Output the name of transcripts that should be filter out (or kept).
    cat expression.tpm | ./filter_with_exp.py - -t 1 -c 1 > filter_out.names
    ./filter_with_exp.py expression.tpm -t 2 -c 2 --write-keep > kept.names
2. Filtering GTF file.
    ./filter_with_exp.py expression.tpm -g input.gtf -t 2 -c 2 > new.gtf
"""

import sys, argparse
import pandas as pd
from dypylib.bio.base import DictFile
from dypylib.bio.seq.base import Gff

def bigger_than_thre(s, thre=1, count=2):
    """ If the number of value that is big than 'thre' is more than
    'count'. """
    return s[s > thre].count() > count

def get_count_df(exp, thre, count):
    """Main function.
    Get bool data: count_df from exp file"""
    # Read expression table
    exp = pd.read_table(exp,index_col=0)
    #bigger_than_thre(exp.iloc[0,:])
    count_df = exp.apply(lambda x:bigger_than_thre(x, thre,\
        count),axis=1)
    return count_df

def filter_gff(mygff, keep, outfile):
    for rec in mygff:
        if rec.tx_id in keep:
            outfile.write(rec.as_str()+"\n")

def main(argv):
    parser = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        description=__doc__)
    parser.add_argument('exp',nargs='?',help="Expression file.")
    parser.add_argument('-o','--outfile',nargs='?',\
        help="output file",default=sys.stdout,\
        type=argparse.FileType('w'))
    parser.add_argument('--write-keep',help=("Write names of "
        "transcipt should be kept, rather than filter out "
        "(default). Conflict with -g").format(),\
        action = 'store_true')
    parser.add_argument('-t','--threshold',nargs='?', \
        help="expression threshold",type=float,default=1)
    parser.add_argument('-c','--count',nargs='?',type=int,default=2,\
        help=('least number of sample a transcript '
        'should expressed in').format())
    parser.add_argument('-g','--GTF',nargs='?',\
        help=('Provide GTF file for filtering.'))
    args = parser.parse_args(argv[1:])

    count_df = get_count_df(args.exp, args.threshold, args.count)
    
    if args.GTF:
        # Output the good transcripts.
        mygtf = Gff(args.GTF)
        filter_gff(mygtf,count_df[count_df],args.outfile)
        return
    
    if args.write_keep:
        args.outfile.write("\n".join(count_df[count_df].index)+"\n")
    else:
        args.outfile.write("\n".join(count_df[~ count_df].index)+"\n")

if __name__ == '__main__':
    main(sys.argv)
