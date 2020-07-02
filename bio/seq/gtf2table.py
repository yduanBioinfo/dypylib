#!/usr/bin/env python3

import sys
from bio.seq.base import Gff, GtfDict

'''
input gtf should only contain exon feature.
stat for each transcript
TranscriptID,GeneID,length,count_of_exons,length_of_exons,count_of_introns,leng_of_introns
'''

def iter_stats(gtfdict):
    """ Get tx stats from GtfDict
    """
    txdict = gtfdict.get_TxDict()
    for key, val in txdict.items():
        tx_id = key
        gene_id = val.gene_id
        length = val.get_length()
        count_of_exons = val.get_count()
        length_of_exons = list(map(lambda x:len(x),val.values()))
        count_of_introns = count_of_exons - 1
        length_of_introns = list(map(lambda x:len(x),val.get_introns()))
        yield tx_id, gene_id, length, count_of_exons, length_of_exons, count_of_introns, length_of_introns

def iter_stats_str(gtfdict):
    """ Concate list element with ","
    """
    def _join(lst,sep=","):
        return sep.join(map(str,lst))
    for ss in iter_stats(gtfdict):
        res = list(ss)
        res[4] = _join(res[4])
        res[6] = _join(res[6])
        res = map(str,res)
        yield res

def gtf2table(gtfdict,outfile,write_header=None,sep="\t"):
    header = ["TranscriptID","GeneID","length","count_of_exons","length_of_exons","count_of_introns","length_of_introns"]
    if write_header:
        outfile.write(sep.join(header)+"\n")
    for line in iter_stats_str(gtfdict):
        outfile.write(sep.join(line)+"\n")

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Get a table description for gtf file.")
    parser.add_argument('infile',nargs='?',help="input gtf file",type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help="outfile.",type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('-t','--gff-type',nargs='?',default='normal',choices=['normal','gc'])
    parser.add_argument('-n','--none-header',action='store_false')
    args = parser.parse_args(argv[1:])

    mygtfdict = GtfDict(args.infile,fm=args.gff_type)
    gtf2table(mygtfdict,args.outfile,write_header=args.none_header)

if __name__ == '__main__':
    import sys
    main(sys.argv) 
