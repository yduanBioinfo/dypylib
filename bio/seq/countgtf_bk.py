#!/usr/bin/env python3

import sys
from bio.seq.base import Gff

'''
Count transcript and gene of input gffs
input: gffs
output:
gffname  gene_count tx_count mono_exon_tx two_exon_tx more_than_two_exon_tx tx_per_gene ratio_of_mono_exon_tx ratio_of_two_exon_tx ratio_of_more_than_two_exon_tx
'''

def count_tx_and_gene(gff):
    """ Return: count of gene; count of transripts; tx per gene.
    """
    tx = set()
    gene = set()
    for rec in gff:
        if rec.type == "exon":
            tx.add(rec.tx_id)
            gene.add(rec.gene_id)
    return len(gene), len(tx), len(tx)/len(gene)

def count(gffs,outfile,sep = "\t",outnames=[]):
    for i in range(len(gffs)):
        gff=gffs[i]
        if outnames:
            name = outnames[i]
        else:
            name = ".".join(gff.split("/")[-1].split(".")[:-1])
        outfile.write(name+sep+sep.join(map(str,count_tx_and_gene(Gff(gff))))+"\n")

def main(argv):
    import argparse

    parser = argparse.ArgumentParser(description="Count tx and gens in gff")
    parser.add_argument('infiles',nargs='+',help="input gffs")
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-n','--names',nargs='?',help="gff names. Seperated by comma.")
    parser.add_argument('-s','--sep',nargs='?',default="\t")
    args = parser.parse_args(argv[1:])

    if args.names:
        names = args.names.split(",")
    else:
        names = []
    count(args.infiles,args.outfile,args.sep,names)

if __name__ == '__main__':

    import sys
    main(sys.argv)
