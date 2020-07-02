#!/usr/bin/env python3

import sys
from bio.seq.base import Gff, GtfDict

'''
Count transcript and gene of input gffs
input: gffs
output:
gffname  gene_count tx_count mono_exon_tx two_exon_tx more_than_two_exon_tx tx_per_gene ratio_of_mono_exon_tx ratio_of_two_exon_tx ratio_of_more_than_two_exon_tx
'''

def count_tx_and_gene(gff):
    """ Return: count of gene; count of transripts; tx per gene.
    """
    #tx = set()
    #gene = set()
    #for rec in gff:
    #    if rec.type == "exon":
    #        tx.add(rec.tx_id)
    #        gene.add(rec.gene_id)
    #for chrdict in gff.values():
    #    print(k)
    #return len(gene), len(tx), len(tx)/len(gene)
    tx_1 = 0 # tx is constructed of 1 exon
    tx_2 = 0 # tx is constructed of 2 exons
    tx_n = 0 # tx is constructed of more than 3 exons
    tx_dict = gff.get_TxDict()
    for tx in tx_dict.values():
        if tx.get_count() == 1:
            tx_1 += 1
        elif tx.get_count() == 2:
            tx_2 += 1
        else:
            tx_n += 1
    tx_count = len(tx_dict)
    gene_count = len(gff.get_GeneDict())
    # empty gff??
    if tx_count == 0 or gene_count == 0:
        return 0,0,0,0,0,0,0,0,0
    return gene_count, tx_count, tx_count/gene_count, tx_1, tx_1/tx_count, tx_2, tx_2/tx_count, tx_n, tx_n/tx_count

def count(gffs,outfile,fm,sep = "\t",outnames=[],header=False):
    if header:
        outfile.write(sep.join(["gene_count", "tx_count", "tx_per_gene", "tx_1", "tx_1_fraction", "tx_2", "tx_2_fraction", "tx_n", "tx_n_fraction"])+"\n")
    for i in range(len(gffs)):
        gff=gffs[i]
        if outnames:
            name = outnames[i]
        else:
            name = ".".join(gff.split("/")[-1].split(".")[:-1])
        outfile.write(name+sep+sep.join(map(str,count_tx_and_gene(GtfDict(gff,fm=fm))))+"\n")

def main(argv):
    import argparse

    parser = argparse.ArgumentParser(description="Count tx and gens in gff")
    parser.add_argument('infiles',nargs='+',help="input gffs")
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-n','--names',nargs='?',help="gff names. Seperated by comma.")
    parser.add_argument('-t','--gff-type',nargs='?',default='normal',choices=['normal','gc'])
    parser.add_argument('--header',action='store_true',help="Set this option to enable header output.")
    parser.add_argument('-s','--sep',nargs='?',default="\t")
    args = parser.parse_args(argv[1:])

    if args.names:
        names = args.names.split(",")
    else:
        names = []
    count(args.infiles,args.outfile,args.gff_type,args.sep,names,args.header)

if __name__ == '__main__':

    import sys
    main(sys.argv)
