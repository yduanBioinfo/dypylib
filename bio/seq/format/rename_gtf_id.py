#!/usr/bin/env python3

import sys, argparse
from bio.seq.base import GtfDict

def rename_gtf(gtf, outfile, gene_id, tx_id, id_mapping=None):
    if id_mapping:
        genes, idmap = gtf.rename_gene_tx(gene_id, tx_id, get_id = True)
    else:
        genes = gtf.rename_gene_tx(gene_id, tx_id)
    for gid, gene in genes.items():
        for tid, tx in gene.items():
            for exon in tx.values():
                outfile.write(exon.as_str()+"\n")
    for idp in idmap:
        id_mapping.write("\t".join(idp)+"\n")

def main(argv):
    parser = argparse.ArgumentParser(description="Rename Gtf Gene id and tx id. (The groups are based on the old ID.)")
    parser.add_argument('gtf',nargs='?',help="input gtf file.",default = sys.stdin,type=GtfDict)
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-g','--gene_id',nargs='?',default="HBG",help = "Prefix of gene id")
    parser.add_argument('-t','--tx_id',nargs='?',default=None,help = "Prefix of transcript id")
    parser.add_argument('-r','--refer-table',nargs='?',default=None,help ="Output tabular file to record old/new id",type=argparse.FileType('w'))
    args = parser.parse_args(argv[1:])

    rename_gtf(args.gtf,args.outfile,args.gene_id,args.tx_id,args.refer_table)

if __name__ == '__main__':
    main(sys.argv)
