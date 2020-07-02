#!/usr/bin/env python3

import sys
from bio.seq.base import Gff

'''
input gtf should only contain exon feature.
output length, count for exon or transcripts.
'''

def make_feature_dict(mygff,key="transcript_id"):

    data={}
    for rec in mygff:
        data.setdefault(rec.attr.get(key),[]).append(rec)
    return data

def stat_db(mydb):
    #be fit for transcripts survey
    #data:(transcript name, exon counts, transcript length)
    #exon_length:length of each exon

    data=[]
    exon_length = []
    for key,value in mydb.items():
        count=len(value)
        exons_l=list(map(len,value))#length of each exon
        exon_length.extend(exons_l)
        length=sum(exons_l)
        data.append((key,count,length))
    return data,exon_length

def write_db_stat(stat_of_db,outfile,write_header=True,sep="\t"):
    header=sep.join(("Tx","count","length"))
    if write_header:
        outfile.write(header+"\n")
    for stat in stat_of_db:
        outfile.write(sep.join(map(str,stat))+"\n")

def write_exons_length(data,outfile):
    outfile.write("\n".join(map(str,data))+"\n")

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Get description for gtf file.")
    parser.add_argument('infile',nargs='?',help="input gtf file",type=argparse.FileType('r'))
    parser.add_argument('-o','--outprefix',nargs='?',help="output prefix.")
    parser.add_argument('-t','--gff-type',nargs='?',default='normal',choices=['normal','gc'])
    args = parser.parse_args(argv[1:])

    mygff=Gff(args.infile,fm=args.gff_type)
    txdb = make_feature_dict(mygff)
    stat_of_tx, length_of_exons = stat_db(txdb)
    outprefix=args.outprefix
    out_db_stat=open(outprefix+".tx_stat.txt",'w')
    out_exon_length=open(outprefix+".exon_len.txt",'w')
    write_db_stat(stat_of_tx,out_db_stat)
    write_exons_length(length_of_exons,out_exon_length)

if __name__ == '__main__':
    import sys
    main(sys.argv) 
    #mygff=Gff(sys.argv[1],fm='gc')
    #txdb = make_feature_dict(mygff)
    #stat_of_tx, length_of_exons = stat_db(txdb)
    #outprefix=sys.argv[2]
    #out_db_stat=open(outprefix+".tx_stat.txt",'w')
    #out_exon_length=open(outprefix+".exon_len.txt",'w')
    #write_db_stat(stat_of_tx,out_db_stat)
    #write_exons_length(length_of_exons,out_exon_length)
