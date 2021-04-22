#!/usr/bin/env python3

"""
    Usage: 
    cat infile | - ./extract_column.py A,B,C > outfile
    ./extract_column.py infile A,B,C > outfile
"""
import sys,argparse
from bio.base import LineFileSpliter

def extract_by_index(infile,outfile,index,header=[],sep="\t"):
    if header:
        # header is True, output _header of LineFileSpliter
        if not isinstance(header,list):
            header = infile._header
        outfile.write(sep.join(map(lambda x:header[x],index))+"\n")
    for la in infile:
        outfile.write(sep.join(map(lambda x:la[x],index))+"\n")

def extract_by_match(infile,outfile,match_str,write_header=True,sep="\t"):
    index = []
    for s in match_str:
        index.append(infile._header.index(s))
    extract_by_index(infile,outfile,index,write_header,sep)

def main(argv):

    parser = argparse.ArgumentParser(description="Extract columns")
    parser.add_argument('infile',nargs='?',help="file to be filtered, \"-\" for stdin ")
    parser.add_argument('key',nargs='?',help="Which columns to be select. eg 0,3,4 or Name,aa,comment ")
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-s','--sep',nargs='?',default="\t")
    parser.add_argument('-n','--header',action='store_false',help="Don't output header")
    parser.add_argument('-m','--mode',nargs='?',default='match',choices=['index','match'])
    args = parser.parse_args(argv[1:])

    # LineFileSpliter can handle "-"
    myfile = LineFileSpliter(args.infile,has_header=args.header,sep=args.sep)
    if args.mode == "index":
        index = list(map(int,args.key.split(",")))
        header = args.header
        extract_by_index(myfile,args.outfile,index,header,args.sep)
    elif args.mode == "match":
        index = args.key.split(",")
        extract_by_match(myfile,sys.stdout,index,args.header,args.sep)

if __name__ == '__main__':
    import sys
    main(sys.argv)
