#!/usr/bin/env python3

#from os import replace
from shutil import move
from tempfile import NamedTemporaryFile as NTmp
from dypylib.bio.base import DictFile

def replace_lst(lst,mytable):
    return list(map(lambda x:mytable.get(x,x),lst))
    
# Replace only header.
def replace_header(infile,mytable,outfile,sep="\t"):
    header = infile.readline()
    header = replace_lst(header.strip().split(sep),mytable)
    outfile.write(sep.join(header)+"\n")
    outfile.write(infile.read())
    infile.close()
    outfile.close()

# Replace whole file.
def replace_file(infile,mytable,outfile,sep="\t"):
    for line in infile:
        outfile.write(sep.join(replace_lst(line.strip().split(sep),mytable))+"\n")
    infile.close()
    outfile.close()

def main(argv):

    import argparse

    parser = argparse.ArgumentParser(description="Replace elements with the mapping file.")
    parser.add_argument('file1',nargs='?',help="file to be add information, \"-\" for stdin ")
    parser.add_argument('file2',nargs='?',help="dictfile")
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-i','--in-place',action='store_true',help="edit file in place.")
    parser.add_argument('-n','--header-only',action='store_true',help="Replace only the header")
    parser.add_argument('-k','--keypos',default=0,type=int,help="position of key for the dictfile(0-besed)[default 0]")
    parser.add_argument('-v','--valpos',default=1,type=int,help="position of value for the dictfile(0-besed)[default 1]")
    args = parser.parse_args(argv[1:])
    
    if args.file1 == '-':
        file1=sys.stdin
    else:
        file1=open(args.file1)
    file2=DictFile(args.file2,keypos=args.keypos, valuepos=args.valpos)
    outfile=args.outfile
    if args.in_place:
        outfile = NTmp('w+t',delete=False)
        inpath = file1.name
        outpath = outfile.name
    
    if args.header_only:
        replace_header(file1,file2,outfile)
    else:
        replace_file(file1,file2,outfile)

    if args.in_place:
        move(outpath,inpath)

if __name__ == '__main__':

    import sys
    main(sys.argv)
