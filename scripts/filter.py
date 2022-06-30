#!/usr/bin/env python 

'''
[usage] ./filter file1 file2 outfile
    ./filter file1 file2 > outfile
[purpose]
    get lines in file1 which has key words in file2
'''

import sys

class CaseInsensitiveSet(set):
    """Make sense when element is tuple
    Refers to https://stackoverflow.com/questions/27531211/how-to-get-case-insensitive-python-set"""

    def get_lower(self, item):
        """Get lower for string or tuple"""
        if isinstance(item, str):
            return item.lower()
        elif isinstance(item, tuple):
            return tuple([self.get_lower(i) for i in item])
        else:
            return item

    def add(self, item):
        set.add(self, self.get_lower(item))

    def __contains__(self, item):
        return set.__contains__(self, self.get_lower(item))

def get_vals(lst,keycols):
    """ eg. get_vals([a,b,c,d], [0,3]) -> [a,d] """
    return tuple(map(lambda x:lst[x],keycols))

def readKeys(infile, sep, has_header=True, keycol=[0], ignorecase=False):
    if ignorecase:
        keys = CaseInsensitiveSet()
    else:
        keys = set()
    if has_header:
        infile.readline()
    for eachline in infile:
        tmp = eachline.strip("\n").split(sep)
        keys.add(get_vals(tmp, keycol))
    return keys

def _filter(infile,save_keys,sep,has_header,keycol,outfile,mode):
    """ Core function.
    mode:{filter,differ}
    mode_filter: return filter results
    mode_differ: return differ results
    """
    if has_header:
        outfile.write(infile.readline())
    #print(save_keys)
    for eachline in infile:
        tmp = eachline.strip("\n").split(sep)
        if (get_vals(tmp,keycol) in save_keys) and mode=="filter":
            outfile.write(eachline)
        elif (get_vals(tmp,keycol) not in save_keys) and mode=="differ":
            outfile.write(eachline)

def filter(infile,save_keys,sep,has_header,keycol,outfile):
    _filter(infile,save_keys,sep,has_header,keycol,outfile,"filter")

def differ(infile,save_keys,sep,has_header,keycol,outfile):
    _filter(infile,save_keys,sep,has_header,keycol,outfile,"differ")

def parse_keycol(keycol, sep=","):
    """ convert \"1,2,3,5\" into [1,2,3,5] """
    return list(map(int,keycol.strip().split(sep)))

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Filter file1 accoding to one file2. Two mode are avaliable.")
    parser.add_argument('file1',nargs='?',help="file to be filtered, \"-\" for stdin ")
    parser.add_argument('file2',nargs='?',help="file to guide the filtering, \"-\" for stdin(conflict with \"-\" for file1)")
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-s1','--sep1',nargs='?',default="\t")
    parser.add_argument('-s2','--sep2',nargs='?',default="\t")
    parser.add_argument('-k1','--keycol1',nargs='?',default="0")
    parser.add_argument('-k2','--keycol2',nargs='?',default="0")
    parser.add_argument('-n1','--has_header1',action='store_false',help="if file1 hasn't header,set this option")
    parser.add_argument('-n2','--has_header2',action='store_false',help='same with n1 but for file2')
    parser.add_argument('-i','--ignorecase',action='store_true',help='Ignore case')
    parser.add_argument('-m','--mode',nargs='?',default='filter',choices=['filter','differ'])
    args = parser.parse_args(argv[1:])

    sep1=args.sep1
    sep2=args.sep2
    keycol1=parse_keycol(args.keycol1)#key column which define filter limits
    keycol2=parse_keycol(args.keycol2)
    has_header1=args.has_header1
    has_header2=args.has_header2
    
    if args.file1 == '-' and args.file2 == '-':
        raise KeyError("You can not set file1 and file2 to \"-\" at the same time!")
    if args.file1 == '-':
        file1=sys.stdin
    else:
        file1=open(args.file1)
    if args.file2 == '-':
        file2=sys.stdin
    else:
        file2=open(args.file2)
    outfile=args.outfile

    keys = readKeys(file2,sep2,has_header2,keycol2,args.ignorecase)
    if args.mode == "filter":
        filter(file1,keys,sep1,has_header1,keycol1,outfile)
    elif args.mode == "differ":
        differ(file1,keys,sep1,has_header1,keycol1,outfile)

if __name__ == '__main__':

    main(sys.argv)
