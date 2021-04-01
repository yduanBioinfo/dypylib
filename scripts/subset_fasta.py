#!/usr/bin/env python3 

'''
_version=0.0.2
get subset of a fasta file randomly
'''

import sys
import random
from dypylib.bio.seq.base import *
#from bio.seq.base import Fasta, Falist, write_fasta_o

# overflow(true): If sample is large than population, return whole population.
def random_fa(seqs,outfile,output_num,overflow,lth):
    if overflow and len(seqs) < output_num:
        subset = seqs
    # return whole set when output num set to -1
    elif output_num == -1:
        subset = seqs
    else: 
        subset = random.sample(seqs,output_num)

    for seq in subset:
        write_fasta_o(outfile,seq,lth)

def random_fa_with_filter(seqs,keeps,outfile,output_num,overflow=False,filt_type="ref",lth=80,if_verse=False):
    """
    #ref_filt: read fasta head line in refSef style,">ref|refname|xxbd|xxname|" and only keep refname or xxname.
    if_verse: Set it to choose sequence outside keeps list rather than inside.
    """

    filt_seqs = []
    for seq in seqs:
        is_in = filter_keys(seq.name,keeps,filt_type)
        if (not if_verse) and (not is_in):
            continue
        if if_verse and is_in:
            continue
        #if seq.name in keeps:
        filt_seqs.append(seq)
    random_fa(filt_seqs,outfile,output_num,overflow,lth)

def read_keyf(key_file,col=0,sep="\t"):

    outdata = set()
    for line in key_file:
        tmp = line.strip().split(sep)
        if not tmp:
            continue
        outdata.add(tmp[col])
    return outdata

def filter_keys(myname,keeps,mytype="ref"):

    if mytype == "ref":
        return filter_ref_key(get_id_dic(myname),keeps)
    elif mytype == "GenCode":
        return filter_gcd_key(myname,keeps)
    elif mytype == "basic":
        return filter_bsc_key(myname,keeps)

def filter_ref_key(id_dic,keeps,key="ref"):

    myid = id_dic.get(key)
    if not myid:
        return False
    if myid in keeps:
        return True
    else:
        return False

def filter_gcd_key():

    pass

def filter_bsc_key(myname,keeps):

    if myname in keeps:
        return True
    else:
        return False

def get_id_dic(name,sep="|"):
    #get id information from refseq seq name

    tmp = name.split(sep)[:-1]
    key_val = {}
    for i in range(0,len(tmp),2):
        key_val[tmp[i]] = tmp[i+1]
    return key_val

def main(argv):

    import argparse

    parser = argparse.ArgumentParser(description="""
    Get subset of a fasta file randomly (for rand_num > 0)
            
    Example:
        1. Subset fasta file and keep sequence who's ID is in 
        the namelist.
            subset_fasta.py input.fa -k namelist -t basic
        2. Subset fasta file and keep sequence who's ID is not
        in the namelist.
            subset_fasta.py input.fa -k namelist -t basic -v
        3. Random select 100 sequences from fasta file.
            subset_fasta.py input.fa -n 100 -t basic
        4. Random select 100 sequences from fasta file and who's
        ID are in the namelist.
            subset_fasta.py input.fa -k namelist -n 100 -t basic
    """,
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('infile',nargs='?',help="fasta file for sampling ",default=sys.stdin,type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-k','--key_file',nargs='?',help="if key_file provided, only sequence in key_file should be keep",type=argparse.FileType('r'))
    parser.add_argument('-r','--if-trim',action='store_true',\
        help = ("Whether to trim fasta ID by [SPACE]? Set for yes.\n"
            "Note: this option will be replaced by --trim-name in "
            "the upcoming version."))
    parser.add_argument('--trim-name',help="Trim fasta ID by [SPACE]",action='store_true')
    parser.add_argument('-v','--if-verse',help="Choose sequences outside the key_file rather than inside the key_file",action='store_true')
    parser.add_argument('-t','--filt_type',nargs='?',help="filter type,can be \"ref\", \"basic\" or \"GenCode\"(not avaliable now). The default should be altered to basic.",default="ref")
    parser.add_argument('-l','--lowlim',nargs='?',help="low limits of length",type=int,default=-1)
    parser.add_argument('--overflow',dest='over',help="To allow sample be large than population(return whole fasta).",action='store_true')
    parser.add_argument('--row_len',nargs='?',help="sequence length of each row",type=int,default=80)
    parser.add_argument('-u','--uplim',nargs='?',help="up limits of length",type=int,default=float("inf"))
    parser.add_argument('-n','--rand_num',nargs='?',help="how many sequence should be random get",type=int,default=-1)
    args = parser.parse_args(argv[1:])

    if args.key_file:
        keys = read_keyf(args.key_file)
        #print(keys)
        myfa = Fasta(args.infile,l_lowlim=args.lowlim,l_uplim=args.uplim,if_trim=(args.if_trim or args.trim_name))
        #for key in myfa:
        #    print(key.name)
        #sys.exit()
        filt_type = args.filt_type
        random_fa_with_filter(myfa,keys,args.outfile,args.rand_num,args.over,filt_type,lth=args.row_len,if_verse = args.if_verse)#have option ref_filt
    else:
        myfa = Falist(args.infile,l_lowlim=args.lowlim,l_uplim=args.uplim,if_trim=(args.if_trim or args.trim_name)) 
        random_fa(myfa,args.outfile,args.rand_num,args.over,args.row_len)

if __name__ == '__main__':
    
    import sys
    sys.exit(main(sys.argv))
