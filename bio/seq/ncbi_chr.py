#!/usr/bin/env python

NCBI2CHR_hs38={"NC_000010.11":"chr10","NC_000011.10":"chr11","NC_000012.12":"chr12","NC_000013.11":"chr13","NC_000014.9":"chr14","NC_000015.10":"chr15","NC_000016.10":"chr16","NC_000017.11":"chr17","NC_000018.10":"chr18","NC_000019.10":"chr19","NC_000001.11":"chr1","NC_000020.11":"chr20","NC_000021.9":"chr21","NC_000022.11":"chr22","NC_000002.12":"chr2","NC_000003.12":"chr3","NC_000004.12":"chr4","NC_000005.10":"chr5","NC_000006.12":"chr6","NC_000007.14":"chr7","NC_000008.11":"chr8","NC_000009.12":"chr9","NC_012920.1":"chrM","NC_000023.11":"chrX","NC_000024.10":"chrY"}
NCBI2CHR_mm10={"NC_000076.6":"chr10","NC_000077.6":"chr11","NC_000078.6":"chr12","NC_000079.6":"chr13","NC_000080.6":"chr14","NC_000081.6":"chr15","NC_000082.6":"chr16","NC_000083.6":"chr17","NC_000084.6":"chr18","NC_000085.6":"chr19","NC_000067.6":"chr1","NC_000068.7":"chr2","NC_000069.6":"chr3","NC_000070.6":"chr4","NC_000071.6":"chr5","NC_000072.6":"chr6","NC_000073.6":"chr7","NC_000074.6":"chr8","NC_000075.6":"chr9","NC_005089.1":"chrMT","NC_000086.7":"chrX","NC_000087.7":"chrY",}
NCBI2CHR_Gz10={"NC_007121.6":"chr10","NC_007122.6":"chr11","NC_007123.6":"chr12","NC_007124.6":"chr13","NC_007125.6":"chr14","NC_007126.6":"chr15","NC_007127.6":"chr16","NC_007128.6":"chr17","NC_007129.6":"chr18","NC_007130.6":"chr19","NC_007112.6":"chr1","NC_007131.6":"chr20","NC_007132.6":"chr21","NC_007133.6":"chr22","NC_007134.6":"chr23","NC_007135.6":"chr24","NC_007136.6":"chr25","NC_007113.6":"chr2","NC_007114.6":"chr3","NC_007115.6":"chr4","NC_007116.6":"chr5","NC_007117.6":"chr6","NC_007118.6":"chr7","NC_007119.6":"chr8","NC_007120.6":"chr9","NC_002333.2":"chrM"}#GRCz10

def gff_ncbi2chr(infile,outfile,nmDict,sep="\t",filt_chr=False):
    
    try:infile=open(infile)
    except:pass
    try:outfile=open(outfile,'w')
    except:pass

    for line in infile:
        tmp = line.strip().split(sep)
        _ = nmDict.get(tmp[0])
        if filt_chr and not _:
            continue
        tmp[0] = nmDict.get(tmp[0]) or tmp[0]
        outfile.write(sep.join(tmp)+"\n")

def fa_ncbi2chr(infile,outfile,nmDict,sep="\t"):
    #used once. for further ultilize some update must be applied
    try:infile=open(infile)
    except:pass
    try:outfile=open(outfile,'w')
    except:pass

    for line in infile:
        if line[0] == ">":
            line_array=line.strip().split()
            name=line_array[0][1:]
            name=nmDict.get(name) or name
            outfile.write(">"+name+"\n")
        else:
            outfile.write(line)

def fa_ncbi2chr2(infile,outfile,nmDict,sep="\t"):
    #used once.for refSeq genome as gi|xxx|ref|xxxxx| other information
    #for further ultilize some update must be applied
    try:infile=open(infile)
    except:pass
    try:outfile=open(outfile,'w')
    except:pass

    for line in infile:
        if line[0] == ">":
            line_array=line.strip().split()
            name=line_array[0][1:]
            ref_name_array = name.split("|")
            ref_name = ref_name_array[3]
            name=nmDict.get(ref_name) or name
            outfile.write(">"+name+"\n")
        else:
            outfile.write(line)
#import sys
#gff_ncbi2chr(sys.argv[1],sys.stdout)

def main(argv):
    
    import argparse, sys

    parser = argparse.ArgumentParser(description="transform chromsome id from NCBI accssion number to chr* in gtf file")
    parser.add_argument('infile',nargs='?',default=sys.stdin,help="gtf/gff file",type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('--sep',nargs='?',default="\t",help="default tab")
    parser.add_argument('--strict',action='store_true',help="filt out non-chromosome sequence")
    parser.add_argument('-s','--species',nargs='?',default="hs38",help="genome version",choices=["hs38","mm10","Gz10"])
    args = parser.parse_args(argv[1:])

    if len(argv) == 1:
        parser.print_help()
        sys.exit(1)
    if args.species == "hs38":
        gff_ncbi2chr(args.infile,args.outfile,NCBI2CHR_hs38,filt_chr=args.strict)
        #fa_ncbi2chr2(args.infile,args.outfile,NCBI2CHR_hs38)
    elif args.species == "mm10":
        gff_ncbi2chr(args.infile,args.outfile,NCBI2CHR_mm10,filt_chr=args.strict)
    elif args.species == "Gz10":
        #fa_ncbi2chr(args.infile,args.outfile,NCBI2CHR_Gz10)
        gff_ncbi2chr(args.infile,args.outfile,NCBI2CHR_Gz10,filt_chr=args.strict)
    else:
        raise NameError("you specify one unsupport genome version")

if __name__ == '__main__':

    import sys
    sys.exit(main(sys.argv))
