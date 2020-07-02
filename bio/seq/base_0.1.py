#!/usr/bin/env python

import sys, re, string
from collections import OrderedDict

'''
License: GNU General Public License v3.0 (http://www.gnu.org/licenses/gpl-3.0.html)
Author: Mr. You Duan
Email: yduan@outlook.com

classes: Fasta Seq Gff Gff_rec
functions: write_fasta 
'''

class Seq(str):

    alphabet_nucleotide = set(["A","T","G","C","a","t","g","c"])
    alphabet_nucleotide_ext = set(["R","Y","M"])
    alphabet_nucleotide_ext.update(alphabet_nucleotide)

    def __new__(self,string,*args,**kwargs):

        return str.__new__(self,string)

    def __init__(self,string,name=None,*args,**kwargs):

        self.__name = name

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self,value):
        self.__name = value

    def __getattr__(self,name):

        if name == "__name":
            return self.__name
        else:
            return super(str,self).__getattr__(name)

    def getRC(self):

        '''
        get reverse compliment sequence
        '''

        intab = 'AaTtGgCc'
        outtab = 'TtAaCcGg'
        transtab = string.maketrans(intab,outtab)
        return Seq(self[::-1].translate(transtab),name=self.name)

    def isnucleotide(self,strict=True):
        if strict:
            _albet = Seq.alphabet_nucleotide
        else:
            _albet = Seq.alphabet_nucleotide_ext
        for _ in self:
            if _ not in _albet:
                return False
        return True

    def write2fasta(self,outfile,lth=80):
        if not self.__name:
            write_fasta_s(outfile,self,"seqNull",lth)
        else:
            write_fasta_o(outfile,self,lth)

class Fasta(object):
    
    '''
    data = {key:[seq]}
    key = geneID
    info = gene infomation
    val = sequence
    l_uplim, l_lowlim: up/low length limits of sequence
    '''

    def __init__(self,infile,mode=1,l_uplim=float("inf"),l_lowlim=-1):

        self.__length = 0
        self.filename, self.infile = open_file(infile)
        assert(isinstance(l_uplim,(int,float)))
        assert(isinstance(l_lowlim,(int,float)))
        self.l_uplim=l_uplim
        self.l_lowlim=l_lowlim
        if mode == 3:
            self.__data = {}
            # data2 keeps only the first field of name.
            self.__data2 = {}
            self.isiter = False
            self.loadfile(isdic=True)
        if mode == 2:
            self.__data = []
            self.isiter = False
            self.loadfile()
        if mode == 1:
            self.isiter = True#this object can be looked as an iterator
    
    def __len__(self):

        return self.__length

    #def __getitem__(self,key):

    #    return self.__data.get(key)
    def __getitem__(self,key):

        return self.__data[key]

    def __contains__(self,key):

        return key in self.__data

    def __iter__(self):

        #if not self.isiter:
            #should change to mode 1
        #    raise TypeError("Fasta mode should change to 1")

        self.fastaIter=self.iterFasta()
        return self

    def next(self):

        if not self.isiter:
            #should change to mode 1
            raise TypeError("Fasta mode should change to 1")
#        mygenerator = self.__iterFasta()
        seq = self.fastaIter.next()
        #length filter
        while len(seq)>self.l_uplim or len(seq)<self.l_lowlim:
            seq = self.fastaIter.next()
            continue
        #if seq:
        #    print 2
        return seq
        #else:
        #    return StopIteration()

    def iterFasta(self):
        #mode 1:yield each record as tuple(name,seq)
        #mode 2:return all data as dic

        name = ""
        seq = []
        for line in self.infile:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if not name:
                    seq = []
                    name = line.lstrip(">")
                    continue
                seq = "".join(seq)
                yield Seq(seq,name)
                seq = []
                name = line.lstrip(">")
            else:
                seq.append(line)
        seq = "".join(seq)
        yield Seq(seq,name)

    def getChr(self,key):

        if self.__data.get(key):
            return self.__data.get(key)
        else:
            return self.__data2.get(key)

    def getSeq(self,scaf,st,ed):
        scaf = self.getChr(scaf)
        if st > ed:
            st, ed = ed, st
            #return self.getRC(scaf[st-1: ed])
            return scaf[st-1: ed].getRC()
        else:
            return scaf[st-1:ed]

    def getNames(self):

        return self.__data.keys()

    def loadfile(self,isdic=False):

        for seq in self.iterFasta():
            if len(seq)>self.l_uplim or len(seq)<self.l_lowlim:
                continue
            if isdic:
                self.__data[seq.name] = seq
                self.__data2[seq.name.split()[0]] = seq
            else:
                self.__data.append(seq)
            self.__length += 1

    def getRC1(self,seq):
        '''
        abord
        get reverse compliment sequence
        '''

        intab = 'AaTtGgCc'
        outtab = 'TtAaCcGg'
        transtab = string.maketrans(intab,outtab)
        return seq[::-1].translate(transtab)

    def getRC(self,seq):
        return seq.getRC()
    
    def getData(self):

        return self.__data

class Gff(object):

    gff3_partern = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")

    def __init__(self,infile,mode=1,sep="\t"):

        self.sep = sep
        self.filename, self.infile = open_file(infile)
        
        if mode == 2:
            self.__data = []
            self.isiter = False
            self.loadfile()
        if mode == 1:
            self.isiter = True#this object can be looked as an iterator

    def __iter__(self):

        return self.iterator()

    def iterator(self):

        if not self.isiter:
            #should change to mode 1
            raise TypeError("Fasta mode should change to 1")
        return self.__iterGff()
    
    def next(self):

        return self.iterator().next()

    def __iterGff(self):

        for line in self.infile:
            line=line.strip()
            if line[0]=="#":
                continue
            yield Gff_rec(line.split(self.sep))

    def loadfile(self):
        
        #wait for compete
        raise TypeError("object Gff has not stock mode")
    
    @classmethod
    def Parse_attr(self,attr,attr_partern=gff3_partern):

        '''Parses Gff attribution string and results it as a dictionary
        '''

        #attr_dic = {}
        attr_dic = OrderedDict()
        #tmplst = map(lambda x:attr_partern.match(x),attr.strip().split(";")[:])
        tmplst = map(lambda x:attr_partern.match(x),attr.strip().split(";"))
        if not tmplst[-1]:
            tmplst.pop()#in some case,";" exists as the last charactor, which will cause a None item in tmplst[-1]
        for each in tmplst:
            if not each:
                print("gff type error")
                sys.exit()
            val = each.group(2)
            if val.startswith('"') and val.endswith('"'):
                val = val[1:-1]
            attr_dic.setdefault(each.group(1),val)
        return attr_dic

class Gff_rec(Seq):

    '''gff records
    '''
    def __init__(self,lst,chr=0,sour=1,type=2,start=3,end=4):
    
        self.chr = lst[chr]
        self.source = lst[sour]    
        self.type = lst[type]
        self.start = int(lst[start])
        self.end = int(lst[end])
        self.score = lst[5]
        self.strand = lst[6]
        self.phase = lst[7]
        self.attr = Gff.Parse_attr(lst[8])
        self.lst = lst
        self.seq = ''

    def __getattr__(self,name):

        raise KeyError("%s don\'t have an attribute named %s"% ("self.attr",name))
    
    def __str__(self,sep="\t"):

        return sep.join(self.lst)

    def __len__(self):

        return abs(int(self.end)-int(self.start))+1

    def str(self):

        return self.__str__()

    def get_lst(self):

        return [self.chr,self.source,self.type,self.start,self.end,\
            self.score,self.strand,self.phase,self.attr]

def open_file(infile):
    try:
        filename = infile
        fp = open(infile)
    except TypeError:
        filename = infile.name
        fp = infile
    return filename, fp

def write_fasta_o(outfile,seq,lth=80):
    #seq : Seq object
    #lenth: max sequence length in eachline
    write_fasta_s(outfile,seq,seq.name,lth)

def write_fasta_s(outfile,seq,name,lth=80):
    #lth: max sequence length in eachline

    outfile.write(">"+name+"\n")
    if lth==0:
        outfile.write("".join(seq)+"\n")
    elif lth<0 or not isinstance(lth,int):
        raise ValueError("length of sequence should be int and >=0")
    else:
        outfile.write("\n".join(seq[i*lth:i*lth+lth] for i in range((len(seq)-1)/lth+1))+"\n")

if __name__ == '__main__':

    import sys 
    myf = sys.argv[1]
    x=Fasta(myf,mode=2)
    print(x[-1].name)
