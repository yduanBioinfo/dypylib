#!/usr/bin/env python

import sys, re, string, itertools
from collections import OrderedDict

'''
License: GNU General Public License v3.0 (http://www.gnu.org/licenses/gpl-3.0.html)
Author: Mr. You Duan
Email: yduan@outlook.com

classes: Fasta Seq Gff Gff_rec
functions: write_fasta 
'''

vs=sys.version

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

    def __getitem__(self,key):
        return Seq(super(Seq,self).__getitem__(key),self.name)
    # python version 2.x
    def __getslice__(self,*args):
        return Seq(super(Seq,self).__getslice__(*args),self.name)

    def getRC(self):

        '''
        get reverse compliment sequence
        '''

        intab = 'AaTtGgCc'
        outtab = 'TtAaCcGg'
        if vs < '3.0':
            transtab = string.maketrans(intab,outtab)
        else:
            transtab = str.maketrans(intab,outtab)
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

class Fadict(dict):

    def __new__(self,*args,**kwargs):
        return dict.__new__(self)

    def __init__(self,infile,strict=False,**kwargs):
        faiter = Fasta(infile,**kwargs)
        self.filename, self.l_uplim, self.l_lowlim = \
        faiter.filename, faiter.l_uplim, faiter.l_lowlim
        for seq in faiter:
            if strict and seq.name in self.keys():
                raise KeyError("Double key: %s has being found. Repick your file or set strict to \"False\"" % seq.name)
            self[seq.name] = seq

    def getChr(self,key):
        return self.get(key)

    def getSeq(self,scaf,st,ed):
        scaf = self.getChr(scaf)
        if st > ed:
            st, ed = ed, st
            return scaf[st-1: ed].getRC()
        else:
            return scaf[st-1:ed]

    def getNames(self):
        return self.keys()

class Falist(list):

    def __new__(self,*args,**kwargs):
        return list.__new__(self)

    def __init__(self,*args,**kwargs):
        faiter = Fasta(*args,**kwargs)
        self.filename, self.l_uplim, self.l_lowlim =\
        faiter.filename, faiter.l_uplim, faiter.l_lowlim
        for seq in faiter:
            self.append(seq)

class Fasta(object):
    
    '''
    l_uplim, l_lowlim: up/low length limits of sequence
    '''

    def __init__(self,infile,l_uplim=float("inf"),l_lowlim=-1):
        self.filename, self.infile = open_file(infile)
        assert(isinstance(l_uplim,(int,float)))
        assert(isinstance(l_lowlim,(int,float)))
        self.l_uplim=l_uplim
        self.l_lowlim=l_lowlim

    def __iter__(self):
        self.fastaIter=self.iterFasta()
        return self

    # python2 iterator
    def next(self):
        seq = self.fastaIter.next()
        #length filter
        while len(seq)>self.l_uplim or len(seq)<self.l_lowlim:
            seq = self.fastaIter.next()
            continue
        return seq

    # python3 iterator
    def __next__(self):
        seq = self.fastaIter.__next__()
        #length filter
        while len(seq)>self.l_uplim or len(seq)<self.l_lowlim:
            seq = self.fastaIter.__next__()
            continue
        return seq
        #return self.next()

    def iterFasta(self):
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
        if vs > '3.0':
            tmplst = list(tmplst)
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

# Genome element
class GENT(object):

    def __init__(self,start,end,strand="."):
        self.start = int(start)
        self.end = int(end)
        # if start end has been adjusted, in which end >= start.
        #self.adjsed = False
        self.strand = strand
        self._guess_strand()
        self._adj_pos()

    def __len__(self):
        return abs(self.end - self.start) + 1

    def _guess_strand(self):
        if self.start > self.end:
            if self.strand == "+":
                raise ValueError("Stand is not suit for coordinates.\n")
            self.strand = "-"

    def _adj_pos(self):
        if self.start > self.end:
            self.start, self.end = self.end, self.start

class EXON(GENT):

    pass

# segments genome element
# like: transcripts, genes.
class SGENT(GENT):

    def __init__(self,es=[],perm_op=True):
        # es: elements
        self.start = float("inf")
        self.end = - float("inf")
        self.strand = "."
        self.__length = 0
        self.__data = []
        self.__data.extend(es)
        # Whether to permit elements overlap
        self.perm_op = perm_op

        self.add_es(es)
        #self.adjsed = True

    def add_es(self,es):
        if not self.perm_op:
            es=sorted(es,key=lambda x:x.start)
        for e in es:
            self.add_e(e)

    def add_e(self,e):
        if self.strand == ".":
            self.strand = e.strand
        if self.strand != e.strand:
            raise ValueError("strand of exon and transcript should be identical.\n")
        self.check_overlap(e)
        self.start = min(self.start, e.start)
        self.end = max(self.end, e.end)

        self.__length += len(e)

    def check_overlap(self,e2):
        if self.perm_op:
            return
        st1,ed1 = self.start, self.end
        st2,ed2 = e2.start, e2.end
        if check_overlap(st1,ed1,st2,ed2):
            raise TypeError("Overlapped elements are not permited.\n")

    def __len__(self):
        return self.__length

    def __iter__(self):
        self.iterator=self.iterator()
        return self

    # python2 iterator
    def next(self):
        return self.iterator.next()

    # python3 iterator
    def __next__(self):
        return self.iterator.__next__()

    def iterator(self):
        return iter(self.__data)

# transcript
class Tx(SGENT):

    def __init__(self,*args,**kwargs):
        kwargs["perm_op"] = False
        super(Tx,self).__init__(*args,**kwargs)

    def add_exon(self,exon):
        self.add_e(exon)

    def add_e(self,exon):
        assert isinstance(exon,EXON)
        super(Tx,self).add_e(exon)

    def get_exons(self):
        return self.__data

class Gene(SGENT):

    def add_gene(self,tx):
        self.add_e(tx)

    def add_e(self,tx):
        assert isinstance(tx,Tx)
        super(Gene,self).add_e(tx)

    def get_txs(self):
        return self.__data

    # Get representative transcript of one gene, usally the longest one.
    # Useful when doing gene function annotation.
    def get_present_tx(self):
        pass

class GeneList(SGENT):

#    def __new__(self,*args,**kwags):
#        return list.__new__(self)

    def add_gene(self,gene):
        self.add_e(tx)

    def add_e(self,gene):
        assert isinstance(gene,Gene)
        super(GeneList,self).add_e(gene)

    def getGenes(self):
        return self.__data

    def sort(self):
        self.__data.sort(key=lambda x:x.start)

def check_overlap(st1,ed1,st2,ed2):
    def co(b1,a2):
        if b1 >= a2:
            return True
        else:
            return False
    if st1 < st2:
        return co(ed1,st2)
    else:
        return co(ed2,st1)

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
        if vs < '3.0':
            outfile.write("\n".join(seq[i*lth:i*lth+lth] for i in range((len(seq)-1)/lth+1))+"\n")
        else:
            outfile.write("\n".join(seq[i*lth:i*lth+lth] for i in range(int((len(seq)-1)/lth)+1))+"\n")
        
if __name__ == '__main__':

    import sys 
    myf = sys.argv[1]
    x=Fasta(myf,mode=2)
    print(x[-1].name)
