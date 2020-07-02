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
    """
        Gff iterator.
    """

    gff3_partern = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")

    def __init__(self,infile,mode=1,sep="\t",fm="normal"):

        self.sep = sep
        self.fm = fm
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
            #yield Gff_rec(line.split(self.sep))
            yield Gff_rec(line.split(self.sep),self.fm)

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

class GffRec(Seq):

    def __init__(self, lst, fm):
        if fm == "normal":
            return Gff_rec(self, lst)

class Gff_rec(Seq):

    '''gff records
    '''
    def __init__(self,lst,fm,chr=0,sour=1,type=2,start=3,end=4):
    
        self.chr = lst[chr]
        self.Chr = self.chr
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
        self.format_attr(fm)

    def __getattr__(self,name):
        raise KeyError("%s don\'t have an attribute named %s"% ("self.attr",name))
    
    def __str__(self,sep="\t"):
        return sep.join(self.lst)

    def __len__(self):
        return abs(int(self.end)-int(self.start))+1

    def str(self):
        return self.__str__()

    def format_none(self):
        self.gene_id = ""
        self.tx_id = ""

    def format_normal(self):
        self.gene_id = self.attr["gene_id"]
        self.tx_id = self.attr["transcript_id"]

    # Grass carp gff format
    def format_gc(self):
        self.gene_id = self.attr["ID"].rstrip(".EXON")
        self.tx_id = self.gene_id

    def format_attr(self,fm):
        if fm == "normal":
            self.format_normal()
        elif fm == "gc":
            self.format_gc()
        elif fm == "none":
            self.format_none()
        else:
            raise KeyError("Wrong file format code: %s" % fm)

    def get_lst(self):
        return [self.chr,self.source,self.type,self.start,self.end,\
            self.score,self.strand,self.phase,self.attr]

# Genome element
class GENT(object):

    def __init__(self,start,end,strand="."):
        self.start = int(start)
        self.end = int(end)
        # if start end has been adjusted, in which end >= start.
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

    def __init__(self,*args,**kwargs):
        if isinstance(args[0],Gff_rec):
            o = args[0]
            self.init(o.start,o.end,o.strand,o.Chr,o.gene_id,o.tx_id)
        else:
            self.init(*args,**kwargs)

    def init(self,start,end,strand,Chr,gene_id,tx_id):
        self.Chr = Chr
        self.gene_id = gene_id
        self.tx_id = tx_id
        super(EXON,self).__init__(start,end,strand)

# segments genome element
# like: transcripts, genes.
# When parse gff record, basic information like start, end, strand will be load first.
# While comes across gtf file, basic infromation get through add_exon methods.
class SGENT(dict,GENT):

    def __new__(self,*args,**kwags):
        return dict.__new__(self)

    def __init__(self,es=[],perm_op=True,static=False,start=float("inf"),end=-float("inf"),strand=".",rec=None):
        # es: elements
        # perm_op: permit elements overlap
        # static: when static is True, start, end, strand, length option should not be update by add_exon.
        if rec == None:
            self._init_from_GENTs(es,perm_op,static,start,end,strand)
        elif isinstance(rec,Gff_rec):
            self._init_from_rec(rec,perm_op)
        else:
            sys.stderr.write("rec type is wrong\n")
            sys.exit(1)

    def __len__(self):
        return self.get_length()

    def _init_from_GENTs(self,es,perm_op,static,start,end,strand):
        self.start = start
        self.end = end
        self.strand = strand
        if static:
            self._init_static()
        else:
            self.static = False
            self._length = 0
        self._count = 0
        #self.extend(es)
        # Whether to permit elements overlap
        self.perm_op = perm_op

        # Add single element
        if type(es) != list:
            es = [es]
        self.add_es(es)

    def _init_static(self):
        self._adj_pos(self)
        self._length = self.start - self.end
        self.static = True

    def _init_from_rec(self,rec,perm_op):
        # Bugs!!!
        self.static = True
        self.start = rec.start
        self.end = rec.end
        self.strand = rec.strand
        # maybe only gc gff
        self.ID = rec.attr["ID"]
        self.name = rec.attr["Name"]
        self._count = 0
        self.perm_op = perm_op

    def add_es(self,es):
        for e in es:
            self.add_e(e)

    def _update_info(self,e):
        if self.static:
            sys.stderr.write("Wrong update calling.\n")
        self._length += len(e)
        self.start = min(self.start, e.start)
        self.end = max(self.end, e.end)

    def add_e(self,e):
        self.check_strand(e)
        self.check_overlap(e)
        self.add_data(e)
        self.add_len(e)
        if not self.static:
            self._update_info(e)

    def check_strand(self,e):
        if self.strand == ".":
            self.strand = e.strand
        if self.strand != e.strand:
            raise ValueError("strand of exon and transcript should be identical.\n")

    # Append
    def add_data(self,e, ID = None):
        if not ID:
            try:
                ID = e.ID
            except:
                ID = self.get_count()
        if ID in self.keys():
            raise KeyError("Reapeat ID, will cause cover issue.\n")
        self[ID] = e

    def add_len(self,e):
        # non-overlap
        self._count += 1

    def check_overlap(self,e2):
        if self.perm_op:
            return
        if check_overlap_e(self,e2):
            sys.stderr.write("%s, %s, %s\n"%(e2.Chr,e2.gene_id,e2.tx_id))
            sys.stderr.write("start1, end1, start2, end2 = %d,%d,%d,%d\n"%(self.start,self.end,e2.start,e2.end))
            raise TypeError("Overlapped elements are not permited.\n")

    def any_overlap(self):
        #print(self)
        tmp = sorted(self.values(),key = lambda x:x.start)
        #print(len(tmp))
        #print(self.get_count())
        for i in range(self.get_count()-1):
            if check_overlap_e(tmp[i],tmp[i+1]):
                return True
        return False

    # Get transcript start site.
    def get_tss(self):
        st = self.start
        ed = self.end
        if st > ed:
            st, ed = ed, st
        if self.strand == "-":
            return ed
        else:
            return st

    # Get transcript end site.
    def get_tes(self):
        st = self.start
        ed = self.end
        if st > ed:
            st, ed = ed, st
        if self.strand == "-":
            return st
        else:
            return ed

    def get_count(self):
        return self._count

    def get_length(self):
        return self._length

    def get_sorted_values(self,key = lambda x: x.start):
        return sorted(self.values(),key = key)

# transcript
class TxDict(SGENT):

    def __init__(self,*args,**kwargs):
        kwargs["perm_op"] = False
        super(TxDict,self).__init__(*args,**kwargs)

    def add_exon(self,exon):
        self.add_e(exon)

    def add_e(self,exon):
        assert isinstance(exon,EXON)
        self.ID = exon.tx_id
        self.Chr = exon.Chr
        self.gene_id = exon.gene_id
        super(TxDict,self).add_e(exon)

class GeneDict(SGENT):

    def add_exon(self,exon):
        if exon.tx_id in self.keys():
            # transcript add exon
            self[exon.tx_id].add_exon(exon)
            super(GeneDict,self)._update_info(exon)
        else:
            self.add_e(TxDict([exon]))

    def add_tx(self,tx):
        self.add_e(tx)

    def add_e(self,tx):
        assert isinstance(tx,TxDict)
        self.ID = tx.gene_id
        self.Chr = tx.Chr
        super(GeneDict,self).add_e(tx)

    # Get representative transcript of one gene, usally the longest one.
    # Useful when doing gene function annotation.
    def get_present_tx(self):
        pass

class ChrDict(SGENT):

    def __init__(self,es=[],perm_op=True):
        # es: elements
        # perm_op: permit elements overlap
        self.start = float("inf")
        self.end = - float("inf")
        self._length = 0
        self._count = 0
        #self.extend(es)
        # Whether to permit elements overlap
        self.perm_op = perm_op

        # Add single element
        if type(es) != list:
            es = [es]
        self.add_es(es)

    def add_exon(self,exon):
        if exon.gene_id in self.keys():
            self[exon.gene_id].add_exon(exon)
        else:
            self.add_e(GeneDict(TxDict([exon])))

    def add_tx(self,tx):
        if tx.gene_id in self.keys():
            self[tx.gene_id].add_tx(tx)
        else:
            self.add_e(GeneDict([tx]))

    def add_gene(self,gene):
        self.add_e(tx)

    def add_e(self,gene):
        assert isinstance(gene,GeneDict)
        self.ID = gene.Chr
        #super(ChrDict,self).add_e(gene)
        self.add_data(gene)
        self.add_len(gene)

    def sort(self):
        self.sort(key=lambda x:x.start)

    def __len__(self):
        return self.get_count()

class GenomeDict(SGENT):

    def add_exon(self,exon):
        if exon.Chr in self.keys():
            self[exon.Chr].add_exon(exon)
        else:
            self.add_e(ChrDict(GeneDict(TxDict([exon]))))

    def add_e(self,gene):
        assert isinstance(gene,ChrDict)
        super(GenomeDict,self).add_e(gene)

    def add_gene(self,gene):
        self.add_e(ChrDict(gene))

    def add_gene_rec(self,rec):
        gene = GeneDict(rec=rec)

    def check_strand(self,e):
        pass

    def __len__(self):
        return self.get_count()

class GffDict(GenomeDict):
    
    """ Read gff file, convert to GenomeDict data.
        Only contain exon, gene now. Others should be added.
    """
    def __init__(self,infile,sep="\t",keeps=set(["exon","gene"]),fm="normal"):
        super(GffDict,self).__init__([])
        self.loadfile(infile,sep,keeps,fm)

    def loadfile(self,infile,sep,keeps,fm):
        for rec in Gff(infile,sep=sep,fm=fm):
            # filter type
            if (keeps != None) and (rec.type not in keeps):
                continue
            if rec.type == "exon":
                self.add_exon(EXON(rec))
            elif rec.type == "gene":
                self.add_gene_rec(rec)

class GtfDict(GenomeDict):
    
    """ Read gtf file, convert to GenomeDict data.
        Only exon are used.
    """
    def __init__(self,infile,sep="\t",keeps=set(["exon"]),fm="normal"):
        super(GffDict,self).__init__([])
        self.loadfile(infile,sep,keeps,fm)

    def loadfile(self,infile,sep,keeps,fm):
        for rec in Gff(infile,sep=sep,fm=fm):
            # filter type
            #print(type(rec))
            if (keeps != None) and (rec.type not in keeps):
                continue
            self.add_exon(EXON(rec))

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

def check_overlap_e(e1,e2):
    st1,ed1 = e1.start, e1.end
    st2,ed2 = e2.start, e2.end
    return check_overlap(st1,ed1,st2,ed2)

# Check list of elements overlap
def check_overlap_l(data):
    data = sorted(data,key=lambda x:x.start)
    for i in range(len(data)-1):
        if check_overlap_e(data[i],data[i+1]):
            print("gene1,gene2 = %s, %s" % (data[i].ID,data[i+1].ID))
            print("s1,e1,s2,e2 = %d, %d, %d, %d" % (data[i].start,data[i].end,data[i+1].start,data[i+1].end))
            return True
    return False

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
