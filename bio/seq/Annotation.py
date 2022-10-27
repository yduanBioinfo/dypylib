#!/usr/bin/env python

'''
License: GNU General Public License v3.0 (http://www.gnu.org/licenses/gpl-3.0.html)
Author: Mr. You Duan
Email: yduan@outlook.com

classes: Fasta Seq Gff Gff_rec
functions: write_fasta 
'''

import sys, re, string, itertools, copy
from abc import ABC
from _collections_abc import Mapping
import gffutils, sqlite3
from collections import OrderedDict as Ordic
try:
    from pygtftk.gtf_interface import GTF as tkGTF
    from pygtftk.Line import Feature as tkFeature
    import pysnooper
except:
    pass

try:
    from case_insensitive_dict import CaseInsensitiveDict
except:
    pass

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

    def get_freq(self,alpha="g"):
        """ Get the frequency of the alpha.
        """
        ts = self.lower()
        alpha = alpha.lower()
        freq = ts.count(alpha)/len(ts)
        return freq

    def get_gc(self):
        gc = self.get_freq("g") + self.get_freq("c")
        return gc

    def to_fasta(self,outfile,lth=80):
        if not self.__name:
            write_fasta_s(outfile,self,"seqNull",lth)
        else:
            write_fasta_o(outfile,self,lth)

    def write2fasta(self,outfile,lth=80):
        sys.stderr.write("write2fasta will retired, use to_fasta replace\n")
        self.to_fasta(outfile,lth)

class Gff(object):
    """
        Gff iterator.
    """

    gff3_partern = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")

    def __init__(self,infile,mode=1,sep="\t",fm="normal",check_sort=False):
        self.sep = sep
        self.fm = fm
        self.filename, self.infile = open_file(infile)
        # Use to check duplicated gene/tx.
        self.check_sort = check_sort
        if self.check_sort:
            self.itered_tx = set()
            self.curr_tx = ''
        
        if mode == 2:
            self.__data = []
            self.isiter = False
            self.loadfile()
        if mode == 1:
            self.isiter = True#this object can be looked as an iterator

    def __iter__(self):
        return self.iterator()

    def check_exist(self, name, name_set):
        if name in name_set:
            return True
        name_set.add(name)
        return False

    def if_gff_sorted(self, rec):
        # Tx_id information is not available.
        if self.fm == 'none':
            return 
        if rec.tx_id == self.curr_tx:
            return
        if self.check_exist(rec.tx_id, self.itered_tx):
            # gff is not sorted
            error_info = ('Tx: "{}" is not continuous '
                'in file\n').format(rec.tx_id)
            sort_help=('Input gff/gtf is not sorted\nYou can sort gtf'
                'with the following order:\ncat input.gtf|cgat'
                'gtf2gtf --method=sort --sort-order "'
                'gene+transcript" > sorted.gtf\n')
            error_info += sort_help
            raise IndexError(error_info)
        else:
            # Meet a new tx.
            # Update current transcript.
            self.curr_tx = rec.tx_id

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
            if len(line)==0 or line[0]=="#":
                continue
            #yield Gff_rec(line.split(self.sep))
            rec = Gff_rec(line.split(self.sep),self.fm)
            if self.check_sort:
                self.if_gff_sorted(rec)
            yield rec

    def loadfile(self):
        #wait for compete
        raise TypeError("object Gff has not stock mode")
    
    @classmethod
    def Parse_attr(self,attr,attr_partern=gff3_partern):
        '''Parses Gff attribution string and results it as a dictionary
        '''

        #attr_dic = {}
        attr_dic = Ordic()
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

    #def __getattr__(self,name):
    #    raise KeyError("%s don\'t have an attribute named %s"% ("self.attr",name))
    
    def __str__(self,sep="\t",version=2):
        data = self.get_lst()
        data[-1] = self.package_attr(version)
        return sep.join(map(str,data))

    def __len__(self):
        return abs(int(self.end)-int(self.start))+1

    def as_str(self,sep="\t",version=2):
        return self.__str__(sep,version)

    def format_none(self):
        self.gene_id = ""
        self.tx_id = ""

    def format_normal(self):
        # gene records don't have gene_id attribution
        self.gene_id = self.attr.get("gene_id")
        # gene records don't have transcript_id attribution
        self.tx_id = self.attr.get("transcript_id")

    # Grass carp gff format
    def format_gc(self):
        # Warning: rstrip should be replaced.
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

    # Refers to Ensembl GTF
    def pack_attr_gtf(self):
        data = []
        for key, val in self.attr.items():
            data.append(key+" \""+val+"\"")
        return "; ".join(data)+";"

    def edit_gene_id(self,val):
        self.gene_id = val
        # Raise unknown results when format is not normal.
        self.attr["gene_id"] = val

    def edit_tx_id(self,val):
        self.tx_id = val
        # Raise unknown results when format is not normal.
        self.attr["transcript_id"] = val

    # Refers to Ensembl GFF
    def pack_attr_gff(self):
        data = []
        for key, val in self.attr.items():
            data.append(key+"="+val)
        return "; ".join(data)

    def package_attr(self,gff_version=2):
        if gff_version == 2:
            return self.pack_attr_gtf()
        elif gff_version == 3:
            return self.pack_attr_gff()
        else:
            raise KeyError("Version can be only 2 or 3.\n")

    def get_lst(self):
        return [self.chr,self.source,self.type,self.start,self.end,\
            self.score,self.strand,self.phase,self.attr]

class GENT(object):
    """ Genomic element.
        To-do:
        Try:
        class GffutilsGENT(object, metaclass=Test_main):
        TypeError: metaclass conflict: the metaclass of a derived class must be a (non-strict) subclass of the metaclasses of all its bases
    """

    def __len__(self):
        return abs(self.end - self.start) + 1

# Genome element
class dyGENT(GENT):
    """ Genomic element for dypylib
    """

    def __init__(self,start,end,strand=".",attr={}):
        self.engine = "dypylib"
        self.start = int(start)
        self.end = int(end)
        # if start end has been adjusted, in which end >= start.
        self.strand = strand
        self._guess_strand()
        self._adj_pos()
        self.attr=attr

    def _guess_strand(self):
        if self.start > self.end:
            if self.strand == "+":
                raise ValueError("Stand is not suit for coordinates.\n")
            self.strand = "-"

    def _adj_pos(self):
        if self.start > self.end:
            self.start, self.end = self.end, self.start

class GffutilsGENT(GENT):
    """ Genomic element for Gffutils
    """
    def __init__(self, *args, **kwargs):
        """Should make a function to ensure self.db = self.db,
        self.engine = self.engine when create from GffutilsGENT
        """
        self.engine = default_param(kwargs, 'engine', 'dypylib')
        if self.engine == 'gffutils':
            self._gffutils__init__(*args, **kwargs)
        elif self.engine == 'dypylib':
            self._dy__init__(*args, **kwargs)

    def __getattr__(self, attr):
        # Try to return attributs definded in database.
        if self.engine == 'gffutils':
            return self._gffutils__getattr__(attr)

    def _gffutils__init__(self, db, *args, name = "",\
            **kwargs):
        self.db = db
        self.name = name

    def _gffutils__getattr__(self, attr):
        # Try to return attributs definded in database.
        return getattr(self.db[self.name], attr)

    def _guess_strand(self):
        """
        """
        if self.start > self.end:
            if self.strand == "+":
                raise ValueError("Stand is not suit for coordinates.\n")
            self.strand = "-"

    def _adj_pos(self):
        """dy methods, Warning: gffutils might be affected.
        """
        if self.start > self.end:
            self.start, self.end = self.end, self.start

class BaseExon(object):
    """Base class of Exon"""

class dyExon(dyGENT, BaseExon):

    def __init__(self,*args,**kwargs):
        """ Bugs:
            When the attibution has been altered, the correponding value in self.rec keeps original.
            Need a update_rec function to solve this problem.
        """
        if isinstance(args[0],Gff_rec):
            o = args[0]
            self.init(o.start,o.end,o.strand,o.Chr,o.gene_id,o.tx_id)
            self.rec = o
        else:
            self.init(*args,**kwargs)
            self.rec = None

    def init(self,start,end,strand,Chr,gene_id,tx_id):
        self.Chr = Chr
        self.gene_id = gene_id
        self.tx_id = tx_id
        super(dyExon,self).__init__(start,end,strand)

    def edit_gene_id(self,val):
        """ Alter gene_id
        """
        self.gene_id = val
        self.rec.edit_gene_id(val)

    def edit_tx_id(self,val):
        """ Alter tx_id
        """
        self.tx_id = val
        self.rec.edit_tx_id(val)

    def update_rec(self):
        """ Under development...
        """
        pass

    def as_str(self):
        return self.rec.as_str()

    def as_gtf(self):
        return self.as_str()

class GffutilsExon(GffutilsGENT, BaseExon):
    pass

class BaseCDS(object):
    """Base class of CDS"""

class dyCDS(dyGENT, BaseCDS):
    pass

class GffutilsCDS(GffutilsGENT, BaseCDS):
    pass

class BaseIntron(object):
    """Base class of CDS"""

class dyIntron(dyGENT, BaseIntron):
    """
    """

    def __init__(self,start,end,strand,Chr,gene_id = None,tx_id = None):
        self.Chr = Chr
        self.chrom = Chr
        self.gene_id = gene_id
        self.tx_id = tx_id
        super(dyIntron,self).__init__(start,end,strand)

    def __hash__(self):
        return hash((self.start, self.end, self.strand, self.Chr, \
                self.gene_id, self.tx_id))

    def __eq__(self, other):
        try:
            return (self.Chr == other.Chr) and (self.strand == other.strand) and \
                   (self.start == other.start) and (self.end == other.end)
        except:
            return False

class GffutilsIntron(GffutilsGENT, BaseIntron):
    pass

class dySGENT(dict,dyGENT):
    """Segments genome element
    like: transcripts, genes.
    When parse gff record, basic information like start, end, strand will be load first.
    While comes across gtf file, basic infromation was acquired through add_exon methods.
    """
    def __new__(self,*args,**kwargs):
        return dict.__new__(self)

    def __init__(self,es=[],perm_op=True,static=False,start=float("inf"),end=-float("inf"),strand=".",rec=None,engine="dypylib"):
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

    def _init_from_gffutils(self, db, key):
        """Init Object from a gffutils database"""
        pass

    def _init_from_rec(self,rec,perm_op):
        # Bugs!!!
        self.static = True
        self.start = rec.start
        self.end = rec.end
        self.strand = rec.strand
        # maybe only gc gff
        self.ID = rec.attr.get("ID")
        self.name = rec.attr.get("Name")
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
        tmp = sorted(self.values(),key = lambda x:x.start)
        for i in range(self.get_count()-1):
            if check_overlap_e(tmp[i],tmp[i+1]):
                return True
        return False

    # Get transcript start site.
    def get_tss(self):
        return self.get_5p_end()

    def get_tes(self):
        return self.get_3p_end()

    def get_5p_end(self):
        return self._get_tss()

    def get_3p_end(self):
        return self._get_tes()

    def _get_tss(self):
        st = self.start
        ed = self.end
        if st > ed:
            st, ed = ed, st
        if self.strand == "-":
            return ed
        else:
            return st

    # Get transcript end site.
    def _get_tes(self):
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

class GffutilsSGENT(Mapping, GffutilsGENT):
    """Mapping object for biology purpose under Gffutils engine.
    Data are deposited in database.

        get_parents -> get_parents_id like get_children_id
    """
    def __new__(self,*args,engine = "gffutils",**kwargs):
        return Mapping.__new__(self)

    def __init__(self, *args, engine="gffutils", **kwargs):
        super(GffutilsSGENT, self).__init__(*args,engine=engine, **kwargs)

    def __getitem__(self, key):
        if self.engine == 'gffutils':
            return self._gffutils__getitem__(key)

    def __contains__(self, key):
        return key in self.children_id

    def _gffutils__getattr__(self, attr):
        """Don't init children and parents in __int__,
        rather find it from datbase when is needed.

        This option can avoid loading every relation in
            python once construct an object.
        """
        if attr == 'children':
            self.children = self.get_children()
            return self.children
        elif attr == 'parents':
            self.children = self.get_parents()
            return self.children
        elif attr == 'children_id':
            self.children_id = self.get_children_id()
            return self.children_id
        else:
            return getattr(self.db[self.name], attr)

    def __iter__(self):
        """Self.keys() are generated here."""
        keys = self.children_id
        return iter(keys)

    def __len__(self):
        return len(self.children_id)

    def get_children_id(self):
        if self.engine == 'gffutils':
            return self._gffutils_get_children_id()
        else:
            raise TypeError("Only gffutils have get_children_id function")

    def _gffutils_get_children_id(self):
        children = self.db.children(self.name)
        return [i.id for i in children]

    def get_children(self):
        return self._gffutils_get_children()

    def _gffutils_get_children(self):
        return self.db.children(self.name)

    def get_parents(self):
        return self._gffutils_get_parents()

    def _gffutils_get_parents(self):
        return self.db.parents(self.name)

    def _gffutils__getitem__(self, key):
        if key not in self:
            info = "{} is not contained in {}, please check it.\
                    ".format(key, self.name)
            raise Keself.dbyError(info)
        return self.db[key]

class BaseTranscript(object):
    """Base class of Transcript"""

# transcript
class TxDict(dySGENT, BaseTranscript):

    def __init__(self,*args,engine="dypylib",**kwargs):
        kwargs["perm_op"] = False
        super(TxDict,self).__init__(*args,**kwargs)

    def add_exon(self,exon):
        self.add_e(exon)

    def add_e(self,exon):
        assert isinstance(exon,Exon)
        self.ID = exon.tx_id
        self.Chr = exon.Chr
        self.gene_id = exon.gene_id
        super(TxDict,self).add_e(exon)

    def edit_gene_id(self,val):
        """ Alter gene_id
        """
        self.gene_id = val
        for exon in self.values():
            exon.edit_gene_id(val)

    def edit_tx_id(self,val):
        """ Alter tx_id
        """
        self.ID = val
        for exon in self.values():
            exon.edit_tx_id(val)

    def as_gtf(self):
        """ Get string in GTF format.
        """
        return "\n".join(map(lambda x:x.as_gtf(),self.values()))

    def get_introns(self):
        """ Return a list of introns of this transcript.
        """
        res = []
        exons = sorted(list(self.values()),key=lambda x:x.start)
        for i in range(len(exons)-1):
            e_0 = exons[i]
            e_1 = exons[i+1]
            _start = e_0.end + 1
            _end = e_1.start - 1
            if _start > _end:
                raise ValueError("The end of the intron must be bigger than the start")
            res.append(Intron(_start,_end,e_0.strand,e_0.Chr,e_0.gene_id,e_0.tx_id))
        return res

class GffutilsTranscript(GffutilsSGENT, BaseTranscript):
    """Chromosome object
    """

    def _gffutils__getitem__(self, key):
        if key not in self:
            info = "{} is not contained in {}, please check it.\
                    ".format(key, self.name)
            raise KeyError(info)
        feature = self.db[key]
        if feature.featuretype == 'exon':
            return Exon(self.db, name=feature.id, engine='gffutils')
        elif feature.featuretype == 'cds':
            return CDS(self.db, name=feature.id, engine='gffutils')
        else:
            return GffutilsGENT(self.db, name=feature.id, engine='gffutils')

class GeneDict(dySGENT):

    def add_exon(self,exon):
        if exon.tx_id in self.keys():
            # transcript add exon
            self[exon.tx_id].add_exon(exon)
            super(GeneDict,self)._update_info(exon)
        else:
            self.add_e(Transcript([exon]))

    def add_tx(self,tx):
        self.add_e(tx)

    def add_e(self,tx):
        assert isinstance(tx,Transcript)
        self.ID = tx.gene_id
        self.Chr = tx.Chr
        super(GeneDict,self).add_e(tx)

    def edit_gene_id(self,val):
        """ Alter gene_id
        """
        self.ID = val
        for tx in self.values():
            tx.edit_gene_id(val)

    def merge(self,e,new_id=None):
        """ merge self and SGENT e.
        """
        out = copy.deepcopy(self)
        out.start = min(self.start,e.start)
        out.end = max(self.end,e.end)
        out.update(e)
        if new_id:
            out.edit_gene_id(new_id)
        return out

    def as_gtf(self):
        """ For GTF
        """
        return "\n".join(map(lambda x:x.as_gtf(),self.values()))

    def get_represent_tx(self):
        """ Get representative transcript of one gene, usally the longest one.
            Useful when doing gene function annotation."""
        max_len = 0 
        for tx in self.values():
            if len(tx) > max_len:
                rep_tx = tx
                max_len = len(tx)
        return rep_tx

class ChrDict(dySGENT):

    def __init__(self,es=[],perm_op=True,engine="dypylib"):
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
            self.add_e(Gene(Transcript([exon])))

    def add_tx(self,tx):
        if tx.gene_id in self.keys():
            self[tx.gene_id].add_tx(tx)
        else:
            self.add_e(Gene([tx]))

    def add_gene(self,gene):
        self.add_e(tx)

    def add_e(self,gene):
        assert isinstance(gene,Gene)
        self.ID = gene.Chr
        #super(ChrDict,self).add_e(gene)
        self.add_data(gene)
        self.add_len(gene)

    def sort(self):
        self.sort(key=lambda x:x.start)

    def __len__(self):
        return self.get_count()

class GenomeDict(dySGENT):
    """to-do: Remove this object"""

    def add_exon(self,exon):
        if exon.Chr in self.keys():
            self[exon.Chr].add_exon(exon)
        else:
            self.add_e(Chr(Gene(Transcript([exon]))))

    def add_e(self,gene):
        assert isinstance(gene,Chr)
        super(GenomeDict,self).add_e(gene)

    def add_gene(self,gene):
        self.add_e(Chr(gene))

    def add_gene_rec(self,rec):
        gene = Gene(rec=rec)

    def check_strand(self,e):
        pass

    def __len__(self):
        return self.get_count()

class GTF(object):

    """
    A wraper for tkGTF

    """

    def __init__(self, input_obj, ID=None, *argv, **kwargv):
        """
        The constructor

        Since tkGTF can be constructed with both File, (tk)GTF and string object, the input_obj can be any of them, theoretically.
        Because tkGTF is slow, it is time comsuming to reconsctruct a GTF object, while a wrapper is more efficient.

        :param input_obj: tkGTF object.
        "init_children bool: Set to True to call self.__init_children.
        """
        self.wrapped_class = input_obj
        self.ID = ID

    def __getattr__(self,attr):
        orig_attr = self.wrapped_class.__getattribute__(attr)
        # why if callable???
        #if callable(orig_attr):
        #    return orig_attr
        if(self == self.wrapped_class):
            print("self == wrapped")
        return orig_attr

    def __getitem__(self,key):
        return self.get(key)

    def __str__(self):
        """
        Make the output of print looks like:

        - A GTF object containing 14 lines:
        - File connection: tmp/tttt_met_eb.gtf
        
        CI01000000	0,2,2	gene	252526	259879	1000.00	+	.	gene_id	"CIWT.10";
        CI01000000	0,2,2	transcript	252526	259876	1000.00	+	.	transcript_id "CIWT.10.1"; gene_id "CIWT.10";
        CI01000000	0,2,2	exon	252526	252935	1000.00	+	.	transcript_id "CIWT.10.1"; gene_id "CIWT.10";
        CI01000000	0,2,2	exon	253671	253903	1000.00	+	.	transcript_id "CIWT.10.1"; gene_id "CIWT.10";
        CI01000000	0,2,2	exon	256634	256857	1000.00	+	.	transcript_id "CIWT.10.1"; gene_id "CIWT.10";
        CI01000000	0,2,2	exon	257469	257660	1000.00	+	.	transcript_id "CIWT.10.1"; gene_id "CIWT.10";

        Otherthan:
        <bio.seq.base_md.GTFTx object at 0x7fd251e3e358>
        """
        return str(self.wrapped_class)

    @staticmethod
    def get_first_feature(input_obj):
        """
        Get the first tkFeature object of a tkGTF
        using iteration.

        :param input_obj: A tkFeature object
        """
        assert isinstance(input_obj, tkGTF)
        for f in input_obj:
            return f

    @staticmethod
    def get_key_feature_obj(input_obj, key):
        """
        Help to get a primary Feature object.
        For-loop of python is much more faster than get_key_feature of tkGTF.

        :param input_obj: A tkFeature object
        :param key: feature name, like transcript, gene and exon.
        """
        if isinstance(input_obj, tkGTF):
            for rec in input_obj:
                if rec.ft_type == key:
                    return rec
        elif isinstance(input_obj,tkFeature):
            return input_obj
        elif isinstance(input_obj, GTF):
            for rec in input_obj.wrapped_class:
                if rec.ft_type == key:
                    return rec
        else:
            print(type(input_obj))
            raise TypeError

    @staticmethod
    def _get_key_feature_obj(input_obj, key):
        """
        Help to get a primary Feature object.

        :param input_obj: A tkFeature object
        :param key: feature name, like transcript, gene and exon.
        """
        if isinstance(input_obj,tkFeature):
            return input_obj
        mygtf = input_obj.select_by_key('feature',key)
        return GTF.get_first_feature(mygtf)
    
class iGTF(GTF):
    """
    Subclass of GTF, the iterable version of GTF.
    """

    def __init__(self, input_obj, ID=None, init_children=True,as_dict=False,*argv,**kwags):
        super(iGTF,self).__init__(input_obj, ID=ID)
        if init_children:
            self._init_children()
            self._length = len(self._keys)

        # If as_dict, store hash relationship.
        # Else, call value by dynamic search.
        self.is_dict = False
        if as_dict:
            self._init_dict()

    def _init_children(self):
        """
        Crucial for iteration.
        Called when init_children is set to True in __init__
        """
        self._keys = self.get_chroms(nr=True)

    def _init_dict(self):
        """
        Init a static dictionary, which can accelerate the extracting operation.
        """
        self._dict = Ordic()
        for key, val in self.items():
            self._dict[key] = val
        self.is_dict = True

    def __len__(self):
        return self._length

    def __iter__(self):
        if self.is_dict:
            return self._dict.values()
        else:
            return iter(self._keys)

    def items(self):
        if self.is_dict:
            return self._dict.items()
        
        for key in self:
            yield key, self[key]

    def keys(self):
        if self.is_dict:
            return self._dict.keys()
        else:
            return self._keys

    def values(self):
        if self.is_dict:
            return self._dict.values()

        for key in self:
            yield self[key]

class fGTF(GTF):
    """
    Featured GTF.
    A fGTF object always has a representative feature line in GTF(ensembl).

    :Example:
    GTFGene
    GTFTx
    GTFExon
    GTFCds
    """
    _fdic = {'chr':'chrom','chrom':'chrom','seqname':'chrom','seqid':'chrom',
             'source':'src','src':'src',
             'type':'ft_type','feature':'ft_type','ft_type':'ft_type',
             'start':'start',
             'end':'end',
             'score':'score',
             'strand':'strand',
             'phase':'frame','frame':'frame'}

    def __init__(self, input_obj, primary_rec=None, ID=None, *argv, **kwargv):
        """
        :param input_obj: Typically GTF object, sometimes Feature object.
        :param primary_rec: The represent Feature of GTF object, the first line will be used when None is given.
        :param ID: The same with GTF __init__.
        """
        super(fGTF,self).__init__(input_obj,ID, *argv, **kwargv)

        # Set primary_rec to call some attribution in _fdic using "." method.
        # When the input_obj is tkGTF, choose the first record as p_r when
        # the primary_rec isn't provided in the parameters.
        if isinstance(input_obj,tkFeature):
            self.primary_rec = input_obj
        elif isinstance(input_obj,tkGTF):
            self.primary_rec = primary_rec if primary_rec else GTF.get_first_feature(input_obj)
        else:
            raise TypeError("Only tkGTF and tkFeature is supported.")

    def __getattr__(self,attr):
        if attr in fGTF._fdic:
            return getattr(self.primary_rec, fGTF._fdic[attr])
        return self.wrapped_class.__getattribute__(attr)

    def __len__(self):
        return self.end - self.start

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

class eGTF(fGTF):

    """
    Base element GTF.
    A eGTF object always has a representative feature line in GTF(ensembl).

    :Example:
    GTFExon
    GTFCds
    """
    #def __len__(self):
    #    return self.end - self.start
    pass

class GTFCds(eGTF):pass

class GTFExon(eGTF):pass

class GTFTx(fGTF):
    """
    Object bind closely to GTFGene

    :Relationship:

    :Parent class: GTFGene
    :Child class: EXON

    """
    def __init__(self, input_obj, ID=None):
        """
        GTF.get_key_feature_obj function is very slow due to select_by_key.
        Hence, only use get_first_feature instead
        """
        super(GTFTx,self).__init__(input_obj,
                                   primary_rec = GTF.get_key_feature_obj(input_obj, 'transcript'),
                                   #primary_rec = GTF.get_first_feature(input_obj),
                                   ID=ID)

    def __len__(self):
        _length = 0
        for e in self:
            _length += len(e)
        return _length

    def __iter__(self):
        return iter(self.get_exons())

    def get_exons_old(self):
        """
        Select_by_key is extrimly slow, avoid using!
        """
        return map(GTFExon,self.select_by_key('feature','exon'))

    def get_cds_old(self):
        return map(GTFCds,self.select_by_key('feature','cds'))

    def _iter_type(self, key):
        for rec in self.wrapped_class:
            if rec.ft_type == key:
                yield rec

    def get_exons(self):
        for rec in self._iter_type('exon'):
            yield GTFExon(rec)

    def get_cds(self):
        for rec in self._iter_type('cds'):
            yield GTFCds(rec)

    def get_introns(self):pass

class GTFGene(iGTF,fGTF):
    """
    Object bind closely to GTFChr

    :Relationship:

    :Parent class: GTFChr
    :Child class: GTFTx

    """
    def __init__(self, input_obj, ID=None):
        super(GTFGene,self).__init__(input_obj,
                                     primary_rec = GTF.get_key_feature_obj(input_obj, 'gene'),
                                     #primary_rec = GTF.get_first_feature(input_obj, 'gene'),
                                     ID=ID,
                                     as_dict= True)

    def _init_children(self):
        """
        Crucial for iteration.
        Called when init_children is set to True in __init__
        """
        self._keys = self.get_tx_ids(nr=False)

    def get(self,key):
        """
        Core function for bunch of methods, including
        __getitem__, iteration, items(), keys(), values()
        
        Divergent child class which represent certain biology segement
        can be created by simply rewriting this function.
        """
        return GTFTx(self.select_by_key("transcript_id",key),key)

class GTFChr(iGTF):
    """
    Object bind closely to GTFGenome

    :Relationship:

    :Parent class: GTFGenome
    :Child class: GTFGene

    """
    def __init__(self, input_obj, ID=None):
        #super(GTFChr,self).__init__(input_obj,ID=ID,as_dict=True)
        super(GTFChr,self).__init__(input_obj,ID=ID)

    def _init_children(self):
        """
        Crucial for iteration.
        Called when init_children is set to True in __init__
        """
        self._keys = self.get_attr_value_list(key="gene_id")

    def get(self,key):
        """
        Core function for bunch of methods, including
        __getitem__, iteration, items(), keys(), values()
        
        Divergent child class which represent certain biology segement
        can be created by simply rewriting this function.
        """
        return GTFGene(self.select_by_key("gene_id",key),key)

    def get_tx_feature(self):
        return list(map(GTFTx,self.select_by_key('feature','transcript')))

    def get_txs(self):
        ids = self.get_tx_ids(nr=True)
        return list(map(lambda x:GTFTx(self.select_by_key('transcript_id',x)), ids))

class GTFGenome(iGTF):
    """
    Toplevel bio feature handle for GTF
    The constructor.

    :Relationship:

    :Depending class: tkGTF

    """
    def __init__(self,gtf):
        """
        :param gtf: File, tkGTF, or any other object accepted by tkGTF.
        """
        super(GTFGenome,self).__init__(tkGTF(gtf))

    def get(self,key):
        """
        Core function for bunch of methods, including
        __getitem__, iteration, items(), keys(), values()
        
        Divergent child class which represent certain biology segement
        can be created by simply rewriting this function.
        """
        return GTFChr(self.select_by_key("chr",key),key)

    def get_chr_to_x(self,e_type,strand=False):
        """ Get a list of {chr:[transcripts]}
            :param: e_type, element type, can be one of [tx, gene].
            :param: strand, set strand to true to include strand in
                key. The key will be (chr, strand) rather than chr.
        """
        if e_type == "tx":
            _key = "transcript_id"
            myclass = GTFTx
        elif e_type == "gene":
            _key = "gene_id"
            myclass = GTFGene

        my_dict = self.extract_data("chrom,"+_key,
                                    as_dict_of_merged_list=True,
                                    nr=True,
                                    hide_undef=True)
        outdict = Ordic()
        for key, val in my_dict.items():
            for v in val:
                e = myclass(self.select_by_key(_key,v),v)
                #if strand == True:
                outkey = (key,e.strand) if strand else key
                _vals = outdict.setdefault(outkey, [])
                _vals.append(e)
            
        return outdict

    def get_chr_to_tx(self,if_strand=False):
        return self.get_chr_to_x("tx",if_strand)

    def get_chr_to_gene(self,if_strand=False):
        return self.get_chr_to_x("gene",if_strand)

    #def get_chr_to_gene1(self,strand=False):
    #    my_dict = self.extract_data("chrom,gene_id",
    #                                as_dict_of_merged_list=True,
    #                                nr=True,
    #                                hide_undef=True)
    #    outdict = Ordic()
    #    for key, val in my_dict.items():
    #        _vals = []
    #        for v in val:
    #            _vals.append(GTFTx(self.select_by_key('gene_id',v),v))
    #        outdict[key] = _vals
    #        
    #    return outdict

    #def get_chr_to_tx1(self,strand=False):
    #    """ Get a list of {chr:[transcripts]}
    #        :param: strand, set strand to true to include strand in
    #            key. The key will be (chr, strand) rather than chr.
    #    """
    #    my_dict = self.extract_data("chrom,transcript_id",
    #                                as_dict_of_merged_list=True,
    #                                nr=True,
    #                                hide_undef=True)
    #    outdict = Ordic()
    #    for key, val in my_dict.items():
    #        for v in val:
    #            tx = GTFTx(self.select_by_key('transcript_id',v))
    #            if strand==True:
    #                key = (key,tx.strand)
    #            _vals = outdict.setdefault(key, [])
    #            _vals.append(tx)
    #        
    #    return outdict

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
            # Some illegle records.
            if (rec.gene_id is None or
                rec.tx_id is None):
                continue
            if rec.type == "exon":
                self.add_exon(dyExon(rec))
            elif rec.type == "gene":
                self.add_gene_rec(rec)

    # Convert gffdict to {tx_id:tx_obj}
    def get_TxDict(self):
        data = Ordic()
        for mychr in self.values():
            for gene in mychr.values():
                for tx in gene.values():
                    data[tx.ID] = tx
        return data

    def _iter_chr_to_genes(self):
        for _chrid, _chr in self.items():
            _genes = _chr.values()
            yield _chrid, _genes

    def _iter_chr_to_txs(self):
        for _chrid, _genes in self._iter_chr_to_genes():
            _txs = []
            for _gene in _genes:
                _txs.extend(_gene.values())
            yield _chrid, _txs

    def get_chr_to_x(self,e_type,strand=False):
        """ Get a list of {chr:[transcripts/genes]}
            :param: e_type, element type, can be one of [tx, gene].
            :param: strand, set strand to true to include strand in
                key. The key will be (chr, strand) rather than chr.
        """
        if e_type == "tx":
            _iter_func = self._iter_chr_to_txs
        elif e_type == "gene":
            _iter_func = self._iter_chr_to_genes

        data = Ordic()

        # for the saving of time
        if not strand:
            for _chrid, _xs in _iter_func():
                data[_chrid] = _xs
            return data

        for _chrid, _xs in _iter_func():
            for _x in _xs:
                outkey = (_chrid, _x.strand) if strand else _chrid
                _vals = data.setdefault(outkey, [])
                _vals.append(_x)
        return data
            
    def get_chr_to_gene(self,if_strand=False):
        return self.get_chr_to_x("gene",if_strand)

    def get_chr_to_tx(self,if_strand=False):
        return self.get_chr_to_x("tx",if_strand)

    def conv_gffdict2tx_level(self):
        return self.get_TxDict()

class dyGenome(dySGENT):
    
    """ Read gtf file, convert to GenomeDict data.
        Only exon are used.

        Example(Maybe correct):
        genomedict = dyGenome(input_gtf_file)
        # Get a geneDict object for gene named "EXMP"
        dict_of_gene = genomedict.get_GeneDict()
        gene = dict_of_gene.get("EXMP")
        # Get information of the gene
        print(gene.start)
        print(gene.end)
        print(gene.strand)
        # Length of gene (sum of exon length), might be wrong, Fix in next version
        print(len(gene))
        # How many transcript in a gene
        print(gene.count)

        to-do: Rename this object as Genome,
        and keep this object as alias of Genome

    """
    def __init__(self,infile,sep="\t",keeps=set(["exon"]),fm="normal", engine="dypylib"):
        """ fm[normal,gc,none]: format of gtf/gff file.
                normal, in most situation, normal format are suitable.
                gc, grass carp gtf V1.
                none, don't parse geneID and txID.
        """
        super(dyGenome,self).__init__([])
        self.loadfile(infile,sep,keeps,fm)

    def loadfile(self,infile,sep,keeps,fm):
        for rec in Gff(infile,sep=sep,fm=fm):
            # filter type
            #print(type(rec))
            if (keeps != None) and (rec.type not in keeps):
                continue
            # Some illegle records without 
            # gene_id or transcript_id attribution.
            if (rec.gene_id is None or
                rec.tx_id is None):
                continue
            self.add_exon(Exon(rec))

    def __len__(self):
        return self.get_count()

    def add_exon(self,exon):
        if exon.Chr in self.keys():
            self[exon.Chr].add_exon(exon)
        else:
            self.add_e(Chr(Gene(Transcript([exon]))))

    def add_e(self,gene):
        assert isinstance(gene,Chr)
        super(dyGenome,self).add_e(gene)

    def add_gene(self,gene):
        self.add_e(Chr(gene))

    def add_gene_rec(self,rec):
        gene = Gene(rec=rec)

    def check_strand(self,e):
        pass

    def get_GeneDict(self):
        """ Get a dict of gene_id - geneDict pair.
            {gene_id:GeneDict}
        """
        data = Ordic()
        for mychr in self.values():
            for gene in mychr.values():
                data[gene.ID] = gene
        return data

    # Convert gffdict to {tx_id:tx_obj}
    def get_TxDict(self):
        """ Get a dict of tx_id - txDict pair.
            {tx_id:TxDict}
        """
        data = Ordic()
        for gene in self.get_GeneDict().values():
            for tx in gene.values():
                data[tx.ID] = tx
        return data

    def rename_gene_tx(self, gene_id, tx_id = None, get_id = False):
        """ Rename the gene id with gene_id, transcript id with tx_id.
            When tx_id is not provide, name it after gene_id.
            example:
                1.
                gene_id = HBG
                tx_id = None
                xxxx ... transcript_id "HBG.1.1"; gene_id "HBG.1"
                xxxx ... transcript_id "HBG.1.2"; gene_id "HBG.1"
                xxxx ... transcript_id "HBG.2.1"; gene_id "HBG.2"
                2.
                gene_id = HBG
                tx_id = HBT
                xxxx ... transcript_id "HBT.1.1"; gene_id "HBG.1"
                xxxx ... transcript_id "HBT.1.2"; gene_id "HBG.1"
                xxxx ... transcript_id "HBT.2.1"; gene_id "HBG.2"
            Return a tuple (A, B)
            A is a dict, which is same as the return of get_GeneDict().
            B is a list, records old/new id.
            (old id, new id, "g"/"t") is the element of B, where "g" means gene id and "t" means transcript id.
            If get_id set to False, only return A.
        """
        data = Ordic()
        id_mapping = []
        if tx_id == None:
            tx_id = gene_id
        # 1-based results
        gene_count = 1
        geneDict = self.get_GeneDict()
        for ogid, txs in geneDict.items():
            # gene suffix
            gsf = "." + str(gene_count)
            gid = gene_id + gsf
            outgene = data.setdefault(gid,Ordic())
            id_mapping.append((ogid,gid,"g"))
            # 1-based results
            tx_count = 1
            for otid, exons in txs.items():
                # transcript suffix
                tsf = "." + str(tx_count)
                tid = tx_id + gsf + tsf
                outtx = outgene.setdefault(tid,Ordic())
                id_mapping.append((otid,tid,"t"))
                # 0-based exon-count
                exon_count = 0
                for exon in exons.values():
                    exon.edit_gene_id(gid)
                    exon.edit_tx_id(tid)
                    outtx[exon_count] = exon
                    exon_count += 1
                tx_count += 1
            gene_count += 1
        if get_id:
            return data, id_mapping
        else:
            return data

class GtfDict(dyGenome):

    def __init__(self, *args, **kwargs):
        print("Warning, Please use Genome rather than GtfDict.\
 The later one will retired soon.")
        super(GtfDict, self).__init__(*args, **kwargs)
                    
# Get direction of two element
# "u" e1 locate in up (5')
# "d" e1 locate in down (3')
# "o" two element are overlaped.
def get_dir(st1,ed1,st2,ed2):
    def co(b1,a2,f = "u"):
        if b1 >= a2:
            return "o"
        else:
            return f
    if st1 < st2:
        return co(ed1,st2,"u")
    else:
        return co(ed2,st1,"d")

def get_dir_e(e1,e2):
    st1,ed1 = e1.start, e1.end
    st2,ed2 = e2.start, e2.end
    return get_dir(st1,ed1,st2,ed2)

def check_overlap(st1,ed1,st2,ed2):
    if get_dir(st1,ed1,st2,ed2) == "o":
        return True
    else:
        return False

# Abord 19/1/12
def check_overlap1(st1,ed1,st2,ed2):
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
            #print("gene1,gene2 = %s, %s" % (data[i].ID,data[i+1].ID))
            #print("s1,e1,s2,e2 = %d, %d, %d, %d" % (data[i].start,data[i].end,data[i+1].start,data[i+1].end))
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
        
#class Dictsql(object):
#    """A dictionary based on sqlit database"""
#    def __init__(self, db = None, name = ""):
#        self.db = db
#        self.name = name


# The following four class is an example classes
#   that prevent create object from destiny class (
#   Last_of_us).
#   In a nother way, the instance of Last_of_us
#   must be create through the factory class (Test)
class Test_main(type):
    def __call__(self, *args, **kwargs):
        if len(args) == 0 or args[0] != 'factory':
            raise TypeError("Can't instantiate directly")
        else:
            return super(Test_main, self).__call__(*args, **kwargs)

class Last_of_us(metaclass=Test_main):
    def __init__(self, *args, **kwargs):
        print("init Last_of_us")
        print(type(self))

class Uncharted(object):
    def __init__(self):
        print("init Uncharted")
        print(type(self))

class Test(ABC):
    games = {'last_of_us': Last_of_us, 'uncharted': Uncharted}
    def __new__(cls, name, *args, **kwargs):
        if name in cls.games:
            # Add 'factory' as the first parameter to 
            #   indicate this class are created by Factory class.
            myclass = cls.games[name]('factory', *args, **kwargs)
            cls.register(cls.games[name])
            return myclass
        else:
            return PSGame()

class GffutilsGene(GffutilsSGENT):
    """Chromosome object
    """

    def _gffutils__getitem__(self, key):
        if key not in self:
            info = "{} is not contained in {}, please check it.\
                    ".format(key, self.name)
            raise KeyError(info)
        return Transcript(self.db, name=key, engine='gffutils')

    def _gffutils_get_children_id(self):
        children = []
        records = self.db.children(self.name, featuretype = 'transcript')
        for rec in records:
            children.append(rec.id)
        return children

class GffutilsChr(GffutilsSGENT):
    """Chromosome object
    """

    def _gffutils__getitem__(self, key):
        if key not in self:
            info = "{} is not contained in {}, please check it.\
                    ".format(key, self.name)
            raise KeyError(info)
        return Gene(self.db, name=key, engine='gffutils')

    def _gffutils_get_children_id(self):
        keys = []
        order = "select id from features where seqid = '%s' AND featuretype = 'gene'" \
                % self.name
        c = self.db.execute(order)
        for i in set(c.fetchall()):
            keys.append(i['id'])
        return keys

class GffutilsGenome(GffutilsSGENT):
    """Genome object.
    
    The children can be chromosomes, schaffolds, or any other
    features previous defined in the annotation file.

    to-do: Rewrite SGENT.
    """

    def _gffutils__getitem__(self, key):
        if key not in self:
            info = "{} is not contained in {}, please check it.\
                    ".format(key, self.name)
            raise KeyError(info)
        return Chr(self.db, name=key, engine='gffutils')

    def _gffutils_get_children_id(self):
        keys = []
        c = self.db.execute("select DISTINCT seqid from features")
        for i in c.fetchall():
            keys.append(i['seqid'])
        return keys

    def convert_feature(self, feature):
        return self._gffutils_convert_feature(feature)

    def _gffutils_convert_feature(self, feature):
        """Convert gffutil object to dypylib object"""
        if feature.featuretype == 'gene':
            return Gene(self.db, name = feature.id, engine='gffutils')
        elif feature.featuretype == 'transcript':
            return Transcript(self.db, name = feature.id, engine='gffutils')
        elif feature.featuretype == 'exon':
            return Exon(self.db, name = feature.id, engine='gffutils')
        elif feature.featuretype == 'cds':
            return CDS(self.db, name = feature.id, engine='gffutils')
        else:
            return feature

    def search(self, key):
        return self._gffutils_search(key)
    
    def _gffutils_search(self, key):
        """Search by feature id.
        Except for chromosome, all other ids like gene/transcript id
            are stored at database.
        """
        # Return a chromsome
        if key in self:
            return self[key]
        else:
            return self.convert_feature(self.db[key])

class GenomeFeature(ABC):
    """Factory Class for Genomic Features"""
    engines = {}
    def __new__(cls, *args, **kwargs):
        engine = default_param(kwargs, 'engine', 'dypylib')
        if engine in cls.engines:
            myclass = cls.engines[engine](*args, **kwargs)
            myclass.init_from_factory = True
            cls.register(cls.engines[engine])
            return myclass
        else:
            raise TypeError("\"{}\" is not a valid engine.\
                    Engines in {} is supported".format(engine,\
                    ", ".join(cls.engines.keys())))

class Genome(GenomeFeature):
    engines = {'dypylib': dyGenome, 'gffutils': GffutilsGenome}

class Chr(GenomeFeature):
    engines = {'dypylib': ChrDict, 'gffutils': GffutilsChr}

class Gene(GenomeFeature):
    engines = {'dypylib': GeneDict, 'gffutils': GffutilsGene}

class Transcript(GenomeFeature):
    engines = {'dypylib': TxDict, 'gffutils': GffutilsTranscript}

class Exon(GenomeFeature):
    engines = {'dypylib': dyExon, 'gffutils': GffutilsExon}

class CDS(GenomeFeature):
    engines = {'dypylib': dyCDS, 'gffutils': GffutilsCDS}

class Intron(GenomeFeature):
    engines = {'dypylib': dyIntron, 'gffutils': GffutilsIntron}
    
def create_genome_using_gffutils(infile, dbfn = ':memory:'):
    """Create Genome object using gffutils
        1. Create gff_db
        2. Create Genome object with gff_db

        to-do: Write method for creating database, which would help
            to efficient create biology objects.
    """
    try:
        gff_db = gffutils.create_db(infile, dbfn)
    # gff_db has already created and saved in memory
    except sqlite3.OperationalError:
        gff_db = gffutils.interface.FeatureDB(dbfn)
    return create_genome_from_gffutils_db(gff_db)

def create_genome_from_gffutils_db(db):
    return Genome(db, engine = 'gffutils')

def test_genome_object(test_gtf):
    mygenome = create_genome_using_gffutils(test_gtf)
    print("")

    print(mygenome)

    print("Test keys():")
    for i in mygenome.keys():
        print(i)

    print("Test for-loop")
    for i in mygenome:
        print(i)

    print("Test items")
    for k,v in mygenome.items():
        print(k)
        print(v)

def default_param(kwargs, key, default):
    return default if key not in kwargs else kwargs[key]

if __name__ == '__main__':

    import sys 
    test_gtf = "../../tests/data/test.gtf"
    test_genome_object(test_gtf)
    ###x = BioMapping(db = '111')
    ###for i in x:
    ### print(i)
    #gff_db = gffutils.create_db(test_gtf,':memory:')
    #gff_db = gffutils.interface.FeatureDB(':memory:')
    #c = gff_db.cursor()
    #c = gff_db.execute("select seqid from features")

    #for i in c.fetchall():
    #    print(i['seqid'])
    #print(gff_db.featchall())

    #print(gffutils.constants.SCHEMA)
