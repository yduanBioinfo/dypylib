#!/usr/bin/env python

from ..base import Seq, Fasta


'''
Fasta format for Ensembl database.
Ensembl header:
>ENST00000347977 ncrna:miRNA chromosome:NCBI35:1:217347790:217347874:-1 gene:ENSG00000195671 gene_biotype:ncRNA transcript_biotype:ncRNA
   ^             ^     ^     ^                                          ^                    ^                           ^ 
   ID            |     |  LOCATION                            GENE: gene stable ID       GENE: gene biotype           TRANSCRIPT: transcript biotype   
                 |   STATUS
              SEQTYPE
Note: biotype seems to be ncrna not ncrna:miRNA.

Fasta format for refseq database.
RefSeq header:
>XM_010970720.1 PREDICTED: Camelus bactrianus transmembrane protein 11 (TMEM11), transcript variant X1, mRNA
>NM_001134752.1 Mus musculus small integral membrane protein 17 (Smim17), mRNA

    XM_010970720.1: accession id
    PREDICTED: status
    Camelus ... X1: Name
    mRNA: sequence type
'''
class EsbSeq(Seq):

    def __init__(self,string,name):
        super(EsbSeq,self).__init__(string,name)
        self.parse_name(name)

    def parse_name(self,name):
        def get_value(string):
            tmp = string.split(":")
            return tmp[1]
        tmparray = name.strip().lstrip(">").split()
        self.ID = tmparray[0]
        self.seqtype = tmparray[1]
        self.location = tmparray[2]
        self.gene_id = get_value(tmparray[3])
        self.gene_type = get_value(tmparray[4])
        self.tx_type = get_value(tmparray[5])
        
class EsbFa(Fasta):

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
                yield EsbSeq(seq,name)
                seq = []
                name = line.lstrip(">")
            else:
                seq.append(line)
        seq = "".join(seq)
        yield EsbSeq(seq,name)

class RefSeq(Seq):
    """
    Attributes:
        ID
        db
        species
        seq_name
    """

    def __init__(self,*args):
    # def __init__(self,string,name):
    # args: string, name
        if args[0] == Seq:
            print "aa"
        else:
            super(RefSeq,self).__init__(*args)
            self.parse_name(args[1])

    def parse_name(self,name):
        tmparray = name.strip().lstrip(">").split("PREDICTED:")
        if len(tmparray) > 2:
            raise ValueError("%s has more than one \":\"" % ":".join(tmparray))
        elif len(tmparray) == 2:
            self._predicted_name(tmparray)
        else:
            self._normal_name(tmparray[0])
        self.db = self.ID.split("_")[0]

    def _predicted_name(self,tmparray):
        self.ID = tmparray[0]
        tmparray2 = tmparray[1].split()
        self.species = " ".join(tmparray2[:2])
        self.seq_name = " ".join(tmparray2[2:])

    def _normal_name(self,name):
        tmparray = name.split()
        self.ID = tmparray[0]
        self.species = " ".join(tmparray[1:3])
        self.seq_name = " ".join(tmparray[3:])

class RefFa(Fasta):

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
                yield RefSeq(seq,name)
                seq = []
                name = line.lstrip(">")
            else:
                seq.append(line)
        seq = "".join(seq)
        yield RefSeq(seq,name)
