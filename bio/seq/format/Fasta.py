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
NM_	mRNA	Protein-coding transcripts (usually curated)
NR_	RNA	Non-protein-coding transcripts
XM_c	mRNA	Predicted model protein-coding transcript
XR_c	RNA	Predicted model non-protein-coding transcript
https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
RefSeq header:
Case 1
>XM_010970720.1 PREDICTED: Camelus bactrianus transmembrane protein 11 (TMEM11), transcript variant X1, mRNA
Case 2
>XM_007508807.1 Bathycoccus prasinos PREDICTED: zinc finger protein 347-like (Bathy15g00730), partial mRNA
Case 3
>NM_001134752.1 Mus musculus small integral membrane protein 17 (Smim17), mRNA

    XM_010970720.1: accession id
    PREDICTED: status
    Camelus ... X1: Name
    mRNA: sequence type
'''
class EsbSeq(Seq):

    def __init__(self,*args):
    # args: string, name
        if args[0].__class__ == Seq:
            self.parse_name(args[0].name)
        else:
            super(EsbSeq,self).__init__(*args)
            self.parse_name(args[1])

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

    def __iter__(self):

        self.fastaIter=self.iterFasta1()
        return self

    def iterFasta1(self):

        for seq in self.iterFasta():
            yield EsbSeq(seq)

class RefSeq(Seq):
    """
    Attributes:
        ID
        db
        species
        seq_name
    """

    def __init__(self,*args):
    # args: string, name
        if args[0].__class__ == Seq:
            self.parse_name(args[0].name)
        else:
            super(RefSeq,self).__init__(*args)
            self.parse_name(args[1])
        self.ID = self.ID.strip()

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
        myID = tmparray[0].split()
        tmparray2 = tmparray[1].split()
        # case 2
        if len(myID)>1:
            self.ID = myID[0]
            self.species = " ".join(myID[1:])
            if len(myID) > 3:
                sys.stderr.write("Warning: %s is not particular case 2" % tmparray[0])
        # case 1
        else:
            self.ID = tmparray[0]
            self.species = " ".join(tmparray2[:2])
        self.seq_name = " ".join(tmparray2[2:])

    def _normal_name(self,name):
        tmparray = name.split()
        self.ID = tmparray[0]
        self.species = " ".join(tmparray[1:3])
        self.seq_name = " ".join(tmparray[3:])

class RefFa(Fasta):

    def __iter__(self):

        self.fastaIter=self.iterFasta1()
        return self

    def iterFasta1(self):

        for seq in self.iterFasta():
            yield RefSeq(seq)
