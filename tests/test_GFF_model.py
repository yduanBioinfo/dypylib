# Test GFF module

import filecmp, pytest
from tempfile import NamedTemporaryFile
from bio.seq.Annotation import create_genome_using_gffutils
from bio.seq.Annotation import Chr, Gene, Transcript

test_gtf = "tests/data/test.gtf"
mygenome = create_genome_using_gffutils(test_gtf)
mychr = mygenome['CI01000023']
mygene = mychr['CIWT.8140']

def test_load_GFF():
    target_keys = ['CI01000023', 'CI01000025']
    assert len(mygenome.keys()) == 2
    assert len(mygenome) == 2
    # Test __contains__
    for i in target_keys:
        assert i in mygenome.keys()
    assert sorted(list(set(mygenome.keys()))) == target_keys
    for k,v in mygenome.items():
        assert k in target_keys

def test_Genome_getitem():
    assert isinstance(mychr, Chr)
    assert mychr.name == 'CI01000023'
    mychildren = mychr.get_children_id()
    # Test first and last gene_id is in chromosome
    assert 'CIWT.8168' in mychildren
    assert 'CIWT.8138' in mychildren

def test_Genome_search():
    assert isinstance(mygenome.search('CI01000023'),Chr)
    assert isinstance(mygenome.search('CIWT.8168'),Gene)
    assert isinstance(mygenome.search('CIWT.8138'),Gene)
    assert isinstance(mygenome.search('CIWT.8140.8'),Transcript)
    # Test wrong id
    from  gffutils.exceptions import FeatureNotFoundError
    wrong_id = 'CIWT.814094'
    with pytest.raises(FeatureNotFoundError, match=wrong_id):
        mygenome.search(wrong_id)

def test_Gene():
    assert isinstance(mygene, Gene)
    for gene in mychr.values():
        for tx in gene.values():
            assert isinstance(tx, Transcript)
            assert isinstance(tx.name, str)

def test_Tx():
    assert 'CIWT.8140.asdfui' not in mygenome
    assert 'CIWT.8140.asdfui' not in mygene
    assert 'CIWT.8140' not in mygene
    assert mygene.id == 'CIWT.8140'

    with pytest.raises(KeyError):
        mygenome['CIWT.8140.uioui']

    for tx in mygene.values():
        for i in tx.values():
            assert i.featuretype in ['exon','cds']
