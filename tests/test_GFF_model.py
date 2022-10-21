# Test GFF module

import filecmp, pytest
from tempfile import NamedTemporaryFile
from bio.seq.Annotation import create_genome_using_gffutils
from bio.seq.Annotation import Genome, Chr, Gene, Transcript
from bio.seq.Annotation import GtfDict

test_gtf = "tests/data/test.gtf"
global mygenome, mychr, mygene

def test_load_GFF():
    global mygenome, mychr, mygene
    mygenome = create_genome_using_gffutils(test_gtf)
    mychr = mygenome['CI01000023']
    mygene = mychr['CIWT.8140']

def test_on_genome_object():
    target_keys = ['CI01000023', 'CI01000025']
    assert len(mygenome.keys()) == 2
    assert len(mygenome) == 2
    # Test __contains__
    for i in target_keys:
        assert i in mygenome.keys()
    assert sorted(list(set(mygenome.keys()))) == target_keys
    for k,v in mygenome.items():
        assert k in target_keys

def test_on_chr_object():
    assert isinstance(mychr, Chr)
    assert mychr.name == 'CI01000023'
    mychildren = mychr.get_children_id()
    # Test first and last gene_id is in chromosome
    assert 'CIWT.8168' in mychildren
    assert 'CIWT.8138' in mychildren

def test_on_Gene_object():
    assert isinstance(mygene, Gene)
    for gene in mychr.values():
        for tx in gene.values():
            assert isinstance(tx, Transcript)
            assert isinstance(tx.name, str)

def test_on_Tx_object():
    assert 'CIWT.8140.asdfui' not in mygenome
    assert 'CIWT.8140.asdfui' not in mygene
    assert 'CIWT.8140' not in mygene
    assert mygene.id == 'CIWT.8140'

    # Test wrong key
    with pytest.raises(KeyError):
        mygenome['CIWT.8140.uioui']

    for tx in mygene.values():
        for i in tx.values():
            assert i.featuretype in ['exon','cds']

def test_search_method_of_Genome():
    assert isinstance(mygenome.search('CI01000023'),Chr)
    assert isinstance(mygenome.search('CIWT.8168'),Gene)
    assert isinstance(mygenome.search('CIWT.8138'),Gene)
    assert isinstance(mygenome.search('CIWT.8140.8'),Transcript)
    # Test wrong id
    from  gffutils.exceptions import FeatureNotFoundError
    wrong_id = 'CIWT.814094'
    with pytest.raises(FeatureNotFoundError, match=wrong_id):
        mygenome.search(wrong_id)
    from bio.seq.Annotation import GffutilsGenome
    print(GffutilsGenome(mygenome.db, engine="gffutils"))

def test_GtfDict():
    global mygenome, mychr, mygene
    mygenome = Genome(test_gtf, engine="dypylib")
    mychr = mygenome['CI01000023']
    mygene = mychr['CIWT.8140']
    test_on_genome_object()
    #test_on_chr_object()
    #test_on_Gene_object()

def test_metaClass():
    """Test for the examples using metaclass.
    """
    from bio.seq.Annotation import Test, Last_of_us
    b = Test('last_of_us','asdfw','fawew','werf')
    assert isinstance(b, Test)
    assert not isinstance(Last_of_us, Test)
    with pytest.raises(TypeError):
        l = Last_of_us()

