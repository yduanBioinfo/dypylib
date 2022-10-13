# Test GFF module

import filecmp
from tempfile import NamedTemporaryFile
from bio.seq.Annotation import create_genome_using_gffutils

test_gtf = "tests/data/test.gtf"
mygenome = create_genome_using_gffutils(test_gtf)

def test_load_GFF():
    target_keys = ['CI01000023', 'CI01000025']
    for i in mygenome.keys():
        assert i in target_keys
    for k,v in mygenome.items():
        assert k in target_keys

def test_Genome_getitem():
    from bio.seq.Annotation import Chr
    mychr = mygenome['CI01000023']
    assert isinstance(mychr, Chr)
    assert mychr.name == 'CI01000023'
    for gene in mygenome['CI01000023']:
        print(gene)
    print("childern")
    #for s in mygenome.db.children('CI01000023'):
    #    print(s)

#def test_filter():
#    tmpfile = NamedTemporaryFile()
#    res_file = "tests/data/scripts_filter/res_file1_file2_filter.txt"
#    order = ["filter.py","-n1","-n2","-m","filter","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name]
#    filter_main(order)
#    assert filecmp.cmp(tmpfile.name, res_file)

