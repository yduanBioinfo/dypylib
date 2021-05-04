# Test filter_with_exp.py

from dypylib.bio.seq.RNA_Seq.filter_with_exp import get_count_df
from dypylib.bio.seq.RNA_Seq.filter_with_exp import main as f_main
data = "tests/data/tx.tpm"
GTF_file = "tests/data/test.gtf"

def test_filter_name():
    res=["CIWT.8135.1","CIWT.8173.1","CIWT.8547.4"]
    count_df = get_count_df(data,10,2)
    assert res == list(count_df[count_df].index)

def test_filter_with_keep(capsys):
    ref_output="CIWT.8135.1\nCIWT.8173.1\nCIWT.8540.1\nCIWT.8541.4\nCIWT.8547.4\n"
    f_main("../../bio/seq/RNA_Seq/filter_with_exp.py {} -t 5 -c 2 --write-keep".format(data).split())
    out,err = capsys.readouterr()
    assert out == ref_output

def test_filter_with_GTF(capsys):
    ref_output = ('CI01000023	2,0	exon	139260	139654	.	+	.	transcript_id "CIWT.8135.1"; gene_id "CIWT.8135";\n'
    'CI01000023	2,0	exon	158543	158897	.	+	.	transcript_id "CIWT.8135.1"; gene_id "CIWT.8135";\n'
    'CI01000023	2,0	exon	160924	166249	.	+	.	transcript_id "CIWT.8135.1"; gene_id "CIWT.8135";\n'
    'CI01000023	2,0	exon	166473	167511	.	+	.	transcript_id "CIWT.8135.1"; gene_id "CIWT.8135";\n'
    'CI01000023	0,2	exon	439056	439411	1000.00	-	.	transcript_id "CIWT.8173.1"; gene_id "CIWT.8173";\n'
    'CI01000023	0,2	exon	440755	441146	1000.00	-	.	transcript_id "CIWT.8173.1"; gene_id "CIWT.8173";\n'
    'CI01000023	0,2	exon	441775	441862	1000.00	-	.	transcript_id "CIWT.8173.1"; gene_id "CIWT.8173";\n'
    'CI01000023	0,2	exon	443836	443942	1000.00	-	.	transcript_id "CIWT.8173.1"; gene_id "CIWT.8173";\n'
    'CI01000023	0,2	exon	444021	444075	1000.00	-	.	transcript_id "CIWT.8173.1"; gene_id "CIWT.8173";\n'
    'CI01000023	0,2	exon	444158	444432	1000.00	-	.	transcript_id "CIWT.8173.1"; gene_id "CIWT.8173";\n'
    'CI01000025	0,2	exon	125196	125310	1000.00	-	.	transcript_id "CIWT.8547.4"; gene_id "CIWT.8547";\n'
    'CI01000025	0,2	exon	125521	128773	1000.00	-	.	transcript_id "CIWT.8547.4"; gene_id "CIWT.8547";\n').format()
    f_main('../../bio/seq/RNA_Seq/filter_with_exp.py {} -t 10 -c 2 -g {}'.format(data,GTF_file).split())
    out,err = capsys.readouterr()
    assert ref_output == out
