# Test filter_with_exp.py

from dypylib.bio.seq.RNA_Seq.filter_with_exp import get_count_df
data = "tests/data/tx.tpm"

def test_filter_name():
    res=["CIWT.8135.1","CIWT.8173.1","CIWT.8547.4"]
    count_df = get_count_df(data,10,2)
    assert res == list(count_df[count_df].index)
