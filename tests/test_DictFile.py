# Test DictFile class

from bio.base import DictFile
data = DictFile("tests/data/sample.dict")

def test_len():
    assert len(data) == 3
    assert data['a'] == 'mm'
    assert data['b'] == 'cc'
    assert data['d'] == 'ee'
