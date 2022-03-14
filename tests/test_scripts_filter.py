#!/usr/bin/env python3

import filecmp
from tempfile import NamedTemporaryFile
from scripts.filter import main as filter_main

def test_filter():
    tmpfile = NamedTemporaryFile()
    res_file = "tests/data/scripts_filter/res_file1_file2_filter.txt"
    order = ["filter.py","-n1","-n2","-m","filter","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name]
    filter_main(order)
    assert filecmp.cmp(tmpfile.name, res_file)

def test_diff():
    tmpfile = NamedTemporaryFile()
    res_file = "tests/data/scripts_filter/res_file1_file2_diff.txt"
    order = ["filter.py","-n1","-n2","-m","differ","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name]
    filter_main(order)
    assert filecmp.cmp(tmpfile.name, res_file)

def test_null():
    tmpfile = NamedTemporaryFile()
    res_file = "tests/data/scripts_filter/res_file1_file2_filter.txt"
    order = ["filter.py","-n1","-n2","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name]
    filter_main(order)
    assert filecmp.cmp(tmpfile.name, res_file)

def test_filter_ignorecase():
    tmpfile = NamedTemporaryFile()
    res_file = "tests/data/scripts_filter/res_file1_file2_noCase_filter.txt"
    order = ["filter.py","-n1","-n2","-m","filter","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name,"-i"]
    filter_main(order)
    assert filecmp.cmp(tmpfile.name, res_file)

def test_diff_ignorecase():
    tmpfile = NamedTemporaryFile()
    res_file = "tests/data/scripts_filter/res_file1_file2_noCase_diff.txt"
    order = ["filter.py","-n1","-n2","-m","differ","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name,'--ignorecase']
    filter_main(order)
    assert filecmp.cmp(tmpfile.name, res_file)

def test_null_ignorecase():
    tmpfile = NamedTemporaryFile()
    res_file = "tests/data/scripts_filter/res_file1_file2_noCase_filter.txt"
    order = ["filter.py","-n1","-n2","tests/data/scripts_filter/file1.txt","tests/data/scripts_filter/file2.txt","-o",tmpfile.name,'--ignorecase']
    filter_main(order)
    assert filecmp.cmp(tmpfile.name, res_file)
