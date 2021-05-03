#!/usr/bin/env python3

from collections import OrderedDict as Ordic
from dypylib.bio.seq.base import Gff
from random import random, lognormvariate

gene_exp=Ordic()
tx_exp=Ordic()
for rec in Gff(open("test.gtf")):
    #for i in range(10):
    tx_id = rec.tx_id
    gene_id = rec.gene_id

    # It's not the first time meet the tx_id
    # (This is a multi-exon transcript)
    if tx_id in tx_exp:
        continue
    exp = [lognormvariate(1,3) for i in range(3)]
    tx_exp[tx_id] = exp
    gene_old_exp = gene_exp.setdefault(gene_id, [0,0,0])
    for i in range(len(gene_old_exp)):
        gene_old_exp[i] += exp[i]

# output
sep = "\t"
samples = ('mm','cc','ee')
# tx
f1 = open("tx.tpm",'w')
# write header
f1.write("Tx_ID"+sep)
f1.write(sep.join(samples)+"\n")
# write body
for key, val in tx_exp.items():
    f1.write(key+sep+sep.join(map(str,val))+"\n")
f1.close()
# gene
f2 = open("genes.tpm",'w')
# write header
f2.write("Gene_ID"+sep)
f2.write(sep.join(samples)+"\n")
# write body
for key, val in gene_exp.items():
    f2.write(key+sep+sep.join(map(str,val))+"\n")
f2.close()
