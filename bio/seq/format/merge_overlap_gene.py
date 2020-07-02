#!/usr/bin/env python3

import sys, argparse
from bio.seq.base import GtfDict, check_overlap_e, check_overlap

def seperate_by_strand(genes):
    """ Seperate_genes by strand
    """
    pstv = [] # strand +
    ngtv = [] # strand -
    us = [] # strand ., unstrand
    for gene in genes:
        if gene.strand == "+":
            pstv.append(gene)
        elif gene.strand == "-":
            ngtv.append(gene)
        elif gene.strand == ".":
            us.append(gene)
        else:
            raise KeyError("undefind strand")
    return pstv, ngtv, us

def merge_overlap1(genes,pre="IHBMG",id_count=0):
    genes = sorted(genes, key = lambda x: x.start)
    curr_gene = genes[0]
    ifmerge = False
    for i in range(1,len(genes)):
        if check_overlap_e(curr_gene,genes[i]):
            curr_gene = curr_gene.merge(genes[i],pre+str(id_count))
            ifmerge = True
        # not overlap, so output.
        else:
            yield curr_gene
            curr_gene = genes[i]
            if ifmerge:
                id_count += 1
                ifmerge = False
    yield curr_gene
    id_count += 1
    #return id_count

def merge_genes(genes,new_id=None):
    curr_gene = genes[0]
    for i in range(1,len(genes)):
        curr_gene = curr_gene.merge(genes[i])
    if new_id:
        curr_gene.edit_gene_id(new_id)
    return curr_gene

def detect_fusion(data):
    """ input [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (6, 7)]
        output [(0,(1,2,3,4))]
               [(index of fusion gene,index of member gene)]
    """ 

    def _get_over_dic(data):
        # val of over_dic: [index of overlap genes]
        count = {}
        def _add_count(a,b):
            count.setdefault(a,[]).append(b)
       
        for pair in data:
            _add_count(pair[0],pair[1])
            _add_count(pair[1],pair[0])
        return count
    
    def _get_max_count(over_dic):
        mymax = (0,[])
        for key, val in over_dic.items():
            #print("key")
            #print(key)
            #print("val")
            #print(val)
            if len(val) > len(mymax[1]):
               mymax = (key, val)
        return mymax

    def _del_e(data,name):
        outdata = []
        for pair in data:
            if name not in pair:
                outdata.append(pair)
        return outdata

    fusion = [] # Save the results.
    while 1:
        over_dic = _get_over_dic(data)
        mymax = _get_max_count(over_dic)
        if len(mymax[1]) < 2:
            return fusion
        else:
            fusion.append(mymax)
            data = _del_e(data,mymax[0])

def check_overlap_frac(st1,ed1,st2,ed2,thre1=0.6,thre2=0.6):
    if not check_overlap(st1,ed1,st2,ed2):
        return "N"
    overlap = min(ed1,ed2) - max(st1,st2) + 1
    frac1 = overlap/(ed1-st1+1)
    frac2 = overlap/(ed2-st2+1)
    if frac1 > thre1 and frac2 > thre2:
        return "OC"
    else:
        return "OU"

def merge_genes_in_block(genes):
    """ OU: two gene overlap, but not satisfied the threshold
        OC: two gene overlap and satisfied the threshold
        N: two gene not overlap
        results:
            gene cluster, fusion gene --- members.
            gene cluster:[[gene_1,gene_2],[gene_3],[gene_4,gene_5,gene6]]
            fusion:[(gene_1,(m_gene1,m_gene2)),(gene_2,(m_gene3,m_gene5,m_gene6))]

    """

    def _get_fusion_name(fusion):
        return (genes[fusion[0]].ID,tuple(map(lambda x:genes[x].ID,fusion[1])))

    def get_fusion_name(fusions):
        res = []
        for fu in fusions:
            res.append(_get_fusion_name(fu))
        return res

    def _update_st_ed(st1,ed1,st2,ed2):
        st = min(st1,st2)
        ed = max(ed1,ed2)
        return st, ed

    choosed = set() # The index of genes that have been collapsed.
    l = [] # lable record OU
    c = [] # collapsed genes
    for i in range(len(genes)):
        if i in choosed:
            continue
        choosed.add(i)
        curr_st = genes[i].start
        curr_ed = genes[i].end
        curr_collapsed = [i] # The index of genes OC.
        for j in range(len(genes)):
            if j in choosed:
                continue
            res = check_overlap_frac(curr_st,curr_ed,genes[j].start,genes[j].end)
            if res == "OU":
                l.append((i,j))
            elif res == "OC":
                curr_st, curr_ed = _update_st_ed(curr_st, curr_ed, genes[j].start, genes[j].end)
                curr_collapsed.append(j)
                choosed.add(j)
            elif res == "N":
                #c.append(curr_collapsed)
                break
        c.append(curr_collapsed)
    return map(lambda x:map(lambda y: genes[y],x),c), get_fusion_name(detect_fusion(l))

def iter_overlap_block(genes):
    """ Seperate whole annotation to different blocks by overlapping.
    """
    # Emplty block
    if not genes:
        return []
    genes = sorted(genes, key = lambda x: x.start)
    block = []
    curr_gene = genes[0]
    block.append(curr_gene)
    for i in range(1,len(genes)):
        if check_overlap_e(curr_gene,genes[i]):
            curr_gene = curr_gene.merge(genes[i])
            block.append(genes[i])
        # not overlap, so output.
        else:
            yield block
            curr_gene = genes[i]
            block = []
            block.append(curr_gene)
    yield block

def merge_and_write(genes,outfile,pre="IHBMG",id_count=0,ref_file=None):
    """ genes demo: [[gene_a,gene_b],[gene_c],[gene_d,gene_e,gene_f]]
        The first and the last list in demo list will be merged, then named after "pre" and id_count.
        output:
        xxx, pre_0(gene_a_b), xxx
        xxx, gene_c, xxx
        xxx, pre_1(gene_d_e_f), xxx
    """
    def _write_ref_file(genes,c_id,outfile,sep="\t"):
        if not outfile:
            return
        for gene in genes:
            outfile.write(gene.ID+sep+c_id+"\n")

    for sub_genes in genes:
        sub_genes = list(sub_genes) # convert map object to list.
        # collapse
        if len(sub_genes) > 1:
            c_id = pre+str(id_count)
            outfile.write(merge_genes(sub_genes,c_id).as_gtf()+"\n")
            id_count += 1
            _write_ref_file(sub_genes,c_id,ref_file)
        # output directly
        else:
            outfile.write(sub_genes[0].as_gtf()+"\n")
    return id_count

#def purish_fusion(ref_table,raw_fusion,out_fusion):
#    """ Read ref_table in to dict.
#        Repace fusion names with the dict and filter out those 
#    """
#    pass

def output_fusion(genes,outfile,sep="\t"):
    if not outfile:
        return
    for fu_mem in genes:
        for m in fu_mem[1]:
            outfile.write(fu_mem[0]+sep+m+"\n")

def merge_gtf(gtf,outfile,pre,ref_file,fusion_file):
    id_count = 0
    for genes in gtf.values():
        _genes = seperate_by_strand(genes.values())
        for strand in _genes:
            for block in iter_overlap_block(strand):
                gene_cluster, fusion_pairs = merge_genes_in_block(block)
                #### debug####
                #for i in gene_cluster:
                #    for j in i:
                #        print(j.ID)
                #### ########
                # merge and output
                id_count = merge_and_write(gene_cluster,outfile,id_count=id_count,ref_file=ref_file)
                # output fusion genes
                output_fusion(fusion_pairs,fusion_file)

def main(argv):
    parser = argparse.ArgumentParser(description="Merge overlaped gene.")
    parser.add_argument('gtf',nargs='?',help="input gtf file.",default = sys.stdin,type=GtfDict)
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-g','--gene-id',nargs='?',default="IHBMG",help = "Prefix of gene id")
    parser.add_argument('-r','--refer-table',nargs='?',default=None,help ="Output tabular file to record old/new id",type=argparse.FileType('w'))
    parser.add_argument('-f','--fusion-gene',nargs='?',default=None,help ="Output fushion gene id and its target.",type=argparse.FileType('w'))
    args = parser.parse_args(argv[1:])

    merge_gtf(args.gtf,args.outfile,args.gene_id,args.refer_table,args.fusion_gene)

if __name__ == '__main__':
    main(sys.argv)
