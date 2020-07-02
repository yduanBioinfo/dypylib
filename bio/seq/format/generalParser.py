#!/usr/bin/env python3

import os, re
import pandas as pd

import sys
"""
    APIs:
    read_Tbout_raw
    read_Tbout
    read_Rfam_family
    read_spe_exp (deprecated)
    read_spe_exps (deprecated)
    read_ioe_file (deprecated): ioe file of suppa: hq_lq_lnc_b_all_strict.event.ioe
"""

def _tbout_get_header(src, comments='#'):
    """
        Estimate the number of field of tboutput from description lines.
        Assume the following format.
        example of header for hmmscan tboutput.
        #                                                                        --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
        # target name        accession  query name                    accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        #------------------- ----------          -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        :param:src: file handle.
    """
    sep = re.compile(r'-\s')
    theHeader = 2
    theNumber = 3
    counts = 0
    out_header = []
    for line in src:
        line = line.strip()
        if line[0] == comments:
            line = line.lstrip(comments)
            counts += 1
        if counts == theHeader:
            header = line
        if counts == theNumber:
            # Get index according to the third line.
            idx = [0]
            idx.extend([i.end() for i in sep.finditer(line)])
            for i in range(len(idx)-1):
                out_header.append(header[idx[i]:idx[i+1]].strip())
            return out_header

def _convert_tbout_to_tab_delimiter(src,comments='#',raw_length=18):
    """ The number of elements in each line is constant equal to raw_length.
        No space was found in elements except the last one.
    """
    for line in src:
        line = line.strip()
        if line[0] == comments:
            continue
        line_array = line.split()
        last_element = " ".join(line_array[raw_length-1:])
        line_array2 = line_array[:raw_length-1]
        line_array2.append(last_element)
        yield line_array2

def _pick_best_one(df):
    """ While multiple alignment can be found for one certain sequence,
        a representative record should be chosen as the main feature of the source sequence.
        Criterion is based on E-value.
    """
    output = []
    for k,gp in df.groupby("query_name"):
        # Get min Evalue rows,
        # which could be more than one where two row both have minimal E-value and maximum score).
        gp = gp[gp.E_value == gp.E_value.min()]
        gp = gp[gp.score == gp.score.max()]
        if gp.shape[0] > 1:
            gp = gp.iloc[[1],:]
        output.append(gp)
    return pd.concat(output)

def read_Tbout_raw(source):
    """ Parser for tbout file without filtering step.
    """
    # for cmmscan
    #header =["target_name", "target_accession", "query_name", "query_accession", "mdl", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc", "pass", "gc", "bias", "score", "E_value", "inc", "description_of_target"]
    header = _tbout_get_header(open(source))
    df = pd.DataFrame.from_records(_convert_tbout_to_tab_delimiter(open(source), raw_length=len(header)))
    df.columns=header
    return df

def read_Tbout(source,threshold=0.0001,only_one=False):
    """ Parser for tbout file of cmscan program.
        The output file will be read as DataFrame object of pandas,
        following by filter out untrustable alignment accoding to the threshold on E-value.
        When only_one has set to True, Only one record with the smallest E-value and the
        largest score is retained per query.
    """
    df = read_Tbout_raw(source)
    df = df.astype({'E-value':'float32'})
    # Filter out untrustable alingments with a relative strict threshold.
    if threshold != 0:
        df = df[df['E-value']<threshold]
    if only_one:
        df = _pick_best_one(df)
    return df

def read_Rfam_family(source="/home/yduan/database/Rfam/Rfam_14.1/family.txt"):
    df = pd.read_csv(source,sep='\t',encoding = "ISO-8859-1",header=None)
    names = list(df.columns)
    names[0] = "accession"
    names[1] = "short_name"
    names[3] = "whole_name"
    names[9] = "description"
    names[18] = "type"
    df.columns = names
    return df

def read_spe_exp(source,drop_names = False):
    df = pd.read_csv(source,sep='\t',header=None,names=['time_type','tx_count','gene_count','tx_names'])
    # Both versions of split work.
    df['time'], df['type'] = df['time_type'].str.split(',').str
    #df = pd.concat((df,df['time_type'].str.split(',',expand=True).rename(columns={0:'time',1:'type'})),1)

    df = df.drop('time_type',1)
    # Reorder columns.
    df = df[['time', 'type', 'tx_count', 'gene_count', 'tx_names']]
    # drop tx_names
    if drop_names:
        df = df.drop('tx_names',1)
    return df

def _spe_exp_parse_name(name):
    """
    hq_lq_lnc_spe_f99.txt --> (hq_lq_lnc, f99)
    hq_lnc1_spe_1.txt --> (hq_lnc1, 1)
    """
    pt = re.compile(r'^(\S+)_spe_([^.]+)\.txt$')
    name = os.path.basename(name)
    return pt.match(name).groups()

def _read_spe_exp(source):
    count = read_spe_exp(source, drop_names = True)
    count['merge'], count['exp_thre'] = _spe_exp_parse_name(source)
    return count

def read_spe_exps(sources):
    """ Read specific expression files into one dataframe.
        And add two columns that are merge(hq_lq_lnc/hq_lnc1/hq_lnc2) and exp_thre(1/5/10...).
        Names are dropped by default.
    """
    data = []
    for s in sources:
        data.append(_read_spe_exp(s))
    return pd.concat(data)

def _get_stru(name):
    """ Get type of event from event_id.
    """
    id_desc = name.split(";")
    mytype = id_desc[1].split(":")[0]
    return mytype

def read_ioe_file(source):
    """Read ioe file and parse type of event."""

    df = pd.read_csv(source, sep='\t')
    df['event_type'] = df['event_id'].apply(_get_stru)
    return df
    
