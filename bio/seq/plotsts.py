#!/usr/bin/env python

import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from matplotlib import rc

class Dist(object):

    '''
    object of distribution data which is convenient for ploting histogram graph.
    attrs {width:xxx,...}
    bins [1,20,44,...]
    '''
    def __init__(self,bins,attrs={}):
        assert (isinstance(bins,list) or isinstance(bins,tuple))
        self.__bins = tuple(map(int,bins))
        self.load_attr(attrs)

    def load_attr(self,attrs):
        if isinstance(attrs,dict):
            self.attrs = attrs
        #assume paired list or tuple.
        else:
            self.attrs = {}
            for key,value in attrs:
                self.add_attr(key,value)

    def add_attr(self,key,value):
        if key in self.attrs.keys():
            print("Warning : %s appears more than once.")
        self.attrs[key] = value

    def get_bins(self):

        return self.__bins

    def get(self,key):

        return self.attrs.get(key)

class Cate(Dist):

    '''
    Category data
    '''

    def __init__(self,counts,labels,attrs={}):
        self.__counts = counts
        self.__labels = labels
        self.load_attr(attrs)

    def get_counts(self):
        return self.__counts

    def get_labels(self):
        return self.__labels

def suffix_end(name,end=".pdf"):

    if name.endswith(end):
        return name
    else:
        return name+end

#Get size of a two dimensional table which have enough room to place elements of a certain length.
#eg1 10 -> (4,3)
#eg2 9 -> (3,3)
#eg3 15 -> (4,4)
def two_dimse(length):

    sq = int(length ** 0.5)
    if sq * sq == length:
        return sq, sq
    elif (sq+1) * sq >= length:
        return sq+1, sq
    elif (sq+1) * sq < length:
        return sq+1, sq+1
    else:
        raise ValueError("Programing Error")

#Iter leaves node(with no child element) of an element object.
def iter_leaves(element):

    for i in element:
        #with no child
        if len(i) == 0:
            yield i

#Get all leaves node(with no child element) of an element object.
def get_leaves_v(element):

    #[(key,value),(key,value),...]
    data = []
    for child in iter_leaves(element):
        data.append(get_raw_stat(child))
    return data

#Remove name space prefix of a given name.
def clean_name(name):

    simbol = "}"
    name = name.split(simbol)[1]
    return name

#Get element name and it's text.
def get_raw_stat(element):

    return clean_name(element.tag), element.text

#Convert an element object to a Dist object.
def get_Dist(element):

    attrs  = get_leaves_v(element)
    elmt_bin = element.find("ns:BinCounts",ns)
    bins = [i.text for i in elmt_bin]
    myDist = Dist(bins,attrs)
    return myDist

#Convert an element object to a Category data.
def get_Cate(element):

    attrs  = get_leaves_v(element)
    counts = [int(i.text) for i in element.find("ns:BinCounts",ns)]
    labels = [i.text for i in element.find("ns:BinLabels",ns)]
    mycate = Cate(counts,labels,attrs)
    return mycate

#Compute GC contant from a root element.
def get_GC(root):

    fractions = root.findall("Stats:TotalBaseFractionPerChannel",ns)
    frac_label = map(lambda x:x.get("Channel"),fractions)
    frac_value = map(lambda x:float(get_raw_stat(list(x)[0])[1]),fractions)
    ifrac_G = frac_label.index("G")
    ifrac_C = frac_label.index("C")
    GC_contant = frac_value[ifrac_G]+frac_value[ifrac_C]

    return GC_contant

#write paired data to a tabular in LaTex format
def write_tab_tex(data):

    #introduction
    intro = r"\begin{tabular}{|c|c|}"
    intro_end = r"\end{tabular}"
    hl = r"\hline "
    rhl = r"\\\hline "
    sep = r"&"
    text = r""
    text += (intro+hl)
    for key, value in data:
        text += str(key)+sep+str(value)+rhl
    text += intro_end

    return text

#plot descriptions in a tabular format.
def plot_desc(pp,des):

    rc('text',usetex=True)
    fig = plt.figure(figsize=(10,10))
    mytext = write_tab_tex(des)
    ax = fig.add_axes([0,0,1,1])
    fig.text(0.15,0.55,mytext,fontsize=20,transform=ax.transAxes)
    fig.suptitle("Basic information",fontsize=30)
    pp.savefig(fig)
    rc('text',usetex=False)

#Plot a histogram with a Dist object at ax level.
def plot_dist(mydist,ax,color="r",plot_subtitle=True,short_subtitle=True):

    mybins=mydist.get_bins()
    width=float(mydist.get("BinWidth"))
    minbin=float(mydist.get("MinBinValue"))
    maxbin=float(mydist.get("MaxBinValue"))
    sampleMean=float(mydist.get("SampleMean"))
    sampleStd=float(mydist.get("SampleStd"))
    numbins = int(mydist.get("NumBins"))

    #get mid value of bin x value.
    x_bins0 = np.arange(minbin,minbin+(numbins+0.5)*width,width)
    x_bins = (x_bins0[1:]+x_bins0[:-1])/2
    midbin = mybins[len(mybins)/2]
    plt.hist(x_bins,len(mybins),weights=mybins,color=color)
    if plot_subtitle:
        subt = mydist.get("MetricDescription")
        if short_subtitle:
            subt = subt.split(":")[1]
        ax.set_title(subt)
    ax.text(0.7,0.7,"mean:%.2f\nstd:%.2f"%(sampleMean,sampleStd),transform=ax.transAxes)
    return ax

##ax level.
#Plot a pie with a Cate object at ax level.
def plot_pie(mycate,ax,colors=[],plot_subtitle=False):

    counts = mycate.get_counts()
    labels = mycate.get_labels()
    total = sum(counts)
    plt.pie(counts,labels=labels,autopct=lambda(p):'{:.2f}% ({:.0f})'.format(p,round(p*total/100)))
    if plot_subtitle:
        ax.set_title(mycate.get("MetricDescription"))
    return ax

#Plot a histogram figure with one Dist object.
def plot_sgl_dist(mydist,fig,header):

    fig.suptitle(header,fontsize=20)
    ax = fig.add_subplot(1,1,1)
    plot_dist(mydist,ax,plot_subtitle=False)

#Plot four histogram in one figure with four Dist object.
def plot_qua_dists(mydists,fig,header):

    fig.suptitle(header,fontsize=20)
    handles = []
    labels = []
    for i in range(4):
        ax = fig.add_subplot(2,2,i+1)
        ax = plot_dist(mydists[i],ax,color=colors[i])

##figure level.
#Plot a pie figure with one Cate object.
def plot_cate_dist(mycate,fig,header=""):

    fig.suptitle(header,fontsize=20)
    ax = fig.add_subplot(1,1,1)
    plot_pie(mycate,ax,plot_subtitle=False)

#plot multi-pie plot in one figure with Catas.
def plot_cate_dist_one(mycates,fig,header=""):

    fig.suptitle(header,fontsize=20)
    size = len(mycates)
    x, y = two_dimse(size)
    for i in range(size):
        ax = fig.add_subplot(x,y,i+1)
        plot_pie(mycates[i],ax,plot_subtitle=True)
    
##data level.
#Put all category plots into one figure, set onePlot=True
def plot_cate_figs(pp,data,onePlot=False,suptitle="Types"):

    if onePlot:
        mycates = []
        fig = plt.figure(figsize=(10,10))
        
    for statn,title in data:
        myCate = get_Cate(root.find("Stats:"+statn,ns))
        if onePlot:
            mycates.append(myCate)
        else:
            fig = plt.figure(figsize=(10,10))
            plot_cate_dist(myCate,fig,title)
            pp.savefig(fig)
    if onePlot:
        plot_cate_dist_one(mycates,fig,suptitle)
        pp.savefig(fig)

#plot single figs(one subplot only).
def plot_sgl_figs(pp,data):

    for statn,title in data:
        fig = plt.figure(figsize=(10,10))
        myDist = get_Dist(root.find("Stats:"+statn,ns))
        plot_sgl_dist(myDist,fig,title)
        pp.savefig(fig)

#plot qua figs in which each contains 4 subplots for different base. Data [(stats name,title),()]
def plot_qua_figs(pp,data):

    for statn,title in data:
        fig = plt.figure(figsize=(10,10))
        myDists = map(get_Dist,root.findall("Stats:"+statn,ns))
        plot_qua_dists(myDists,fig,title)
        pp.savefig(fig)

#input parameters
inFile = "/home/yduan/temp/1_G02/m54061_170910_110737.sts.xml"
outfile = "sts_test.pdf"

def plot_sts(infile,outfile):

    #Set Global parameters 
    global ns 
    ns = {'Stats':"http://pacificbiosciences.com/PacBioPipelineStats.xsd",'xsi':"http://www.w3.org/2001/XMLSchema-instance",'ns':"http://pacificbiosciences.com/PacBioBaseDataModel.xsd"}
    plt.style.use('ggplot')
    global colors
    colors =['m','c','y','g','b','r','b']
    #Parse xml-formated infile
    tree = ET.parse(infile)
    global root 
    root = tree.getroot()

    #Get plot data
    #get informations
    text_desc = get_leaves_v(root)
    GCcontant = get_GC(root)
    text_desc.append(("GCcontant",str(GCcontant)))

    sig_data = [("PausinessDist","Pausiness"),("PulseRateDist","Pulse Rate"),("PulseWidthDist","Pulse Width"),("BaseRateDist","BaseRate"),("BaseWidthDist","Base Width"),("BaseIpdDist","Base Ipd"),("LocalBaseRateDist","Local Base Rate"),("NumUnfilteredBasecallsDist","Unfiltered Base Calls"),("ReadLenDist","Read Length"),("InsertReadLenDist","Insert Read Length"),("MedianInsertDist","Median of Insert"),("HqBaseFractionDist","HQ Base Fraction")]
    qua_data = [("BaselineLevelDist","Base Line Level"),("BaselineStdDist","Base Line Std"),("SnrDist","Signal to Noise ratio"),("HqRegionSnrDist","Signal to Noise ratio(HQ region)"),("HqBasPkMidDist","HQ region PkMid"),("DmeAngleEstDist","DME Angle Estimation")]
    cate_data = [("ProdDist","Productivity"),("LoadingDist","Loading"),("ReadTypeDist","ReadType")]

    #ploting
    pp = PdfPages(outfile)
    plot_desc(pp,text_desc)
    plot_qua_figs(pp,qua_data)
    plot_sgl_figs(pp,sig_data)
    plot_cate_figs(pp,cate_data,onePlot=True)
    pp.close()

def main(argv):

    import argparse

    parser = argparse.ArgumentParser(description="Visualization tool for PacBio sts.xml file.")
    parser.add_argument('infile',nargs='?',default=sys.stdin,help="*sts.xml",type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile',nargs='?',help="output file",default="sts_test.pdf")
    args = parser.parse_args(argv[1:])

    plot_sts(args.infile,suffix_end(args.outfile))

if __name__ == '__main__':

    import sys
    main(sys.argv)

