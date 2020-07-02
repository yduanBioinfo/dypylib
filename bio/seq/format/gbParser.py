#!/usr/bin/env python

'''
This Script is used for parsing GBFF.

_author:You-Duan
_version:0.0.3

'''
import re, sys

def nonepop(mylist):#del last item when its none('')

    if mylist[-1]:
        return mylist
    mylist.pop()
    return mylist
    
def unquota(mystring):

    if mystring[0] == "\"":
        mystring = mystring[1:]
    if mystring[-1] == "\"":
        mystring = mystring[:-1]
    return mystring
    
class Replace(object):

    '''read and save replace file
    
       data:#the data structure is too complex!!!
       {feature1:
       {(subkey1,subkey2,...):(value1,value2,value3,...),
        (subkey1,subkey2,...):(v1,v2,v3,...),...}
        feature2:{},featrue3:{},...}
    '''
    '''can only support recognize and replace in the same feature!!!
    '''
    def __init__(self,replace_file):
    
        self.data = {}
        self.key_names = []#names determine if the feature is finded
        self.target_names = []#names to be changed
        self.load_replace_file(replace_file)
        
    def load_replace_file(self,file,sep="\t"):

        '''replace_file:
           subkey1;subkey2;subkey3;    replace1;replace2;replace3;    feature
           ..;..;..;..;    ..;..;..;    feature2
           ..;..;..;..;    ..;..;..;    feature
           ......
        '''
        firstline = file.readline().split(sep)
        self.key_names = firstline[0].split(";")
        self.target_names = firstline[1].split(";")
        nonepop(self.key_names)
        nonepop(self.target_names)
        for eachline in file:
            tmp = eachline.strip().split(sep)
            area1 = tuple(nonepop(tmp[0].split(";")))
            area2 = tuple(nonepop(tmp[1].split(";")))
            try:
                feature = tmp[2]
            except:
                pass
            self.data.setdefault(feature,{})
            self.data[feature][area1] = area2
            
class Gb_Feature(object):

    '''read and save genebank feature
       
       Public:
       locPattern = location pattern
       eg. "CDS             4071..4523"
       desPattern = description pattern
       eg. "/product="tRNA-Glu""
       deslocPattern = description including location pattern
       eg. "/anticodon=(pos:6702..6704,aa:Gly,seq:tcc)"
       deslocPattern2
       eg. "/transl_except=(pos:6962,aa:TERM)"
       namePattern = ORGANISM name
       eg. "  ORGANISM  Tetrancistrum nebulosi"
       locRdec = description related to location; type:list.
       
       Private:
       fea = feature name; type:string.
       locst = start location; type:int.
       loced = end location; type:int.
       qua = qualifiers; type:dic.
       qua = [key1:[value1,value2,...],key2:[value1,value2,...],...]
       quakeys = qua keys.
    '''
    
    namePattern = re.compile(r'\s{2}ORGANISM\s*(.*)')
    # issue:
    # 1.
    #     tRNA            complement(1064049..1064125)
    # 2.
    #      CDS             complement(join(1310321..1311190,1311192..1311311,
    #                 1311315..1311491))

    locPattern = re.compile(r'\s{5}(\w+)\s+.*')
    locPattern0 = re.compile(r'(\w+)\s+(\d+)\.\.(\d+)\s?')
    # back-up property.
    locPattern0000 = re.compile(r'\s{5}(\w+)\s+(\d+)\.\.(\d+)\s?')
    desPattern = re.compile(r'\s{10,}/(\w+)=(.*)\n')
    desPattern0 = re.compile(r'/(\w+)=(.*)')
    deslocPattern = re.compile(r'(\D+)(\d+)\.\.(\d+)(\D+)')
    deslocPattern2 = re.compile(r'(\D+pos:)(\d+)(\D+)')
    locRdec = []
    
    def __init__(self,file):
    
        self.fea = ""
        self.locst = 1
        self.loced = 1
        self.qua = {}
        self.quakeys = []
        self.loadfile(file)
        
    def get(self,key):
    
        return self.qua.get(key)

    def init_header(self,line):
        tmp = Gb_Feature.locPattern0.match(line)
        self.fea = line.split()[0]
        if tmp:#first line
            self.locst = int(tmp.group(2))
            self.loced = int(tmp.group(3))
        else:
            # not the typical start and end.
            # Don't record.
            pass

    def init_qualifier(self,line):
        tmp = Gb_Feature.desPattern0.match(line)
        if tmp.group(1) not in self.quakeys:
            self.quakeys.append(tmp.group(1))
        #self.qua[tmp.group(1)] = tmp.group(2)
        self.qua.setdefault(tmp.group(1),[]).append(tmp.group(2))
        if tmp.group(2)[0] == "(" and tmp.group(1) not in \
        Gb_Feature.locRdec:#contain position information
            Gb_Feature.locRdec.append(tmp.group(1))

    def init_block(self,data,flag):
        data = "".join(data)
        if flag:
            self.init_qualifier(data)
        else:
            self.init_header(data)

    def loadfile(self,file):
        data = []
        header_loaded = False
        for each in file:
            tmp = Gb_Feature.desPattern.match(each)
            if tmp:#qualifier line
                self.init_block(data, header_loaded)
                data = []
                # The first line is always the location line, therefor,
                # whenever inited, the location must have loaded.
                header_loaded = True
            data.append(each.strip())
        self.init_block(data, header_loaded) 
        
    def change_qua(self,qua_names,replace_dic,rep_names):
    #change qualifiers according to replace_dic
    #qua_names:list or tuple. determine which match situation it is,
    #corresponding to area1 of firstline of replace file.
    #rep_names:list or tuple. determine which items to be replaced,
    #corresponding to area2 of firstline of replace file.
    
        if self.fea not in replace_dic:#none of this feature's bussiness
            return
        my_rep_dic = replace_dic.get(self.fea)
        #my_rep_dic = {
        #(subkey1,subkey2,...):(value1,value2,value3),
        #(subkey1,subkey2,...):(v1,v2,v3,...),
        #...,...,...}
        test_keys = []
        if isinstance(qua_names,str):
            qua_names = [qua_names]
        for each in qua_names:
            try:
                test_keys.append(self.qua.get(each)[-1])#needs fix!!!
            except:
                pass
        tmp = self.mapkeys(test_keys,my_rep_dic)
        if tmp:#Successful matches
            values = my_rep_dic.get(tmp)
            for i in range(len(values)):
                if rep_names[i] in self.qua:#if key not exist,do not creat it.
                    self.qua[rep_names[i]][-1] = "\""+unquota(values[i])+"\""
            
    def mapkeys(self,testkey,mykeys):
    #if testkey in mykeys(to our rules):
    #return the key in mykeys.
    #testkey = (subkey1,subkey2,subkey3) can't be key
    #mykeys = [testkey1,testkey2,testkey3,...]
    #our rules:none type in mykeys matchs every key


        def judge(a,b):#judge subkey
        
            if not b:#b == ''
                return True
            if a == b:
                return True
            return False
        for each in mykeys:
            if not (isinstance(each,list) or isinstance(each,tuple)):
                print("[Error:] please contact the author to solve it.")
                sys.exit()
            for i in range(len(each)):
                if not judge(unquota(testkey[i]),each[i]):
                    break
                if i == len(each)-1:#every testkey[i] matchs each[i]
                    return each
        return None
        
class Genebank_parser(object):
    
    '''read and save genebank file.
       data:
       hddat = head data (never mind); type:list; eachline -> eachelement.
       ftlin = features line; type:string.
       mddat = mid data (features); type:list.
       tldat = tail data (sequence); type:string.
    '''
    
    def __init__(self,file):
    
        self.hddat = []
        self.ftlin = ""
        self.mddat = []
        self.tldat = ""
        self.loadfile(file)
        self.length = len(self.tldat)
        
    def __len__(self):
    
        return self.length
    
    def loadfile(self,file):

        def loadpart(part,file):#decide which part to be loaded
            if part == "A":
                self.loadhdpart(file)
            if part == "B":
                self.loadmdpart(file)
            if part == "C":
                self.loadtlpart(file)
                
        tmpfile = []
        partchanged = False
        currPart = "A"
        for each in file:
            if each[0:8] == "FEATURES":
                nextPart = "B"
                partchanged = True
            if each[0:6] == "ORIGIN":
                nextPart = "C"
                partchanged = True
            if partchanged:
                loadpart(currPart,tmpfile)
                tmpfile = [each]
                partchanged = False
                currPart = nextPart
                continue
            tmpfile.append(each)
        loadpart(currPart,tmpfile)
        
    def loadhdpart(self,file):
        
        for each in file:
            self.hddat.append(each.rstrip())
            
    def loadmdpart(self,file):
        
        self.ftlin = file[0].rstrip()
        currqual = [file[1]]
        for each in file[2:]:
            if Gb_Feature.locPattern.match(each):#location part
                self.mddat.append(Gb_Feature(currqual))
                currqual = [each]
            currqual.append(each)
        self.mddat.append(Gb_Feature(currqual))
            
    def loadtlpart(self,file):
    
        _data = []
        for each in file[1:]:
            tmp = each.strip().split()
            # end of file
            if tmp[0] == "//":
                self.tldat = "".join(_data)
                return
            # The first element of tmp is not sequence.
            # tmp: ['1649941', 'ggcgttaagg', 'cctgcaaaca', 
            # 'ccacgcaggc', 'ttgcggtgcc', 'tgggcatcgg', 
            # 'ctgcagcttc']
            _data.extend(tmp[1:])
                
    def sortmddat(self):
    
        self.mddat = self.mysort(self.mddat,True,False)
        #first sort according to start location
        same_st_indx = currst = 0
        #same_st_indx:start index of the same locations
        #currst = current start location
        currst = self.mddat[0].locst
        for i in range(len(self.mddat)):#seconde sort according to end location
            #apply it only when start locations are same
            if self.mddat[i].locst != currst:#comes to the next start location
                if same_st_indx != i-1:#start locations not unique
                    self.mddat[same_st_indx:i] = \
                    self.mysort(self.mddat[same_st_indx:i],False,True)
                same_st_indx = i
                
    def mysort(self,data,isstart=True,descending=False):
    
        length = len(data)
        if isstart:
            whichloc = ".locst"
        else:
            whichloc = ".loced"
        if descending:
            whichcom = "<"
        else:
            whichcom = ">"
        for i in range(length):
            for j in range(length-i-1):
                if eval("data[j]"+whichloc+whichcom+"data[j+1]"+whichloc):
                    data[j],data[j+1] = data[j+1],data[j]
        return data
        
    def set_site(self,old_site,to_site):
    
        self.set_start(int(old_site)-int(to_site)+1)
        
    def set_start_tl(self,old):
    
        #set old site to 1
        #change the tl part
        self.tldat = self.tldat[old-1:]+self.tldat[:old-1]
        
    def set_start(self,old):
    
        #set old site to 1
        def myfun(loc):#ensure 1 < loc <= length
            if 0 < loc <= self.length:
                return loc
            if loc > self.length:
                return myfun(loc - self.length)
            return (self.length + loc)
        
        old = int(old)
        for each in self.mddat:#each = one Gb_Feature;some times rename it
            if each.fea == "source":
                continue
            each.locst, each.loced = \
            myfun(each.locst-old+1), myfun(each.loced-old+1)
            for i in each.quakeys:#set qualifiers
                if i in Gb_Feature.locRdec:#needs to be repaired!!!
                    tmp = Gb_Feature.deslocPattern.match(each.qua[i][-1])
                    if tmp:
                        tmp = tmp.groups()
                        each.qua[i][-1] = tmp[0]+str(myfun(int(tmp[1])-old+1))\
                        +".."+str(myfun(int(tmp[2])-old+1))+tmp[-1]
                        continue
                    tmp = Gb_Feature.deslocPattern2.match(each.qua[i][-1])
                    if tmp:
                        tmp = tmp.groups()
                        each.qua[i][-1] = tmp[0]+str(myfun(int(tmp[1])-old+1))\
                        +tmp[-1]
                        continue
        self.set_start_tl(old)
        
    def change_qua(self,file):
    
        chgEx = Replace(file)#changeExample
        for each in self.mddat:
            each.change_qua(chgEx.key_names,chgEx.data,chgEx.target_names)
            
    def get_features(self,feature_names=[],filt_feas=[]):
    #get all featrues which name is in feature_names 
    #not in filt_feas and return in a list
    #if feature_names == []:all featrues is needed
    #if filt_feas == []:do not filter
    
        myfeatures = []
        if isinstance(feature_names,str):
            feature_names = [feature_names]
        if isinstance(filt_feas,str):
            filt_feas = [filt_feas]
        tmp = ''
        for each in self.mddat:
            if not feature_names:#all features take into account
                if each.fea not in filt_feas:
                    myfeatures.append(each)
                continue
            if each.fea in feature_names and each.fea not in filt_feas:
                myfeatures.append(each)
        return myfeatures
        
    def get_qualifiers(self,qua_names,fea_names=[],filt_feas=[]):
    
        myquas = []
        if isinstance(qua_names,str):
            qua_names = [qua_names]
        for each in self.get_features(fea_names,filt_feas):#every feature
            for qua_name in qua_names:#return qua_name is first finded
                tmp = each.get(qua_name)
                if tmp:
                    myquas.extend(tmp)
                    break
        return myquas
            
    def writeGBFF(self,file):
    
        #print hddat
        for each in self.hddat:
            file.write(each+"\n")
        #print mddat
        file.write(self.ftlin+"\n")
        for each in self.mddat:
            file.write("     %-16s%d..%d"%(each.fea,each.locst,each.loced)+"\n")
            for eachkey in each.quakeys:
                for i in each.qua[eachkey]:
                    text = "/%s=%s"%(eachkey,i)
                    file.write(("%21s"%"")+("\n%21s"%"").join(text[j*58:(j+1)*58] \
                    for j in range(len(text)//58+1))+"\n")
        #print tldat
        file.write("ORIGIN\n")
        for i in range(0,len(self.tldat),60):
            tmp = self.tldat[i:i+60]
            file.write("%9d "%(i+1))
            file.write(" ".join(tmp[j*10:(j+1)*10] for j in \
            range(min((len(tmp)//10+1),6)))+"\n")
        file.write("//\n")

if __name__ == '__main__':

    import sys, argparse, textwrap
    
    parser = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        prog="gbParser",\
        description="Change Ori_site to New_site",\
        epilog=textwrap.dedent(\
        '''example :
            1.python gbParser.py myGBFF.gb -s 50 -o newGBFF.gb [change site]
            2.python gbParser.py -s 50 -n 3 -o newGBFF.gb myGBFF.gb 
            3.python gbParser.py myGBFF.gb [nothing changed]
            4.python gbParser.py myGBFF.gb -r replace.txt [only replace features]
            5.python gbParser.py myGBFF.gb -g [output a file used to build a tree]
        '''))
    parser.add_argument('infile',help="GBFF file",type=argparse.FileType('r'))
    parser.add_argument('-s','--ori_site',help="the site to be changed",type=int)
    parser.add_argument('-n','--new_site',help="the new site",default=1,type=int)
    parser.add_argument('-r','--rep_file',help="replace file",type=argparse.FileType('r'))
    parser.add_argument('-o','--out_file',help="output gb file [default:stdout]",nargs='?',\
    const=sys.stdout,type=argparse.FileType('w'))
    parser.add_argument('-t','--tree',help="output tree building file [default:tree.txt]",\
    nargs='?',const="tree.txt",type=argparse.FileType('w'))
    myargs = parser.parse_args(sys.argv[1:])
    
    def write_tree(mygb,outfile):
    
        myqua_names = ["gene","product"]
        myfilter = ["CDS"]
        for each in mygb.hddat:#write name 
            tmp = Gb_Feature.namePattern.match(each.rstrip())
            if tmp:
                outfile.write(">"+"_".join(tmp.group(1).split())+"\n")
                break
        outfile.write(" ".join(map(unquota,\
        mygb.get_qualifiers(myqua_names,filt_feas=myfilter))))
    
    def main():
        mygb = Genebank_parser(myargs.infile)
        if myargs.ori_site:
            if myargs.ori_site < 1 or myargs.ori_site > len(mygb):
                print("[Error]:the ori_site you offered is wrong(too big or small),\
try another site")
                sys.exit()
            mygb.set_site(myargs.ori_site,myargs.new_site)
            mygb.sortmddat()
        if myargs.rep_file:
            mygb.change_qua(myargs.rep_file)
        if myargs.tree:
            write_tree(mygb,myargs.tree)
        if myargs.out_file:
            mygb.writeGBFF(myargs.out_file)
    main()
    
