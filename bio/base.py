#!/usr/bin/env python3

import sys
import logging

class BaseFile(object):

    def __init__(self, filename):
        self.filename = filename
        if filename:
            logging.debug("Load file `{0}`".format(filename))

class LineFile(BaseFile, list):
    """
    Generic file parser for line-based files
    """
    def __new__(self,*args,**kwargs):
        return list.__new__(self)

    def __init__(self, filename, comment=None, load=False):
        super(LineFile, self).__init__(filename)
        if load:
            fp = must_open(filename)
            #self.lines = [l.strip() for l in fp if l[0]!=comment]
            for l in fp:
                if l[0]!=comment:
                    self.append(l.strip())
            logging.debug("Load {0} lines from `{1}`.".\
                        format(len(self), filename))

class LineFileIterator(LineFile):
    """
    Iterator for line file.
    """
    def __init__(self, filename, comment=None, has_header=False):

        super(LineFileIterator,self).__init__(filename,comment,load=False)
        self.has_header = has_header
        self.fp = must_open(self.filename)
        self.comment = comment

        if has_header:
            self._readheader()

    def __iter__(self):

        return self.next()

    def _readheader(self):

        for line in self.fp:
            if self.is_comment(line):
                continue
            else:
                self._header = line.rstrip("\n")
                break

    def next(self):

        return self._iterator()

    def _iterator(self):
        
        for line in self.fp:
            line = line.rstrip("\n")
            if self.is_comment(line):
                continue
            else:
                yield line

    def is_comment(self,line):

        #Skip all line, when no comment defined.
        if not self.is_comment:
            return False
        if line[0] == self.comment:
            return True
        else:
            return False

    @property
    def header(self):

        if not self.has_header:
            raise TypeError("While has_header has set to False, can not return file header.\n")
        return self._header

class LineFileSpliter(LineFileIterator):
    """
    Iterator for line file.
    Eachline should be split into list.
    """

    def __init__(self,filename, sep="\t", comment=None, has_header=False):

        self.__sep = sep
        super(LineFileSpliter,self).__init__(filename, comment, has_header)

    def _readheader(self):

        super(LineFileSpliter,self)._readheader()
        self._header = self._header.split(self.__sep)

    def _iterator(self):

        for line in self.fp:
            line = line.rstrip("\n")
            if self.is_comment(line):
                continue
            else:
                yield line.split(self.__sep)

class DictFile(BaseFile, dict):
    """
    Generic file parser for multi-column files, keyed by a particular index.
    If valuepos is a single number, DictFile would be {key,value}.
    If valuepos is a list, DictFile would be {key,[values]}.
    """
    def __init__(self, filename, keypos=0, valuepos=1, delimiter=None,
                       strict=True, keycast=None, cast=None, has_header=False):
        super(DictFile, self).__init__(filename)

        # Load attribute.
        self._fp = must_open(filename)
        self.filename = filename
        self.strict = strict
        self.has_header = has_header
        self.__valuepos = self.__value2list(valuepos)#for muti-value purpose
        #ncols correct it self while reading lines.[Hidden bugs]
        # An error is going to triggered when valuepos is None.
        self.ncols = max(keypos, self.max(self.__valuepos)) + 1
        #if has_header:
        #    self._read_header(delimiter)
        self._read_header(delimiter)
        self._read_body(keypos,delimiter,keycast,cast)

    def _read_header(self,delimiter):
        if self.has_header == False:
            return
        self.header = self._fp.readline().rstrip("\n")
        self.header_lst = self.header.split(delimiter)
        # Only keep correct field header accoding to valuepos.
        self.header_lst = list(map(lambda x: self.header.split(delimiter)[x], self.__valuepos)) if self.__valuepos else self.header_lst

    def _read_body(self,keypos,delimiter,keycast,cast):
        for lineno, row in enumerate(self._fp):
            row = row.rstrip()
            atoms = row.split(delimiter)
            # thiscols: Number of columns of current line.
            thiscols = len(atoms)
            if thiscols < self.ncols:
                self._short_line_handle(row,lineno,self.ncols)    
                continue

            key = atoms[keypos]
            value = list(map(lambda x:atoms[x], self.__valuepos)) if (self.__valuepos is not None) else atoms
            if len(value) == 1:
                value = value[0]
            if keycast:
                key = keycast(key)
            if cast:
                value = cast(value)
            self[key] = value

        assert thiscols, "File empty"
        self.ncols = thiscols
        logging.debug("Imported {0} records from `{1}`.".\
                    format(len(self), self.filename))

    def _short_line_handle(self,row,lineno,ncols):
                
        action = "Aborted" if self.strict else "Skipped"

        msg = "Must contain >= {0} columns.  {1}.\n".format(ncols, action)
        msg += "  --> Line {0}: {1}".format(lineno + 1, row)
        logging.error(msg)
        if self.strict:
            sys.exit(1)

    #Convert single item to list.
    def __value2list(self, value):

        #if not self.parse_value:
        #    return value
        if isinstance(value,list):
            return value
        ### Why return None
        if value == None:
            return value
        return [value]

    #get length of value list[Abord]
    def valueLength(self):

        if self.__valuepos == None:
            return self.ncols
        else:
            return len(self.__valuepos)

    def max(self, values):

        if values == None:
            return values
        return max(values)

class SetFile(BaseFile, set):
    """ Read File into set."""
    def __init__(self, filename, valuepos=0, delimiter=False,
                       keycast=None, cast=None, has_header=False):
        super(BaseFile, self).__init__()

        ## Load attribute.
        self._fp = must_open(filename)
        self.filename = filename
        ## Should move this two lines to BaseFile?
        self.has_header = has_header
        #self.valuepos = valuepos
        #if has_header:
        #    self._read_header(delimiter)
        for line in self._fp:
            la = line.strip()
            # Escape from empty lines.
            if len(la) == 0:
                continue
            if delimiter is False:
                item = la
            else:
                la = la.split(delimiter)
                item = la[valuepos]
            self.add(item)

def must_open(filename, mode="r", checkexists=False, skipcheck=False, \
            oappend=False):
    """
    Accepts filename and returns filehandle.

    Checks on multiple files, stdin/stdout/stderr, .gz or .bz2 file.
    """
    if isinstance(filename, list):
        assert "r" in mode

        if filename[0].endswith(".gz") or filename[0].endswith(".bz2"):
            filename = " ".join(filename)  # allow opening multiple gz/bz2 files
        else:
            import fileinput
            return fileinput.input(filename)

    if filename.startswith("s3://"):
        from jcvi.utils.aws import pull_from_s3
        filename = pull_from_s3(filename)

    if filename in ("-", "stdin"):
        assert "r" in mode
        fp = sys.stdin

    elif filename == "stdout":
        assert "w" in mode
        fp = sys.stdout

    elif filename == "stderr":
        assert "w" in mode
        fp = sys.stderr

    elif filename == "tmp" and mode == "w":
        from tempfile import NamedTemporaryFile
        fp = NamedTemporaryFile(delete=False)

    elif filename.endswith(".gz"):
        if 'r' in mode:
            cmd = "zcat {0}".format(filename)
            fp = popen(cmd, debug=False)
        elif 'w' in mode:
            import gzip
            fp = gzip.open(filename, mode)

    elif filename.endswith(".bz2"):
        if 'r' in mode:
            cmd = "bzcat {0}".format(filename)
            fp = popen(cmd, debug=False)
        elif 'w' in mode:
            import bz2
            fp = bz2.BZ2File(filename, mode)
    else:
        if checkexists:
            assert mode == "w"
            overwrite = (not op.exists(filename)) if skipcheck \
                        else check_exists(filename, oappend)
            if overwrite:
                if oappend:
                    fp = open(filename, "a")
                else:
                    fp = open(filename, "w")
            else:
                logging.debug("File `{0}` already exists. Skipped."\
                        .format(filename))
                return None
        else:
            fp = open(filename, mode)

    return fp

if __name__ == '__main__':

    print("hhahah")
