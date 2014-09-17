import os, sys, gzip
import string

class FileWithHeader:
    def __init__(self, fname = '', mode='r', headerLabels = [], NAchar = 'NA', joinChar = ' ', headerLines = []):
        if mode == 'r':
            self.mode = 'r'
        elif mode == 'w':
            self.mode = 'w'
        else:
            raise NameError('Unrecognized mode')

        self.fname=fname 
        self.joinChar = joinChar

        self.NAchar = NAchar 

        if self.mode == 'w':
            if os.path.exists(fname):
                sys.stderr.write('\033[93m'+'Overwriting '+fname+'\n\033[0m')
            
            if os.path.splitext(fname)[-1] == '.gz':
                self.f = gzip.open(fname, 'w')
            else:
                self.f = open(fname,'w')

            self.headerLabels = headerLabels
       
            if string.lower(fname[-3:]) == "vcf" or string.lower(fname[-6:]) == "vcf.gz":
                for line in headerLines:
                    if line[:2] == '##':
                        self.f.write("%s\n" % line)

                self.f.write('#'+joinChar.join(headerLabels)+'\n')
            else:
                self.f.write(joinChar.join(headerLabels)+'\n')

        elif self.mode == 'r':
            if os.path.splitext(fname)[-1] == '.gz':
                self.f = gzip.open(fname, 'r')
            else:
                self.f = open(fname,'r')
            if string.lower(fname[-3:]) == "vcf" or string.lower(fname[-6:]) == "vcf.gz":
                while True:
                    line = self.f.readline().rstrip("\n").rstrip()
                    headerLines.append(line)
                    if len(line)>=2:
                        if line[0] == '#' and line[1] != '#':
                            headerLine = line[1:]
                            break

            else:
                headerLine = self.f.readline().rstrip("\n").rstrip()
            
            self.headerLabels = headerLine.split(self.joinChar)
        
            # store potential VCF headers

            self.headerLines = headerLines
        
        # setup header

        self.lab_to_col = {}
        for i in range(len(self.headerLabels)):
            self.lab_to_col[self.headerLabels[i]] = i

        self.emptyline = [self.NAchar] * len(self.headerLabels)
        self.nlabnotfound = 0
        self.numLabels = len(self.headerLabels)
    def ltc(self):
        return self.lab_to_col;

    def hl(self):
        return self.headerLabels;

    def writelineFromLine(self, line = ''):
        if self.mode == 'w':
            self.f.write(line+'\n')
        else:
            raise NameError('File not opened for writing')

    def writelineFromCol(self, col = []):
        if self.mode != 'w':
            raise NameError('File not opened for writing')

        if len(col) != len(self.headerLabels):
            raise NameError('Number of columns does not match!')
        dat = [str(v) for v in col] 
        self.f.write(self.joinChar.join(dat)+'\n')

    def writeline(self, data = {}):
        if self.mode != 'w':
            raise NameError('File not opened for writing')

        out = self.emptyline[:]
        for key,value in data.iteritems():
            if not key in self.lab_to_col:
                self.nlabnotfound += 1
            else:
                out[self.lab_to_col[key]]=str(value)
        
        self.f.write(self.joinChar.join(out)+'\n')

    def readline(self):
        if self.mode != 'r':
            raise NameError('File not opened for reading!')
        rline = self.f.readline()
        line = rline.rstrip("\n").rstrip(' ').split(self.joinChar)
        if line == ['']:
            return dict()

        if len(line)!=self.numLabels:
            #print line
            #print self.headerLabels
            sys.stderr.write("Line in file %s does not have the correct number of labels\n" % (self.fname))
            1/0
            return None
        else:
            res = {}
            for idx in range(self.numLabels):
                res[ self.headerLabels[idx] ] = line[idx]
            return res
    def readlineList(self):
        if self.mode != 'r':
            raise NameError('File not opened for reading!')

        line = self.f.readline().rstrip("\n").rstrip().split(self.joinChar)
        if line == ['']:
            return []

        if len(line)!=self.numLabels:
            #print line
            #print self.headerLabels
            sys.stderr.write("Line in file %s does not have the correct number of labels\n" % (self.fname))
            return None
        else:
            return line

    def readlines(self):
        self.f.seek(0)
        line = self.f.readline()
        while True:
            dat = self.readline()
            if dat == {}:
                break
            yield dat



    def close(self):
        self.f.close()
        if self.mode == 'w' and self.nlabnotfound>0:
            sys.stderr.write('WARNING: Could not find label: '+str(self.nlabnotfound)+' times in '+self.fname+'\n')





def sortFile(fname = '', foutname='', col = 1, splitChar = ''):
    if os.path.splitext( fname )[-1] == '.gz':
        isZipped = 1
    else:
        isZipped = 0

    data = {}

    

    if isZipped == 1:
        f1 = gzip.open(fname,'r')
        fout = gzip.open(foutname,'w')
    else:
        f1 = open(fname,'r')
        fout = open(foutname,'w')


    try:
        sortcol = int(col)-1
    except ValueError:
        line = f1.readline()
        if splitChar == '':
            headerLabels = line.rstrip("\n").split()
            sortcol = headerLabels.index(col)
            fout.write(line)
        else:
            headerLabels = line.rstrip("\n").split(splitChar)
            sortcol = headerLabels.index(col)
            fout.write(line)

    for line in f1.readlines():
        if splitChar == '':
            coldata = line.rstrip("\n").split()
        else:
            coldata = line.rstrip("\n").split(splitChar)

        try:    
            val = float(coldata[sortcol])
            if data.has_key(val):
                data[val].append(line)
            else:
                data[val]=[line]
        except ValueError:
            print 'Columns contains non-numeric data'
            sys.exit(2)
        
    keys = data.keys()
    keys.sort()
    
    for key in keys:
        for line in data[key]:
            fout.write(line)

    fout.close()
    f1.close()

    # os.system("mv %s %s" % (fname+'sorted', fname))


def getFirstLine(fname = ''):
    
        isZipped1 = 0 
        if os.path.splitext( fname )[-1] == '.gz':
            isZipped1 = 1

        if isZipped1 == 1:
            f1 = gzip.open(fname,'r')
        else:
            f1 = open(fname,'r')
       
        fl = f1.readline()
        f1.close()
        return fl

