import os, sys, gzip, string



class VCFFile:
    def writeHeader(self, info, format, filter):
        self.headerLines = []
        self.info = {}
        self.format = {}
        self.filter = {}
        wroteFileFormat = False
        for line in self.generalHeaderLines:
            oline = line.replace('#','')
            if oline.find('fileformat') != -1:
                oline = 'fileformat=VCF4'
                wroteFileFormat = True
            self.f.write("##%s\n" % oline)
            self.headerLines.append(oline)

        if not wroteFileFormat:
            self.f.write("##fileformat=VCF4\n")


        for inf_id in info.keys():
            line = "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (info[inf_id]['ID'], info[inf_id]['Number'],  info[inf_id]['Type'],  info[inf_id]['Description'])
            self.f.write(line+'\n')
            self.headerLines.append(line)
            self.info[ info[inf_id]['ID'] ] = info[inf_id].copy()
        for id in filter.keys():
            line = "##FILTER=<ID=%s,Description=\"%s\">" % (filter[id]['ID'],filter[id]['Description'])
            self.f.write(line+'\n')
            self.headerLines.append(line)
            self.filter[ filter[id]['ID'] ] = filter[id].copy()

        for id in format.keys():
            line = "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (format[id]['ID'], format[id]['Number'],  format[id]['Type'],  format[id]['Description'])
            self.f.write(line+'\n')
            self.headerLines.append(line)
            self.format[ format[id]['ID'] ] = format[id].copy()  

        # write fields
        hstr = '#'
        hstr += '\t'.join(self.headerLabels)
        self.f.write(hstr+'\n')



    def __init__(self,fname = '',  mode = 'r', version = 3.3, headerLabels = [], info = {}, filter = {}, format = {}, generalHeaderLines = []):
        self.requiredLabels = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER', 'INFO']

        self.sepChar = '\t'
        if mode == 'w':
            if version != 4:
                raise NameError('Incorrect version, must be 4')

            if headerLabels == []:
                raise NameError('headerLabels is empty')
           
            if headerLabels[:8] != self.requiredLabels:
                raise NameError('Required VCF labels not present or in correct order')
            if len(headerLabels)>8 and headerLabels[8]!='FORMAT':
                raise NameError('9th column should be FORMAT')
            if len(headerLabels)==9:
                raise NameError('Must have more than 9 columns')

            self.lab_to_col = {}
            for idx, lab in enumerate(headerLabels):
                self.lab_to_col[lab]=idx
            self.headerLabels = headerLabels[:]

            self.generalHeaderLines = generalHeaderLines
            self.f = open(fname, 'w')
            self.mode = 'w'
            self.writeHeader(info, format, filter)
        if mode == 'r':
            self.mode = 'r'
            if headerLabels != []:
                raise NameError("Cannot specify headerLabels when reading VCF file)")

            headerLine = []
            self.headerLabels = []
            self.generalHeaderLines = []
            self.headerLines = []
            self.lab_to_col = {}
            self.mainVCFLabels = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
            self.mainLabels = ['CHROM','POS','ID','REF','ALT','QUAL']


            self.info = {} #{'ID': {'ID':'','Number':'','Type':'','Description':''} }
            self.filter = {} #{'ID':{'Description':''}}
            self.format = {} #{'ID': {'ID':'','Number':'','Type':'','Description':''} }
            
            infoFound = 0
            headerFound = 0
            if os.path.splitext(fname)[-1]=='.gz': 
                self.f = gzip.open(fname, 'r')
            else:
                self.f = open(fname,'r')
            while True:
                line = self.f.readline().rstrip('\n')
                if line[:2] == '##':
                    self.headerLines.append(line)
                    if line.find('fileformat')!=-1:
                        if line.find('VCF')==-1 and line.find('vcf')==-1:
                            raise NameError('Cannot determine VCF version')

                        if line.find('v3')!=-1:
                            self.version = 3
                        elif line.find('v4')!=-1 or line.find('VCF4')!=-1:
                            self.version = 4
                        else:
                            raise NameError('Cannot determine VCF version')


                    if line[:7]== '##INFO=':
                        sstr = line[7:]
                        if string.count(sstr, '"')!=2:
                            raise NameError("Incomplete description: %s" % line)

                        if self.version<4:
                            fields1 = sstr.split('"')
                            fields2 = fields1[0].split(',')
                            
                            id = fields2[0]
                            if fields2[1] == '.':
                                number = '.'
                            else:
                                number = int(fields2[1])

                            type = fields2[2]
                            description = fields1[-2]

                            if self.info.has_key(id):
                                raise NameError("VCF file already has info-id %s" % id)

                            self.info[id]= {'Number':number, 'ID':id, 'Type':type, 'Description':description}
                        elif self.version==4:
                            if sstr[0]!='<' or sstr[-1]!='>':
                                raise NameError("Incorrect description")
                            fields = sstr[1:-1].split(',')
                            _id = fields[0].split('=')
                            if string.upper(_id[0].replace(' ','')) != 'ID':                                
                                raise NameError("Could not find ID in: %s" % line)
                            id = _id[1]

                            _number = fields[1].split('=')
                            if string.upper(_number[0].replace(' ','')) != 'NUMBER':
                                raise NameError("Could not find NUMBER in: %s" % line)
                            number = _number[1]

                            _type = fields[2].split('=')
                            if string.upper(_type[0].replace(' ','')) != 'TYPE':
                                raise NameError("Could not find TYPE in: %s" % line)
                            type = _type[1]

                            _descr = fields[3].split('=')
                            if string.upper(_descr[0].replace(' ','')) != 'DESCRIPTION':
                                raise NameError("Could not find DESCRIPTION in: %s" % line)
                            description = _descr[1]+','.join(fields[4:])

                            self.info[id]= {'Number':number, 'ID':id, 'Type':type, 'Description':description.replace('"','')}



                        else:
                            raise NameError("Not supported yet")
                            
                        #sys.stderr.write("INFO field %s %s\n" % (fields[0], fields[3]))
                        infoFound = 1
                    elif line.find('##FILTER=')==0:
                        # filter lines
                        sstr = line[9:]
                        if string.count(sstr, '"')!=2:
                            raise NameError("Incomplete description: %s" % line)

                        if self.version<4:
                            fields1 = sstr.split('"')
                            fields2 = fields1[0].split(',')
                            
                            id = fields2[0]
                            number = 0
                            type = 'NA'

                            description = fields1[-2]

                            if self.filter.has_key(id):
                                raise NameError("VCF file already has info-id %s" % id)

                            self.filter[id]= {'Number':number, 'ID':id, 'Type':type, 'Description':description.replace('"','')}
                        elif self.version==4:
                            if sstr[0]!='<' or sstr[-1]!='>':
                                raise NameError("Incorrect description")
                            fields = sstr[1:-1].split(',')
                            _id = fields[0].split('=')
                            if string.upper(_id[0].replace(' ','')) != 'ID':                                
                                raise NameError("Could not find ID in: %s" % line)
                            id = _id[1]

                            _descr = fields[1].split('=')
                            if string.upper(_descr[0].replace(' ','')) != 'DESCRIPTION':
                                raise NameError("Could not find DESCRIPTION in: %s" % line)
                            description = _descr[1]+','.join(fields[4:])
                        
                            type = 'NA'
                            number = 'NA'
                            self.filter[id]= {'Number':number, 'ID':id, 'Type':type, 'Description':description.replace('"','')}



                    elif line.find('##FORMAT=')==0:
                        # format lines
                        sstr = line[9:]
                        if string.count(sstr, '"')!=2:
                            raise NameError("Incomplete description: %s" % line)

                        if self.version<4:
                            fields1 = sstr.split('"')
                            fields2 = fields1[0].split(',')
                            
                            id = fields2[0]
                            if fields2[1] == '.':
                                number = 0
                            else:
                                number = int(fields2[1])

                            type = fields2[2]
                            description = fields1[-2]

                            if self.format.has_key(id):
                                raise NameError("VCF file already has info-id %s" % id)

                            self.format[id]= {'Number':number, 'ID':id, 'Type':type, 'Description':description.replace('"','')}
                        elif self.version==4:
                            if sstr[0]!='<' or sstr[-1]!='>':
                                raise NameError("Incorrect description")
                            fields = sstr[1:-1].split(',')
                            _id = fields[0].split('=')
                            if string.upper(_id[0].replace(' ','')) != 'ID':                                
                                raise NameError("Could not find ID in: %s" % line)
                            id = _id[1]

                            _number = fields[1].split('=')
                            if string.upper(_number[0].replace(' ','')) != 'NUMBER':
                                raise NameError("Could not find NUMBER in: %s" % line)
                            number = _number[1]

                            _type = fields[2].split('=')
                            if string.upper(_type[0].replace(' ','')) != 'TYPE':
                                raise NameError("Could not find TYPE in: %s" % line)
                            type = _type[1]

                            _descr = fields[3].split('=')
                            if string.upper(_descr[0].replace(' ','')) != 'DESCRIPTION':
                                raise NameError("Could not find DESCRIPTION in: %s" % line)
                            description = _descr[1]+','.join(fields[4:])

                            self.format[id]= {'Number':number, 'ID':id, 'Type':type, 'Description':description.replace('"','')}






                    else:
                        self.generalHeaderLines.append(line[:])

                elif line[0] == '#' and line[1] == 'C':
                    if not infoFound:
                        sys.stderr.write("WARNING: did not detect INFO fields!\n")
                    headerLine = line[:].replace('#','')
                    self.headerLabels = headerLine.split()
                    for i in range(len(self.headerLabels)):
                        self.lab_to_col[self.headerLabels[i]]=i

                    headerFound = 1
                    self.nonMainLabels = list(set(self.headerLabels)-set(self.mainLabels)-set(['INFO','FILTER']))
                    break


            for lab in self.mainLabels:
                if not lab in self.headerLabels:
                    raise NameError("Could not find column %s in header of VCF file!" % lab)

            self.minLen = max(self.lab_to_col.values())


    def parseline(self, line = ''):
        line=line.rstrip("\n");

        if line == '':
            return {}
        out = {}

        col = line.split(self.sepChar)
        if len(col)<self.minLen:
            sys.stderr.write('Cannot parse this line:\n'+line+'\n')
        else:
            for lab in self.mainLabels:
                out [ lab] = col[ self.lab_to_col[lab] ]

            for lab in self.nonMainLabels:
                out [ lab] = col[ self.lab_to_col[lab] ]


            # output info
            out['INFO'] = {}
            out['_INFO'] = col[ self.lab_to_col['INFO'] ]


            infos = col[ self.lab_to_col['INFO'] ].split(';')
            for info in infos:
                info_pair = info.split('=')
                id = info_pair[0]
                if len(info_pair)>1:
                    vals = info_pair[1].split(',')
                    out['INFO'][id] = vals
    
                    out["INFO_%s" % id] = info_pair[1]
                else:
                    out['INFO'][id] = []
                    out["INFO_%s" % id] = []

            # output filters
            out['FILTER'] = col[ self.lab_to_col['FILTER'] ].split(';')

            for filter in self.filter.keys():
                hs = "FILTER_%s" % filter
                if filter in out['FILTER']:
                    out[hs] = '1'
                else:
                    out[hs] = '0'

            if self.version == 3:
                if out['FILTER'] == ['0'] or out['FILTER'] == ['.']:
                    out['PASSED'] = True
                else:
                    out['PASSED'] = False
            elif self.version == 4:
                if out['FILTER'] == ['PASS'] or out['FILTER'] == ['.']:
                    out['PASSED'] = True
                else:
                    out['PASSED'] = False
            else:
                raise NameError('Internal error')

        
        return out

    def readline(self):
        if self.mode == 'w':
            raise NameError("Cannot read VCF in write mode!")
        return self.parseline(self.f.readline())


    def writeline(self, dat = {}, checkTags = True):
        if self.mode != 'w':
            raise NameError("Can only write in write-mode")
        linedata = []
        for lab in self.headerLabels:
            if not dat.has_key(lab):
                raise NameError("Cannot find label %s in input" % lab)
        
        for lab in ['CHROM','POS','ID', 'REF','ALT','QUAL']:
            if lab == 'ID':
                if dat[lab] == []:
                    linedata.append('.')
                else:
                    linedata.append(','.join(dat[lab]))
            else:
                linedata.append(dat[lab])
        
        # filter
        if dat['FILTER'] == [] or dat['FILTER'] == ['PASS']:
            linedata.append('PASS')
        else:
            if checkTags:
                for fi in dat['FILTER']:
                    if not self.filter.has_key(fi):
                        raise NameError('Undefined filter!')
            linedata.append(';'.join(dat['FILTER']))

        # info
        infos = []
        for inf in dat['INFO'].keys():
            if checkTags:
                if not self.info.has_key(inf):
                    raise NameError('Undefined info tag!')
                if self.info[inf]['Type'] == 'Flag':
                    infos.append(inf)
                else:
                    infos.append("%s=%s" % (inf, ','.join(dat['INFO'][inf])))
            else:
                if dat['INFO'][inf] == []:
                    infos.append(inf)
                else:
                    infos.append("%s=%s" % (inf, ','.join(dat['INFO'][inf])))

        linedata.append(';'.join(infos))

        for x in self.headerLabels[8:]:
            linedata.append(dat[x])
        self.f.write('\t'.join(["%s" % (v) for v in linedata])+'\n')





    def close(self):
        self.f.close()




