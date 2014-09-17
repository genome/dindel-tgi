import os,sys
from optparse import OptionParser
utilPath = '/home/caa/release-dindel/test_dindel/utils' # specify absolute path of directory with the Dindel python utils/ subdirectory
if not os.path.exists(utilPath):
    sys.stderr.write("Please specify correct path for the Dindel 'utils/' subdirectory by editing third line of this script.\n")
    sys.exit(1)
sys.path.append(utilPath)


import FileUtils
import Fasta
import Variant
import VCFFile


def getCalls(callFile = ''):

    vcf = FileUtils.FileWithHeader(fname = callFile, mode = 'r', joinChar = "\t")

    calls  = {}
    numcalls = 0
    while True:
        dat = vcf.readline()
        if dat == {}:
            break

        if dat['FILTER'] == "PASS" or (dat['FILTER'] == "q20" and float(dat['QUAL'])>=10):
            chrom = dat['CHROM']
            pos = int(dat['POS'])
            ref = dat['REF']
            alt = dat['ALT']
            if alt.find(',')!=-1:
                raise NameError("Cannot deal with these entries")
            
            var = Variant.Variant4(ref = ref, alt = alt)
            # note, must be zero-based pos! SEE 
            newpos = pos + var.offset - 1
            newstr = var.str

            if not calls.has_key(chrom):
                calls[chrom] = {}
            if not calls[chrom].has_key(newpos):
                calls[chrom][newpos] = {}
            if calls[chrom][newpos].has_key(newstr):
                raise NameError('Multiple same variants?')
            
            calls[chrom][newpos][newstr] = dat.copy()
            numcalls += 1

    vcf.close()
    print "Number of calls imported:",numcalls
    return calls


def emptyBuffer(index = -1, buffer = {}, calls = {}, outputFileHandle = -1, bamfiles = []):

    num_bams = len(bamfiles)
    glfs = buffer[index]

    # get first var string

    dat = glfs[0] # first entry


    if dat['nref_all'] == 'NA':
        del buffer[index]
        return "na-error"

    varstring = "%s %s %s" % (dat['tid'], dat['realigned_position'], dat['nref_all'])

    # see if called
    
    try:
        vcfinf = calls[ dat['tid'] ][ int(dat['realigned_position']) ][ dat['nref_all'] ]
    except KeyError:
        del buffer[index]
        return "notcalled"

    # check if length is ok
    if len(glfs) != num_bams:
        sys.stderr.write("Skipping index %s\n" % index)
        del buffer[index] 
        return "skipped"


    # called, get GLFs for all individuals

    output = []

    for dat in glfs:
        tvs = "%s %s %s" % (dat['tid'], dat['realigned_position'], dat['nref_all'])
        if tvs != varstring:
            return "skipped-inconsistent-glf-lines"

        idx = int(dat['indidx'])
        likstr = dat['glf']
        liks = likstr.split(';')
        gen_to_lik = {}
        for lik in liks:
            ld = lik.split(':')
            gen_to_lik[ ld[0] ] = ld[1]

        outstr = "%s %d %s %s %s %s %s\n" % (dat['tid'], int(dat['realigned_position']), dat['nref_all'], gen_to_lik['0/0'], gen_to_lik['0/1'], gen_to_lik['1/1'], bamfiles[idx])

        output.append(outstr)

    # we got to this point, assume everything is ok.

    for line in output:
        outputFileHandle.write(line)

    del buffer[index]
    return "a-ok"

    
def loadGLFFiles(inputGLFFiles = ''):

    glffiles = []
    fg = open(inputGLFFiles, mode = 'r')
    for line in fg.readlines():
        glffiles.append(line.rstrip("\n").split()[0])
    fg.close()

    fp_to_fname = {}

    for glffile in glffiles:
        if not os.path.exists(glffile):
            sys.stderr.write("File %s does not exist\n" % glffile)
            continue
        fg = FileUtils.FileWithHeader(fname = glffile, mode = 'r')
        while True:
            dat = fg.readline()
            if dat['realigned_position'] != 'NA':
                firstpos = int(dat['realigned_position'])
                if fp_to_fname.has_key(firstpos):
                    raise NameError('Huh?')

                fp_to_fname[firstpos] = glffile
                fg.close()
                break

    newglffiles = []
    for pos in sorted(fp_to_fname.keys()):
        print "pos:",pos,"glffile:", fp_to_fname[pos]
        newglffiles.append(fp_to_fname[pos])

    return newglffiles



    
def makeGLF(inputGLFFiles = '', outputFile = '', callFile = '', bamfilesFile = ''):

    # get VCF calls
    
    sys.stdout.write("Reading VCF file\n")
    calls = getCalls(callFile = callFile)
    sys.stdout.write("done\n")
    sys.stdout.flush()

    # read through glf files

    glffiles = loadGLFFiles(inputGLFFiles = inputGLFFiles)

    # get BAMfiles file

    bamfiles = []
    fb = open(bamfilesFile, 'r')
    for line in fb.readlines():
        bamfiles.append(line.rstrip("\n").split()[0])
    fb.close()

    # check each GLF file

    numwritten = 0

    # open output file

    fout = open(outputFile, 'w')

    for glffile in glffiles:
        sys.stdout.write("Checking %s\n" % glffile)
        
        fg = FileUtils.FileWithHeader(fname = glffile, mode = 'r')

        buffer = {}
        curr_index = '-1'
        while True:
            dat = fg.readline()
            if dat == {}:
                break

            newindex = "%s.%s.%s" % (dat['index'], dat['realigned_position'], dat['nref_all'])
            if not buffer.has_key(newindex):
                buffer[newindex] = []
            buffer[newindex].append(dat)

            if newindex != curr_index:
                if curr_index != '-1':
                    result = emptyBuffer(index = curr_index, buffer = buffer, calls = calls, outputFileHandle = fout, bamfiles = bamfiles)
                
                    if result == "a-ok":
                        numwritten += 1
         
                curr_index = newindex

        result = emptyBuffer(index = curr_index, buffer = buffer, calls = calls, outputFileHandle = fout, bamfiles = bamfiles)

        if result == "a-ok":
            numwritten += 1
        
        fg.close()

        print "Number written:", numwritten
        sys.stdout.flush()

    # finish up

    fout.close()

        


                


            


            


def main(argv):
    parser = OptionParser()
    parser.add_option("-g","--inputGLFFiles", dest = "inputGLFFiles", help = "file with all GLFs")
    parser.add_option("-c","--callFile", dest = "callFile", help = "VCF File with calls ( as produced by mergeOutputPooled.py)")
    parser.add_option("-o","--outputFile", dest = "outputFile", help = "output GLF file")
    parser.add_option("-b","--bamFiles", dest = "bamFiles", help = "File with list of BAM files/IDs. NOTE These should be in the same order as in the file given to dindel --bamFiles.")

    (options, args) = parser.parse_args(argv)

    if options.inputGLFFiles == None or options.callFile == None or options.outputFile == None or options.bamFiles == None:
        sys.stderr.write("Please specify all required options.\n")
        sys.exit(1)

    makeGLF(inputGLFFiles = options.inputGLFFiles, outputFile = options.outputFile, callFile = options.callFile, bamfilesFile = options.bamFiles)


if __name__ == "__main__":
    #try:
    main(sys.argv[1:])
    #except NameError:
    #    sys.stderr.write("NameError exception\n")
    #    sys.exit(-1)

