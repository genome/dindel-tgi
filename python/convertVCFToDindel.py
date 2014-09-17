#!/usr/bin/env python
import os, sys, glob, gzip, math
from optparse import OptionParser
utilPath = 'utils/' # specify absolute path of directory with the Dindel python utils/ subdirectory
if not os.path.exists(utilPath):
    sys.stderr.write("Please specify correct path for the Dindel 'utils/' subdirectory by editing third line of this script.\n")
    sys.exit(1)
sys.path.append(utilPath)
import FileUtils
import Fasta, AnalyzeSequence, Variant
import VCFFile

def convert(inputVCFFile = '', outputVariantFile = '', parameters = {}):

    
    fo = open(outputVariantFile, 'w')
    
    fa = Fasta.Fasta(fname = parameters['refFile'])

    vcffiles = inputVCFFile.split(',')

    for vcffile in vcffiles:
        vcf = VCFFile.VCFFile(fname = vcffile, mode = 'r')

        while True:
            dat = vcf.readline()

            if dat == {}:
                break

            pos = int(dat['POS'])
            chr = dat['CHROM']
            ref = dat['REF']

            rseq = ''.join(fa.get(chr, pos, len(ref)))
            if rseq != ref:
                sys.stderr.write("REFSEQ inconsistency\n")

            if float(dat['QUAL'])>=parameters['minQual']:
                altseq = dat['ALT'].split(',')

                for alt in altseq:
                    if alt != "<DEL>" and len(alt)!=len(ref):
                        var = Variant.Variant4(ref = ref, alt = alt)
                        if var.type == "ins" or var.type == "del":
                            fo.write("%s %d %s\n" % (chr, pos+var.offset-1, var.str))
                        

    fo.close()

    vcf.close()
    

def main(argv):
    parser = OptionParser()

    parser.add_option("-i","--inputFile", dest = "inputFile", help = "VCF File")
    parser.add_option("-o","--outputFile", dest = "outputFile", help = "output file with Dindel-style variants calls") 
    parser.add_option("-r","--refFile", dest = "refFile", help = "reference sequence _indexed_ Fasta file")
    parser.add_option("--minQual", dest = "minQual", help = "minimum mapping quality", default = 1)

    (options, args) = parser.parse_args() 

    if options.inputFile == None:
        sys.stderr.write("Please specify --inputFile\n")
        sys.exit(1)
    if options.outputFile == None:
        sys.stderr.write("Please specify --outputFile\n")
        sys.exit(1)
    if options.refFile == None:
        sys.stderr.write("Please specify --refFile\n")
        sys.exit(1)


    parameters = {
            'refFile':options.refFile,
            'minQual':float(options.minQual)
        }

    convert(inputVCFFile = options.inputFile,  outputVariantFile = options.outputFile, parameters = parameters)




if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except:
        sys.stderr.write("An error occurred!\n")
        raise
        sys.exit(1)


