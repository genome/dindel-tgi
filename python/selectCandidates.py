#!/usr/bin/env python
import os, sys, glob, gzip, math
from optparse import OptionParser


def emptyBuffer(fo = 0, variants = {}, parameters = {}, chrom = ''):
    for p in sorted(variants.keys()):
        for v in variants[p].keys():
            if variants[p][v]>=parameters['minCount'] and variants[p][v]<=parameters['maxCount']:
                fo.write("%s %d %s # %d\n" % (chrom, p, v, variants[p][v]))



def selectCandidates(inputFile = '', outputFile = '', parameters = {}):
    
    
    fi = open(inputFile, 'r')
    fo = open(outputFile, 'w')

    visitedChromosomes = {}

    variants = {}
    currChr = -1
    currPos = -1
    while True:
        line = fi.readline()
        if line == '':
            break

        dat = line.split()
        pos = int (dat[1])
        chr = dat[0]
       

        if chr != currChr:
            if visitedChromosomes.has_key(chr):
                raise NameError("Chromosome already processed. Please sort with respect to chromosome first and then position!")
            currChr = chr
            visitedChromosomes[chr]=1
            sys.stderr.write("Analyzing chromosome %s\n" % chr)
            currPos = -1

        if pos < currPos:
            raise NameError("Variants not sorted with respect to position (second column)!")

        # for this line, get variants and frequencies
        y = dat.index('#')
        vars = dat[2:y]
        freqs = dat[y+1:]

        if pos>currPos:
            if currPos != -1:
                emptyBuffer(fo, variants, parameters, chrom = currChr)
            currPos = pos
            variants = { pos: {} }

        if pos == currPos:
            for idx, v in enumerate(vars):
                if not variants[pos].has_key(v):
                    variants[pos][v]=0
                variants[pos][v] += int(freqs[idx])
    
    emptyBuffer(fo, variants, parameters, chrom = currChr)
    
    fo.close()
    fi.close()


def main(argv):
    parser = OptionParser()

    parser.add_option("-i","--inputFile", dest = "inputFile", help = "VCF File")
    parser.add_option("-o","--outputFile", dest = "outputFile", help = "output file with Dindel-style variants calls") 
    #parser.add_option("-r","--refFile", dest = "refFile", help = "reference sequence _indexed_ Fasta file")
    parser.add_option("--minCount", dest = "minCount", help = "minimum count for realigned indel", default = 2)
    parser.add_option("--maxCount", dest = "maxCount", help = "maximum count for realigned indel", default = 10000000)

    (options, args) = parser.parse_args() 

    if options.inputFile == None:
        sys.stderr.write("Please specify --inputFile\n")
        sys.exit(1)
    if options.outputFile == None:
        sys.stderr.write("Please specify --outputFile\n")
        sys.exit(1)
    #if options.refFile == None:
    #    sys.stderr.write("Please specify --refFile\n")
    #    sys.exit(1)


    parameters = {
#            'refFile':options.refFile,
            'minCount':int(options.minCount),
            'maxCount':int(options.maxCount)
        }


    if not os.path.exists(options.inputFile):
        raise NameError("File %s does not exist!" % options.inputFile)
    selectCandidates(inputFile = options.inputFile,  outputFile = options.outputFile, parameters = parameters)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except:
        sys.stderr.write("An error occurred!\n")
        raise
        sys.exit(1)


