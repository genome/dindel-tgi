#!/usr/bin/env python
from optparse import OptionParser
import sys, os, getopt

class Variant:
    def __init__(self, refPos = 0, varString = ''):
        if len(varString)<2:
            raise NameError("Unrecognized variant: "+varString)

        if varString[0] == '-':
            self.type = "del"
            self.seq = varString[1:]
            self.length = len(varString)-1
            self.refStart = refPos
            self.refEnd = refPos+self.length-1
        elif varString[0] == '+':
            self.type = "ins"
            self.seq = varString[1:]
            self.length = len(varString)-1
            self.refStart = refPos
            self.refEnd = refPos-1
        elif len(varString)==4:
            self.type = "snp"
            self.seq = varString[:]
            self.length = 1
            self.refStart = refPos
            self.refEnd = refPos
        else:
             raise NameError("Unrecognized variant: "+varString)

        self.str = varString


def write_output_candidates(newVariants, outputPrefix = '', variantsPerFile = 1000, hapWidth = 60, startIdx = 0, maxVarPerWindow = 16, chromosomes = []):
    idx = startIdx

    numLineWritten = 10000000
    fileOpen = False
    windows = []
    histWindows = {}
    tot = 0
    totVar = 0

    if chromosomes == []:
        chromosomes = sorted(newVariants.keys())
    for chr in chromosomes:
        for pos in sorted(newVariants[chr].keys()):
            if numLineWritten > variantsPerFile:
                # print yes
                idx += 1
                if fileOpen:
                    fo.close()
                fo = open("%s.%d.txt" % (outputPrefix, idx) , 'w')
                fileOpen = True
                numLineWritten = 0
            else:
                numLineWritten += 1
            
            # create Variant objects
            # determine minimum and maximum reference position
            vars = []
            minRef = 2000000000
            maxRef = 0
            for tup in set(newVariants[chr][pos]):
                vars.append(Variant(tup[0],tup[1]))
                totVar += 1
                if vars[-1].refStart<minRef:
                    minRef = vars[-1].refStart
                if vars[-1].refEnd>maxRef:
                    maxRef = vars[-1].refEnd

            leftPos = minRef - hapWidth
            if leftPos<0:
                leftPos = 0
            rightPos = maxRef + hapWidth

            wsize = rightPos - leftPos

            # TODO if window becomes too big, or has too many variants, split it up into multiple windows
            # results should be joined at the mergeOutput stage

            windows.append(wsize)
            if not histWindows.has_key(wsize):
                histWindows[wsize] = 1
            else:
                histWindows[wsize] += 1

            
            numVar = len(vars)
            numWin = numVar/int(maxVarPerWindow)
            lw = len(vars) % int(maxVarPerWindow)

            vc = 0
            finished = False
            for win in range(numWin+lw):
                fo.write("%s %d %d" % ( chr,leftPos, rightPos ) )
                for winIdx in range(maxVarPerWindow):
                    fo.write(" %d,%s" % (vars[vc].refStart, vars[vc].str) )
                    vc += 1

                    if vc == numVar:
                        finished = True
                        break
                fo.write("\n");

                if finished:
                    break


            tot += 1

    fo.close()

    sys.stderr.write("Number of candidates: %d " % totVar)
    sys.stderr.write("Number of windows: %d " % tot)
    sys.stderr.write("Maximum window size: %d " %  (max(windows)) )
    sys.stderr.write("Mean window size: %d\n" % (sum(windows)/len(windows)))
    if False:
        for ws in sorted(histWindows.keys()):
            print " %d:%d" % (ws, histWindows[ws]),
        print "\n"



    return idx 



def split_and_merge(inputVarFile = '', windowFilePrefix = '', minDist = 10, variantsPerFile = 1000, regions = ''):


# get variants
    variants = {}
    numVar = 0
    fvar = open(inputVarFile, 'r')
    numread = 0
    for line in fvar.readlines():
        dat = line.rstrip("\n").split()
        chr = dat[0]
        pos = int(dat[1])
        numread += 1
        if numread % 10000 == 1:
            print numread,"lines read"
        try:
            chrno = int(chr)
            if chrno>24 or chrno == 0:
                sys.stderr.write("Chromosome number>24. Are you sure this is correct?\n")
        except ValueError:
            pass
        
        if not variants.has_key(chr):
            variants[chr] = {}
        
        if not variants[chr].has_key(pos):
            variants[chr][pos] = []

        i = 2    
        while i<len(dat) and dat[i] != '#':    
            variants[chr][pos].append(dat[i])
            numVar += 1
            i += 1
    fvar.close()
    sys.stderr.write("Total variants read: %d\n" % (numVar))

# output statistics

    sys.stderr.write("Number of chromosomes: %d\n" % (len(variants)))
    for chr in variants.keys():
        sys.stderr.write("\tchr %s: %d\n" % ( chr, len(variants[chr]) ))

# per chromosome: group variants and split

    newVariants = {}
    startIdx = 0
    for chr in variants.keys():
        totVar = 0
        positions = sorted(variants[chr].keys())
        newPosition = positions[:]

        done = False
        while not done:
            done = True
            for p in range(1,len(positions)):
                if newPosition[p] != newPosition[p-1] and newPosition[p]-positions[p-1]<=minDist:
                    newPosition[p]=newPosition[p-1]
                    done = False

        # cluster variants
        newVariants[chr] = {}
        for p in range(0, len(newPosition)):
            newPos = newPosition[p]
            pos = positions[p]

            if not newVariants[chr].has_key(newPos):
                newVariants[chr][newPos] = []

            for var in variants[chr][pos]:
                newVariants[chr][newPos].append( (positions[p], var) )
                totVar += 1

                    
        sys.stderr.write("Chromosome: %s Total lines: %d at minimum distance %d\n" % (chr, totVar, minDist))
        sys.stderr.flush()

        old_idx = write_output_candidates(newVariants, outputPrefix = windowFilePrefix, variantsPerFile = variantsPerFile, startIdx = startIdx)
        startIdx = old_idx
        newVariants[chr] = {}

def main(argv):
    parser = OptionParser()
   
    regions = ''

    parser.add_option("-i","--inputVarFile", dest = "varFile", help = "input file with candidate variants")
    parser.add_option("-w","--windowFilePrefix", dest = "windowFilePrefix", help = "prefix of output files with windows") 
    parser.add_option("-m","--minDist", dest = "minDist", help = "mininum distance between two windows", default = 20, type = "int")
    parser.add_option("-n","--numWindowsPerFile", dest = "numWindowsPerFile", help = "number of windows per file", default = 1000, type ="int")
    (options, args) = parser.parse_args() 

    if options.varFile == None:
        sys.stderr.write("Please specify --inputVarFile\n")
        sys.exit(1)
    if options.windowFilePrefix == None:
        sys.stderr.write("Please specify --windowFilePrefix\n")
        sys.exit(1)

    split_and_merge(inputVarFile = options.varFile, windowFilePrefix = options.windowFilePrefix, minDist = options.minDist, variantsPerFile = options.numWindowsPerFile, regions = regions)


if __name__ == "__main__":
    main(sys.argv[1:])

