import sys, os

class FastaTarget:
    def __init__(self):
        self.tid = ''
        self.len = ''
        self.offset = -1
        self.blen = -1
        self.llen = -1

class FastaIndex:
    def __init__(self, fname=''):
        self.f = open(fname, 'r')

        self.ft = {}
        for line in self.f.readlines():
            dat = line.split()
            if len(dat)==5:
                tid = dat[0]
                self.ft[tid] = FastaTarget()
                self.ft[tid].tid = tid
                self.ft[tid].len = int(dat[1])
                self.ft[tid].offset = int(dat[2])
                self.ft[tid].blen = int(dat[3])
                self.ft[tid].llen = int(dat[4])
                # print 'Fasta: ', tid, int(dat[1]), int(dat[2]), int(dat[3]), int(dat[4])        
        self.f.close()

class Fasta:
    def __init__(self, fname = '/nfs/users/nfs_c/caa/s103/ref/human_b36_male.fa'):
        self.fa = open(fname,'r')
        self.fai = FastaIndex(fname+'.fai')

    def get(self, tid, pos1based, len):
        pos = pos1based - 1
        try:
            idx = self.fai.ft[tid]
        except KeyError:
            print 'KeyError: ', tid
            raise NameError('KeyError')
        fpos = idx.offset+ ( int(pos)/idx.blen)*idx.llen + (int(pos)%idx.blen)
        self.fa.seek(fpos,0)

        numread=0
        seq = []
        while numread<len:
            char = self.fa.read(1)
            if char!='\n':
                seq.append(char)
                numread +=1
        return seq

def getChromosomes(faFile = ''):

    faiFile = "%s.fai" % (faFile)

    if not os.path.exists(faiFile):
        raise NameError("Cannot find fai file for %s" % faFile)
    
    faidx = FastaIndex(faiFile)

    fachr = faidx.ft.keys()

    chromosomes = []
    autosomal = ["%d" % c for c in range(1,23)]
    autosomal.extend(['X','Y'])

    for chrom in autosomal:
        if chrom in fachr:
            chromosomes.append(chrom)

    chromosomes.extend( list (set(fachr) - set(autosomal)))
    return chromosomes


