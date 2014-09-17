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

def getPercentiles(hist = {}, pctiles = [ 1,5,10,25, 50, 75, 90, 95, 99]):
    vals = sorted(hist.keys())
    cum = {}

    prevk = 0
    for idx, k in enumerate(vals):
        cum[k]=hist[k]
        #print "cum:",k,hist[k]
        if idx>0:
            cum[k] += cum[prevk]
        prevk = k

    if not cum.has_key(0):
        cum[0] = 0
    tot = cum[prevk]
    iles = pctiles
    ilidx = 0
    iles_k = [0]*len(iles)
    num_pct = len(iles)
    for idx, k in enumerate(vals):
        if ilidx<num_pct and cum[k]>float(iles[ilidx])/100.0*float(tot):
            iles_k[ilidx]=k
            ilidx += 1

    return iles_k


def getVCFString(glf = {}, fa = None, maxHPLen = 10, addFilters = [], useFRFilter = 'no', filterQual = 0):
    # create VCF string for variant

    #fv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (sampleID))
    filters = []
    if useFRFilter == 'yes':
        if not glf['pass_filterFR']:
            filters.append('fr0')

    pos = int(glf['pos'])
    chr = glf['chr']
    seq = fa.get(chr, pos+1-25,50)
    hplen = AnalyzeSequence.HomopolymerLength(seq=seq, pos = 25)

    if True:
        report_pos = pos

        max_del_len = 0
        
        scanAlleles = list(set(glf['nref_all']))
        for gta in scanAlleles:
            var = Variant.Variant(varString = gta)
            if var.type == "del":
                if var.length>max_del_len:
                    max_del_len = var.length
  
        seqlen = 1+max_del_len# 1 base on either side
        refseq = ''.join(fa.get(chr, report_pos, seqlen))
        
 
        # recode genotype

        altseqs = []
        
        altseq_to_type = {}

        rec_gta = []
        for gta in glf['nref_all']:
            g_code = -1
            vnref = Variant.Variant(varString = gta)
            if vnref.type == 'del':
                g_altseq = refseq[0]+refseq[(1+vnref.length):]
            elif vnref.type == 'ins':
                g_altseq = refseq[0]+vnref.seq+refseq[1:]
            elif vnref.type == 'snp':
                g_altseq = refseq[0]+vnref.seq[0]+refseq[2:]
            elif vnref.type == 'ref':
                g_altseq = refseq[:]
                g_code = 0
            else:
                raise NameError('Unknown allele')


            if g_code == -1:
                if g_altseq not in altseqs:
                    altseqs.append(g_altseq)
                    altseq_to_type[g_altseq] = vnref.type
                
                g_code = altseqs.index(g_altseq)+1

        gtd = glf['genotype'].split(':')
        rec_gt = "%s:%d" % (gtd[0], int(float(gtd[1])))


    # check if there is an indel in the altseqs

    onlySNPs = True
    for alts in altseqs:
        if altseq_to_type[alts] != "snp":
            onlySNPs = False

    # if there are no indels, only SNPs, then move position by one

    if onlySNPs:
        report_pos += 1

        refseq = ''.join(fa.get(chr,report_pos,1))
        tmp_altseqs = altseqs[:]
        altseqs = []
        for alt in tmp_altseqs:
            altseqs.append(alt[1:])
        del tmp_altseqs
	
	
    if hplen>maxHPLen:
        filters.append("hp%d" % maxHPLen)

    if glf['qual']<filterQual:
        filters.append("q%d" % filterQual)
    
    # replace 'D' in alt

    tmp_altseqs = []
    for alt in altseqs:
        if alt.find('D')!=-1:
            tmp_altseqs.append('<DEL>')
        else:
            tmp_altseqs.append(alt)

    altseqs=tmp_altseqs[:]
    del tmp_altseqs

    if addFilters != []:
        filters.extend(addFilters)

    if filters== []:
        filterStr = 'PASS'
    else:
        filterStr = ';'.join(filters)
 
    infoStr = "DP=%d;NF=%d;NR=%d;NRS=%d;NFS=%d;HP=%d" % (int(glf['num_hap_reads']), int(glf['num_cover_forward']), int(glf['num_cover_reverse']), int(glf['num_cover_forward_old']), int(glf['num_cover_reverse_old']),hplen)



    #VCF header string
    #fv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (sampleID))


    rstr = "%s\t%s\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (glf['chr'], report_pos, refseq, ','.join(altseqs), "%s" % (glf['qual']), filterStr, infoStr,"GT:GQ","%s" % (rec_gt))
    
    return (rstr, report_pos)


def processDiploidGLFFile(glfFile = '', variants = {}, refFile = '', maxHPLen = 10, isHomozygous = False, doNotFilterOnFR = False, newVarCov = False, filterQual = 20):


    # setup reference sequence

    fa = Fasta.Fasta(fname = refFile)

    # variants will be added to variants

    fglf = FileUtils.FileWithHeader(fname = glfFile)

    numSkipped = 0 # number of windows that were skipped by Dindel

    # read line by line, aggregate results for identical windows

    prevPos = -1
    prevChr = -1
    prevDat = {}
    while True:
        dat = fglf.readline()
        if dat == {}:
            break
        if True:
            errcode = dat['msg']
            index = dat['index'] # index of window in original variant file

            if errcode != "ok":
                numSkipped += 1
                continue

            if dat['analysis_type'] != 'dip.map':
                continue

            if dat['was_candidate_in_window'] != '1':
                continue

            
            glf = {}
            chrom = dat['tid']
            if chrom != prevChr:
                prevPos = -1
                prevChr = chrom

            glf['chr'] = dat['tid']
            glf['pos'] = dat['realigned_position']
           
            pos = int(glf['pos'])
            
            prevDat = dat
            glf['qual'] = int(float(dat['qual']))
            
            if float(glf['qual'])<1.0:
                continue
     
            glf['nref_all'] = dat['nref_all'].split(',')
            if glf['nref_all'] == ['R=>D']:
                continue
            nfa = dat['var_coverage_forward'].split(',')
            nra = dat['var_coverage_reverse'].split(',')
            ai = 0
            
            glf['num_cover_forward'] = int(nfa[ai])
            glf['num_cover_reverse'] = int(nra[ai])

            glf['num_cover_forward_old'] = int(dat['num_cover_forward'])
            glf['num_cover_reverse_old'] = int(dat['num_cover_reverse'])


            glf['num_hap_reads'] =  dat['num_reads']
            glf['genotype'] = dat['glf']

            (vcf_str,report_pos) = getVCFString(glf = glf, fa = fa, filterQual = filterQual)
            if not variants.has_key(chrom):
                variants[chrom] = {}

            if not variants[chrom].has_key(report_pos):
                variants[chrom][report_pos] = []

            variants[chrom][report_pos].append(vcf_str)

            prevPos = pos

def mergeOutput(glfFilesFile = '', sampleID = 'SAMPLE', refFile = '', maxHPLen = 10, vcfFile = '', isHomozygous = False, newVarCov = False, doNotFilterOnFR = False, filterQual = 20):

    # read list of files and group according to chromosome
    # NOTE assumes at most one chromosome per GLF file is analysed
    
    # contains list of chromosome to file mapping
    chrToFiles = {}

    # check files

    fg = open(glfFilesFile,'r')
    lineidx = 0
    okFiles = []
    for line in fg.readlines():
        lineidx += 1
        dat = line.rstrip("\n").split()
        if len(dat)>1:
            sys.stderr.write("WARNING: additional columns in line %d of file %s were ignored\n" % (lineidx, glfFilesFile))
        
        # determine which chromosome the file has
        fn = dat[0]
        if os.path.exists(fn):
            okFiles.append(fn)
        else:
            sys.stderr.write("File %s does not exist\n. Aborting.\n" % fn)
            sys.exit(1)

    
    fg.close()
            
    sys.stdout.write("Number of non-empty GLF files: %d\n" % len(okFiles))

    # open VCF file

    fv = open(vcfFile, 'w')
    fv.write("##fileformat=VCFv4.0\n")
    fv.write("##source=Dindel\n")
    fv.write("##reference=%s\n" % refFile)
    fv.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total number of reads in haplotype window\">\n")
    fv.write("##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Reference homopolymer tract length\">\n")
    fv.write("##INFO=<ID=NF,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant on forward strand\">\n")
    fv.write("##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant on reverse strand\">\n")
    fv.write("##INFO=<ID=NFS,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant site on forward strand\">\n")
    fv.write("##INFO=<ID=NRS,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant site on reverse strand\">\n")
    fv.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    fv.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n")
    fv.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
    fv.write("##FILTER=<ID=q%d,Description=\"Quality below %d\">\n" % (filterQual, filterQual))
    fv.write("##FILTER=<ID=hp%d,Description=\"Reference homopolymer length was longer than %d\">\n" % (maxHPLen, maxHPLen))
    fv.write("##FILTER=<ID=fr0,Description=\"Non-ref allele is not covered by at least one read on both strands\">\n")
    fv.write("##FILTER=<ID=wv,Description=\"Other indel in window had higher likelihood\">\n")


    fv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (sampleID))


    # get calls
    chromosomes = [str(v) for v in range(1,23)]
    chromosomes.extend(['X','Y'])

    # add misc chromosomes (MT etc)
    variants = {}

    for gf in okFiles:
        sys.stdout.write("Calling variants from GLF file %s\n" % (gf))
        processDiploidGLFFile(glfFile = gf, variants = variants, refFile = refFile, maxHPLen = maxHPLen, isHomozygous = isHomozygous, newVarCov = newVarCov, doNotFilterOnFR = doNotFilterOnFR, filterQual = filterQual)
    
    this_chr = chromosomes[:]
    for chr in variants.keys():
        if chr not in this_chr:
            this_chr.append(chr)

    for chr in this_chr:
        if variants.has_key(chr):
            for pos in sorted(variants[chr].keys()):
                for vcfLine in variants[chr][pos]:
                    fv.write("%s\n" % (vcfLine))

    fv.close()

def main(argv):
    parser = OptionParser()

    parser.add_option("-i","--inputFiles", dest = "inputFiles", help = "file that contains list of Dindel '.glf.txt' files that should be merged")
    parser.add_option("-o","--outputFile", dest = "outputFile", help = "output VCF file with variant calls") 
    parser.add_option("-s","--sampleID", dest = "sampleID", help = "sampleID to be recorded in VCF file [only for --type diploid]", default = "SAMPLE")
    parser.add_option("-r","--refFile", dest = "refFile", help = "reference file")
    parser.add_option("--maxHPLen", dest = "maxHPLen", help = "maximum length of homopolymer run to call an indel", default = 10, type ="int")
    parser.add_option("-f", "--filterQual", dest = "filterQual", help = "variants below this quality will be filtered", default = 20, type = "int")

    (options, args) = parser.parse_args() 

    if options.inputFiles == None:
        sys.stderr.write("Please specify --inputFiles\n")
        sys.exit(1)
    if options.outputFile == None:
        sys.stderr.write("Please specify --outputFile\n")
        sys.exit(1)
    if options.refFile == None:
        sys.stderr.write("Please specify --refFile\n")
        sys.exit(1)

    mergeOutput(glfFilesFile = options.inputFiles, sampleID = options.sampleID, maxHPLen = options.maxHPLen, refFile = options.refFile, vcfFile = options.outputFile, filterQual = int(options.filterQual))


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except:
        sys.stderr.write("An error occurred!\n")
        raise
        sys.exit(1)


