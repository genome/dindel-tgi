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

def parseGLFFileLine(dat={}, isHomozygous = False):
        pos = dat [ 'realigned_position' ];
        alleles = ['*'];
        alleles.extend(dat[ 'nref_all'].split(','));
        glfstr = dat ['glf'].split(',');
        indidx = dat ['indidx'];
        
        chr = dat [ 'tid' ];

        # print 'alleles: ', alleles

        genLiks = {}

        for genstr in glfstr:
            gd = genstr.split(':');
            if len(gd)!=2:
                return None


            # print 'gd: ', gd[0], gd[1]
            all = gd[0].split('/');
            # print 'all: ', all[0], all[1], len(all)
            if len(all)!=2 or int(all[0])>=len(alleles) or int(all[1])>=len(alleles):
                return None

            genotype = str(alleles[int(all[0]) ])+'/' + str(alleles[int(all[1])])
            genLiks[genotype]=gd[1];

        return GenotypeLikelihood.GenotypeLikelihood(chr=chr, pos=pos, genLiks=genLiks, miscData =
                dat, isHomozygous = isHomozygous)

def getVCFString(glf = {}, fa = None, maxHPLen = 10, addFilters = [], useFRFilter = 'yes', filterQual = 20):
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

    
    #if glf['vartype']=="snp":
    #    altseq = glf['nref_all'][3]
    #    refseq = fa.get(chr,pos+1,1)[0]
    #    report_pos = pos+1
    #
    #else:
    if True:
        nref_all = glf['nref_all']
        if nref_all.find(',')!=-1:
            raise NameError('Internal error')
        report_pos = pos

        gtAlleles = glf['genotype'].split('/')
        max_del_len = 0
        
        scanAlleles = gtAlleles[:]
        scanAlleles.append(nref_all)
        for gta in scanAlleles:
            var = Variant.Variant(varString = gta)
            if var.type == "del":
                if var.length>max_del_len:
                    max_del_len = var.length
  
        seqlen = 1+max_del_len# 1 base on either side
        refseq = ''.join(fa.get(chr, report_pos, seqlen))
        
 
        vnref = Variant.Variant(varString = nref_all)
        if vnref.type == 'del':
            altseq = refseq[0]+refseq[(1+vnref.length):]
        elif vnref.type == 'ins':
            altseq = refseq[0]+vnref.seq+refseq[1:]
        elif vnref.type == 'snp':
            altseq = refseq[0]+vnref.seq[0]+refseq[2:]
        else:
            raise NameError('Unknown allele')

        # recode genotype

        altseqs = [altseq]
        
        altseq_to_type = {}
        altseq_to_type[altseq] = vnref.type


        rec_gta = []
        for gta in gtAlleles:
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


            rec_gta.append("%d" % g_code)

        rec_gt = '/'.join(rec_gta)


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

    if float(glf['post_prob_all'])<filterQual:
        filters.append('q%d' % filterQual)

    if addFilters != []:
        filters.extend(addFilters)

    if filters== []:
        filterStr = 'PASS'
    else:
        filterStr = ';'.join(filters)
 
    infoStr = "DP=%d;NF=%d;NR=%d;HP=%d" % (int(glf['num_hap_reads']), int(glf['num_cover_forward']), int(glf['num_cover_reverse']), hplen)



    #VCF header string
    #fv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (sampleID))


    rstr = "%s\t%s\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (glf['chr'], report_pos, refseq, ','.join(altseqs), "%1.2f" % (glf['post_prob_all']), filterStr, infoStr,"GT:GQ","%s:%d" % (rec_gt, int(glf['post_prob_geno'])))
    
    return rstr 


def processPooledGLFFiles(bamFilesFile = '', glfFilesFile = '', refFile = '', outputVCFFile = '', maxHPLen = 10, minForwardReverse = 1, minDist = 10, dbSNPWindow = 50, newVarCov = False, doNotFilterOnFR = False, filterQual = 20, numSamples = 1, numBamFiles = 1):
    coverageRange = [20, 10000]

    # read file with glf files
    allFiles = []
    headerLabels = []

    f = open(glfFilesFile,'r')
    idx = 0
    for line in f.readlines():
        idx += 1
        dat = line.rstrip("\n").split()
        for gf in dat:
            if not os.path.exists(gf):
                sys.stderr.write("WARNING: GLF file %s does not exist.\n" % gf)
            else:
                if os.path.splitext(gf)[-1]=='.gz':
                    fgf = gzip.open(gf,'r')
                else:
                    fgf = open(gf,'r')
                line = fgf.readline()
                if line == '':
                    sys.stderr.write("WARNING: GLF file %s is empty.\n" % gf)
                else:
                    d = line.rstrip("\n").split()
                    if headerLabels == []:
                        headerLabels = d[:]
                        allFiles.append(gf)
                    else:
                        if d != headerLabels:
                            sys.stderr.write("Inconsistent header in GLF file %s\n" % gf)
                        else:
                            allFiles.append(gf)
     
                fgf.close()

    
    f.close()

    fa = Fasta.Fasta(fname = refFile)

    # read precall files
    # make hash table [pos][variant][fname]

    numInds = numSamples
    minFreq = 1.0/(float(2*numInds)*5)
    
    nf = 0

    try:
        realpos_col = headerLabels.index('realigned_position')
        var_col = headerLabels.index('nref_all')


        # apply filters across individuals

        tcFilter = "tc%d" % minDist

        col_num_reads = headerLabels.index('num_reads')
            
        col_num_forward_old = headerLabels.index('num_cover_forward')
        col_num_reverse_old = headerLabels.index('num_cover_reverse')
        
        col_num_forward = headerLabels.index('var_coverage_forward')
        col_num_reverse = headerLabels.index('var_coverage_reverse')

        col_post_prob = headerLabels.index('post_prob_variant')
        chr_col = headerLabels.index('tid')
        idx_col = headerLabels.index('indidx')
        ana_col = headerLabels.index('analysis_type')   
    except ValueError:
        raise NameError("GLF files are corrupt. Could not find all required columns.")

    pass_filters = {}
    varStat = {}
    nr = 0
    num_pass = 0
   
    # read depth histo
    rdhist = {}
    for glffile in allFiles:
        fglf = FileUtils.FileWithHeader(fname = glffile, mode = 'r', joinChar = ' ')
        print "Reading", glffile
        done = False
        while True:
            pos = -1
            var = ''
            nr += 1
                     
            if nr % 10000 == 9999:
                print "Number of lines read:",nr+1
            
            num_ind_with_data = 0
            tot_coverage = 0
            tot_num_forward = 0
            tot_num_reverse = 0

            tot_num_forward_old = 0
            tot_num_reverse_old = 0

            skip = False

            for fidx in range(0,numBamFiles):
                try:
                    dat = fglf.readlineList()
                except IOError:
                    sys.stderr.write("WARNING: IOError in %s\n" % glffile)
                    done = True
                    break

                    
                if dat == []:
                    done = True
                    break
                if dat[realpos_col] == 'NA':
                    skip = True
                    break
                if dat[ana_col] != "singlevariant":
                    skip = True
                    break

                if dat[idx_col] != 'NA' and int(dat[idx_col])>= numBamFiles:
                    raise NameError('Error. Is the number of BAM files correctly specified?')

                if pos == -1:
                    pos = int(dat[realpos_col])
                    var = dat[var_col]
                    chr = dat[chr_col]
                else:
                    if int(dat[realpos_col])!=pos:
                        raise NameError('Inconsistent glf files! Is the number of BAM files correctly specified?')
                
                if int(dat[idx_col]) != fidx:
                    sys.stderr.write("Error reading this variant: %s %d %s in %s\n" % (chr,pos,var, glffile))
 
                tot_num_forward_old += int(dat[col_num_forward_old])
                tot_num_reverse_old += int(dat[col_num_reverse_old])
                
                if fidx == 0:
                    # only record for first individual
                    tot_num_forward = int(dat[col_num_forward])
                    tot_num_reverse = int(dat[col_num_reverse])
                    
                numreads = int(dat[col_num_reads])
                if numreads>0:
                    num_ind_with_data += 1
            
                tot_coverage += numreads
            if skip:
                continue
            if done:
                break


            prob = float(dat[col_post_prob])
            freq = float(dat[headerLabels.index('est_freq')])
            if rdhist.has_key(tot_coverage):
                rdhist[tot_coverage] += 1
            else:
                rdhist[tot_coverage] = 1

  
            if prob>0.20:
                if not varStat.has_key(chr):
                    varStat[chr] = {}
                if not varStat[chr].has_key(pos):
                    varStat[chr][pos] = {}
                # hplen
                seq = fa.get(chr, pos+1-25,50)
                hplen = AnalyzeSequence.HomopolymerLength(seq=seq, pos = 25)
                    
               
                 
                varStat[chr][pos][var] = {'QUAL':prob,'NF':tot_num_forward, 'NR':tot_num_reverse, 'NFS':tot_num_forward_old, 'NRS':tot_num_reverse_old, 'DP':tot_coverage, 'NS':num_ind_with_data,'AF':freq,'HP':hplen}

            del dat  
        # finished reading this one
        fglf.close()
    #print "Number of variants passing filters:", num_pass
    
    # apply haplotype coverage and other filters
    
    coverageRange = getPercentiles(rdhist, [1,99])

    fqp = 1.0 - math.pow(10.0, -float(filterQual)/10.0)
    fqp_str = "q%d" % filterQual

    for chr in varStat.keys():
        for pos in varStat[chr].keys():
            for varseq, var in varStat[chr][pos].iteritems():
                
                filters = [] 
                prob = var['QUAL']
                num_ind_with_data = var['NS']
                hplen = var['HP']
                freq = var['AF']
                tot_coverage = var['DP']
                tot_num_forward = var['NF']
                tot_num_reverse = var['NR']
                if prob<fqp:
                    filters.append(fqp_str)
                if (tot_num_forward < minForwardReverse or tot_num_reverse < minForwardReverse) and not doNotFilterOnFR:
                    filters.append('fr0')
                if tot_coverage < coverageRange[0] or tot_coverage>coverageRange[1]:
                    filters.append('ocr')
                if num_ind_with_data<numInds/2:
                    filters.append('s50')
                if hplen>maxHPLen:
                    filters.append("hp%d" % (maxHPLen))
                if freq<minFreq:
                    filters.append("mf")

                if filters == []:
                    if not pass_filters.has_key(chr):
                        pass_filters[chr]={}
                    if not pass_filters[chr].has_key(pos):
                        pass_filters[chr][pos]=[]
                    pass_filters[chr][pos].append(varseq)
                    num_pass += 1

                if filters == []:
                    varStat[chr][pos][varseq]['filter'] = ''
                else:
                    varStat[chr][pos][varseq]['filter'] = ';'.join(filters)    

 
    # now visit each chromosome and apply closeness filter
    chromosomes = [str(c) for c in range(1,23)]
    chromosomes.extend(['X','Y'])

    other_chr = list(set(varStat.keys())-set(chromosomes))
    chromosomes.extend(other_chr)

    # create VCF file
    print "Writing VCF"

    fv = open(outputVCFFile, 'w')
    fv.write("##fileformat=VCFv4.0\n")
    fv.write("##source=Dindel\n")
    fv.write("##reference=%s\n" % refFile)
    fv.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n")
    fv.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total number of reads in haplotype window\">\n")
    fv.write("##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Reference homopolymer tract length\">\n")
    fv.write("##INFO=<ID=NFS,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant site on forward strand\">\n")
    fv.write("##INFO=<ID=NRS,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant site on reverse strand\">\n")
    fv.write("##INFO=<ID=NF,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant on forward strand\">\n")
    fv.write("##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant on reverse strand\">\n")
    fv.write("##INFO=<ID=AF,Number=-1,Type=Float,Description=\"Allele frequency\">\n")
    fv.write("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership build 129 - type match and indel sequence length match within %d bp\">\n" % dbSNPWindow)
    fv.write("##FILTER=<ID=q%d,Description=\"Quality below %d\">\n" % (filterQual, filterQual))
    fv.write("##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n")
    fv.write("##FILTER=<ID=tc%d,Description=\"Indel site was closer than %d base pairs from another site with higher posterior probability\">\n" % (minDist, minDist))
    fv.write("##FILTER=<ID=hp%d,Description=\"Reference homopolymer length was longer than %d\">\n" % (maxHPLen, maxHPLen))
    if not doNotFilterOnFR:
        fv.write("##FILTER=<ID=fr0,Description=\"Non-ref allele is not covered by at least one read on both strands\">\n")
    fv.write("##FILTER=<ID=ocr,Description=\"Number of reads in haplotype window outside coverage range %d %d\">\n" % (coverageRange[0], coverageRange[1]))
    fv.write("##FILTER=<ID=mf,Description=\"Too low non-ref allele frequency\">\n")

    fv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


    for chr in chromosomes:
        if not pass_filters.has_key(chr):
            continue
        # filter out variants that are too close
        totSites = 0
        positions = sorted(pass_filters[chr].keys())
        newPosition = positions[:]

        done = False
        while not done:
            done = True
            for p in range(1,len(positions)):
                if newPosition[p] != newPosition[p-1] and newPosition[p]-positions[p-1]<=minDist:
                    newPosition[p]=newPosition[p-1]
                    done = False

        newSites = {}
        for p in range(0, len(newPosition)):
            newPos = newPosition[p]
            pos = positions[p]

            if not newSites.has_key(newPos):
                newSites[newPos] = {}

            if not newSites[newPos].has_key(pos):
                newSites[newPos][pos]=[]

            for var in varStat[chr][pos].keys():
                newSites[newPos][pos].append(var)

        print "New number of sites:", len(newSites.keys())
        print "Number of sites filtered:",len(pass_filters[chr].keys())-len(newSites.keys())

        # select best call for double sites

        filtered = []
        for newPos in newSites.keys():
            old = newSites[newPos].keys()
            
            pos_probs = []
            pos_vars = []
            pos_pos = []
            for oldPos in old:
                probs = []
                vars = []
                max_prob = -1.0
                max_var = ''
                for var in newSites[newPos][oldPos]:
                    prob=varStat[chr][oldPos][var]['QUAL']
                    if prob>max_prob:
                        max_prob = prob
                        max_var =var
                pos_probs.append(max_prob)
                pos_vars.append(max_var)
                pos_pos.append(oldPos)



            idx = pos_probs.index(max(pos_probs))
            okpos = pos_pos[idx]
            filtered.append(pos_pos[idx])

            for duppos in set(old)-set([okpos]):
                for var in varStat[chr][duppos].keys():
                    
                    if varStat[chr][duppos][var]['filter'] == '':
                        varStat[chr][duppos][var]['filter'] == tcFilter
                    else:
                        varStat[chr][duppos][var]['filter']+=';'+tcFilter

                        
        print "Number of indel sites:",len(filtered)
        
    
        for pos in sorted(varStat[chr].keys()):
            for var in varStat[chr][pos].keys():
                

                indel_report_pos = pos
                #refall = fa.get(chr, pos+1, 1)
                qual = -int(10.0*math.log10(max(1.0-float(varStat[chr][pos][var]['QUAL']),1e-10)))
                infofield = []
                for tag in ['AF','NS','DP','HP','NF','NR','NFS','NRS']:
                    val = (varStat[chr][pos][var][tag])
                    infofield.append("%s=%s" % (tag,val))

                vnref = Variant.Variant(varString = var)
                max_del_len = 0
                if vnref.type == "del":
                    if vnref.length>max_del_len:
                        max_del_len = vnref.length

                seqlen = 1 + max_del_len
                refseq = ''.join(fa.get(chr, indel_report_pos, seqlen))
                if vnref.type == "del":
                    altseq = refseq[0]+refseq[(1+vnref.length):]
                elif vnref.type == "ins":
                    altseq = refseq[0]+vnref.seq+refseq[1:]
                elif vnref.type == "snp":
                    indel_report_pos += 1
                    refseq = refseq[1]
                    altseq = vnref.seq[0]



                infostr = ';'.join(infofield)
                filterstr = varStat[chr][pos][var]['filter']
                if filterstr == '':
                    filterstr = 'PASS'
                id = '.'
                outstr = "%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\n" % (chr, indel_report_pos, id, refseq, altseq, qual, filterstr, infostr) 
                fv.write(outstr)
    fv.close()

def main(argv):
    parser = OptionParser()
    
    regions = ''
    parser.add_option("-i","--inputFiles", dest = "inputFiles", help = "file that contains list of Dindel '.glf.txt' files that should be merged")
    parser.add_option("-o","--outputFile", dest = "outputFile", help = "output VCF file with variant calls") 
    parser.add_option("-r","--refFile", dest = "refFile", help = "reference sequence _indexed_ Fasta file")
    #parser.add_option("-t","--type", dest = "type", help = "type of input files. Choices for TYPE: diploid (single diploid samples), haploid (results in homozygous genotypes) or pool (for pools)", choices= typeOptions)
    parser.add_option("--numSamples", dest = "numSamples", help = "number of samples")
    parser.add_option("--numBamFiles", dest = "numBAMFiles", help = "number of BAM files", default = 1)

    parser.add_option("--maxHPLen", dest = "maxHPLen", help = "maximum length of homopolymer run to call an indel", default = 10, type ="int")
    #parser.add_option("--newVarCov", dest = "newVarCov", action = "store_true", default = False)
    parser.add_option("--filterFR", dest = "filterFR", help = "filter on forward/reverse count of reads (stringent)", action = "store_true", default = False)
    parser.add_option("--filterQual", dest = "filterQual", help = "quality below which variants are filtered", default = 20)

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

    if options.numSamples == None:
        sys.stderr.write("Please specify --numSamples option\n")
        sys.exit(1)

    processPooledGLFFiles(glfFilesFile = options.inputFiles, maxHPLen = options.maxHPLen, refFile = options.refFile, outputVCFFile = options.outputFile, doNotFilterOnFR = (not options.filterFR), filterQual = int(options.filterQual), numSamples = int(options.numSamples), numBamFiles = int(options.numBAMFiles))




if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except:
        sys.stderr.write("An error occurred!\n")
        raise
        sys.exit(1)


