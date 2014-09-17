# dindel-tgi

This is a fork of [dindel-1.01](https://www.sanger.ac.uk/resources/software/dindel/)
by [The Genome Institute at Washington University](http://genome.wustl.edu).
It is licensed under the [GNU GPLv3](http://www.gnu.org/copyleft/gpl.html).

[Dindel](https://www.sanger.ac.uk/resources/software/dindel/) was originally
developed by Cornelis Albers together with Gerton Lunter ([Wellcome Trust Centre for Human Genetics, University of Oxford](http://www.well.ox.ac.uk/home))
and Richard Durbin ([Wellcome Trust Sanger Institute](https://www.sanger.ac.uk)).  

> Dindel is a program for calling small indels from short-read sequence data
('next generation sequence data'). It is currently designed to handle only
Illumina data.
> 
> Dindel takes BAM files with mapped Illumina read data and enables researchers
to detect small indels and produce a VCF file of all the variant calls. It has
been written in C++ and can be used on Linux-based and Mac computers (it has
not been tested on Windows operating systems).
>
> ---
>
> Dindel requires a BAM file containing the read-alignments as input. It then
extracts candidate indels from the BAM file, and realigns the reads to
candidate haplotypes consisting of these candidate indels in windows of ~120
bp. If there is sufficient evidence for an alternative haplotype to the
reference, it will call an indel.
>
> Dindel can test candidate indels discovered with other methods, for instance
longer deletions found by split-read methods or indels obtained through
assembly methods. Dindel will then realign both mapped and unmapped reads to
see if the candidate indel is supported by the reads.
>
> Dindel produces a [VCF file](http://vcftools.sourceforge.net/) with the indel
calls. Genotype likelihoods can be obtained from intermediate files generated
by Dindel.
>
> There is basic support for outputting a realigned BAM file for each
realignment-window. These realigned BAM files can be used to call SNPs near
(candidate) indels.

You can find more information, including a [manual](ftp://ftp.sanger.ac.uk/pub/resources/software/dindel/manual-1.01.pdf)
and the original [source tarball](ftp://ftp.sanger.ac.uk/pub/resources/software/dindel/source_code/dindel-1.01-src.tar.gz),
on the Wellcome Trust Sanger Institute's [Dindel page](https://www.sanger.ac.uk/resources/software/dindel/).
You may also be interested in the *Genome Research* 2010 article **Dindel:
Accurate indel calls from short-read data** (DOI: [10.1101/gr.112326.110](http://dx.doi.org/10.1101%2Fgr.112326.110)).

This fork currently contains a patch to allow sequence data that includes the
`0x800` flag (supplementary alignment) from newer versions of the [SAM](http://samtools.github.io)
specification. It also removes the included copy of Boost and adds [samtools-0.1.19](https://github.com/samtools/samtools/tree/0.1.19)
as a submodule.

## Compiling

This software depends on [Boost](http://www.boost.org), which should probably
be installed using your packager manager of choice.

This software also depends on [SAMtools](http://samtools.github.io) for SAM
format support, and [SeqAn](http://www.seqan.de) for aligning candidate
haplotypes to the reference sequence with the Needleman-Wunsch algorithm. Both
of these dependencies are included in this repo.

If libboost is installed and the samtools submodule is in place, you
should be able to compile the `dindel` binary using `make`.
