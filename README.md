# Cameo: Cancer Analysis of Methylation

## Installing cameo

Cameo is available as a [statically linked binary](https://github.com/tobiasrausch/cameo/releases/), as a minimal [docker container](https://hub.docker.com/r/trausch/cameo/) or as a [singularity containter (SIF file)](https://github.com/tobiasrausch/cameo/releases/).

## Building from source

Cameo can be built from source using a recursive clone and make. Cameo depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/cameo.git`

`cd cameo/`

`make all`

## Pileup tables of 5mC and 5hmC

Cameo assumes coordinate sorted BAM files with MM and ML tags for DNA modifications. To identify 5mC and 5hmC modifications by strand from such a BAM file:

`cameo pileup -g hg38.fa input.bam`

`cameo pileup -o out.table.tsv -g hg38.fa input.bam`

To aggregate modifications only at CpG sites and independent of strand:

`cameo pileup -cp -g hg38.fa input.bam`

To calculate methylation over intervals such as promoter methylation:

`cameo pileup -cp -f promoter.bed -g hg38.fa input.bam`

To identify (de)methylated regions in a BAM file:

`cameo pileup -cpk -g hg38.fa input.bam`

## Credits

[HTSlib](https://github.com/samtools/htslib) is heavily used for alignment processing and [Boost](https://www.boost.org/) for various data structures and algorithms.

