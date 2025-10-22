You can build a [cameo](https://github.com/tobiasrausch/cameo) singularity container (SIF file) using

`sudo singularity build cameo.sif cameo.def`

Once you have built the container you can run analysis using

`singularity exec cameo.sif cameo pileup -g ref.fa input.bam`
