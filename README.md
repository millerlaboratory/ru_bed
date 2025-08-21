# Adaptive Sampling Utilities

## ru_bed

`Usage: bash ru_bed_v3.sh [options]`

Supply a list of genes separated by ',' and a prefix for the name of your target files. 

Each gene is buffered (by default 50kb) and overlapping regions are merged.

By default outputs a bedfile for use with adaptive sampling, a tab separated bed file with gene names for downstream use with samtools, and a bed file containing unbuffered gene targets.

Default behavior saves a copy of named and unbuffered files to /n/alignments/bed_targets/ -- disable this with flag `-s`.

Targets that could not be located are output to ${prefix}.errors.txt and an error is printed to STDOUT

If more than one target region matches the supplied genenames, both are saved and a warning is printed.

Targets are no longer rounded to the nearest 50kb -- to enable this, use flag `-R`

Total target size and percent of diploid human genome (estiamted at 3.1 GB) is output to screen for both buffered and unbuffered input


### Example:

``` bash ru_bed_v3.sh -t F8,F9,SNURF,SNRPN -n testcase -S```

outputs files named `testcase.targets.bed`, `testcase.named.targets.bed`, and `testcase.naked.named.targets.bed`

```
using default controls



total unbuffered target size (Mb): .46
total unbuffered target genome percentage: .01 %
total final target size (Mb): 2.07
total final target genome percentage: .06 %

```

```
> cat testcase.named.targets.bed

chr15   23585400        23787305        NDN
chr15   24723637        25078723        SNRPN.SNURF
chr17   50084101        50301632        COL1A1
chrX    139430739       139663459       F9
chrX    147811919       148051125       FMR1
chrX    154735788       155126940       F8

> cat testcase.targets.bed

chr15   23585400        23687305        NDN     .       +
chr15   23685400        23787305        NDN     .       -
chr15   24723637        24978723        SNRPN.SNURF     .       +
chr15   24823637        25078723        SNRPN.SNURF     .       -
chr17   50084101        50201632        COL1A1  .       +
chr17   50184101        50301632        COL1A1  .       -
chrX    139430739       139563459       F9      .       +
chrX    139530739       139663459       F9      .       -
chrX    147811919       147951125       FMR1    .       +
chrX    147911919       148051125       FMR1    .       -
chrX    154735788       155026940       F8      .       +
chrX    154835788       155126940       F8      .       -

> cat testcase.naked.named.targets.bed 

chr15   23685400        23687305        NDN
chr15   24823637        24978723        SNRPN
chr15   24954987        24977850        SNURF
chr17   50184101        50201632        COL1A1
chrX    139530739       139563459       F9
chrX    147911919       147951125       FMR1
chrX    154835788       155026940       F8

```
### Required Options

```
        -t TARGETS
                genes to target, separated by ','

        -n NAME
                output prefix/name
```

### Options

```
        -c CONTROLS
                Controls genes to output, separated by ','. default COL1A1,FMR1,NDN

        -b BUFFER
                buffer to add to each side of each target. default is 100kb.


        -B CONTROLBUFFER
                Specify control buffer length, if different from gene buffer. default is 100kb.

        -h
                Show help message and exit

        -N Naked run
                No buffer is used and no controls are used. Overwrites other options

        -s networksave
                do not save files to network locations

        -S stranded
                apply buffer offset on strands (buffer applied upstream on + strand and downstream on - strand)

        -R round
                round to nearest 50kb when printing buffered targets

        -m  manual targets
                add targets  that can not be found in the ensembl database (e.g. gene clusters)

                enter each as a coordinate,name pair and separate targets by semicolons.

                these regions must not contain dashes in the regionnames, please replace with underscores

                example: -m chrX:12345-23456,random;chrY:23455-23456,random2"

```     

## makeNamedFile.sh

`Usage: bash makeNamedFile.sh [options]`

Use this to generate analysis files from ONT's six column bedfiles

Input bedfiles are expected to have two targets for each gene, one on the positive strand and one on the negative strand.

### Required Options

```
        -b BEDFILE
                Path to befile (required)

        -n NAME
                output prefix/name
       
```

### Options

```
        -h
                Show help message and exit


        -s networksave
                do not save files to network locations

```    

### Example:

``` bash makeNamedFile.sh -b ontBedFile.bed -n testcase```

outputs files named `testcase.named.targets.bed`, and `testcase.naked.named.targets.bed`

For original file input:
```
chr1    149773075       149853125       ENSG00000301450,H3C14,H3C15     .       +
chr1    149793075       149873125       ENSG00000301450,H3C14,H3C15     .       -
chr1    18610846        18748866        PAX7    .       +
chr1    18630846        18768866        PAX7    .       -
chr1    16998664        17054151        SDHB    .       +
chr1    17018664        17074151        SDHB    .       -
chr1    15816095        15940456        SPEN    .       +
chr1    15836095        15960456        SPEN    .       -
```

The contents of `testcase.named.targets.bed` are collapsed intervals that preserve overlap:

```
chr1    15816095        15960456        SPEN
chr1    16998664        17074151        SDHB
chr1    18610846        18768866        PAX7
chr1    149773075       149873125       ENSG00000301450,H3C14,H3C15
```

The contents of `testcase.naked.named.targets.bed` are the ensembl coordinate for each gene:

```
chr1    15836095        15940456        SPEN
chr1    17018664        17054151        SDHB
chr1    18630846        18748866        PAX7
chr1    149793075       149825368       ENSG00000301450
chr1    149840687       149841208       H3C14
chr1    149852608       149853125       H3C15
```
