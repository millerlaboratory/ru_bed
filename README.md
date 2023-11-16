# ru_bed

`Usage: ru_bed.R [options] genelist prefix`

Supply a list of genes separated by '-' and a prefix for the name of your target files.

By default outputs a bedfile for use with adaptive sampling and a tab separated bed file with gene names for downstream use with samtools.

Output is also printed to stdout.

Example:
``` ./makeRUTarget.R F8-F9-NDN testcase```
outputs testcase.targets.bed and testcase.named.targets.bed

Options:
        -c CONTROLS, --controls=CONTROLS
                Controls genes to output, separated by '-'. default COL1A and FMR.

        -b BUFFER, --buffer=BUFFER
                buffer to add to each side of each target. default is 100kb.

        -g, --nonamesave
                Omit saving a copy of bed file with gene names, default false.

        -e ENSEMBL, --ensembl=ENSEMBL
                Ensembl library file (csv). Including one increases performance speed. default included.

        --controlBuffer=CONTROLBUFFER
                Specify control buffer length, if different from gene buffer. default is 50kb.

        -h, --help
                Show help message and exit