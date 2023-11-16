# ru_bed

`Usage: ru_bed.R [options] genelist prefix`

Supply a list of genes separated by '-' and a prefix for the name of your target files.

By default outputs a bedfile for use with adaptive sampling and a tab separated bed file with gene names for downstream use with samtools.

Output is also printed to stdout.

Example:
``` Rscript ru_bed.R F8-F9-NDN testcase``` **OR** ```./ru_bed.sh -t F8-R9-NDN -n testcase```
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

If running with shell script rather than Rscript, do not supply positional arguments, use the following flags:

        -t names of genes to target, separated by '-'
        
        -n prefix for named files

The long form of flags are also not available in the shell script (use -c, not --controls).
    