# ru_bed

`Usage: Rscript ru_bed.R [options] genelist prefix`

Supply a list of genes separated by '-' and a prefix for the name of your target files.

By default outputs a bedfile for use with adaptive sampling and a tab separated bed file with gene names for downstream use with samtools.

Output is also printed to stdout, along with the total target size and percent of diploid human genome (estiamted at 3.1 GB)

Example:
``` Rscript ru_bed.R F8-F9 testcase``` **OR** ```./ru_bed.sh -t F8-R9 -n testcase```
outputs testcase.targets.bed and testcase.named.targets.bed, as well as:

```
  external_gene_name chromosome_name start_position end_position
1                 F9            chrX      139400000    139700000
2               FMR1            chrX      147850000    148050000
3                 F8            chrX      154700000    155150000
4             COL1A1           chr17       50100000     50300000

[1] "Total target size is 1.15 MB"
[1] "Total percent of diploid human genome is 0.04 percent"
```

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
    