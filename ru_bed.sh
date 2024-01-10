#!/bin/bash

module load Rtools/1.0

Help()
{
  echo ""
  echo "A convenience wrapper for ru_bed.R"
  echo "Supply a list of genes separated by '-' and a prefix for the name of your target files."
  echo "By default outputs a bedfile for use with adaptive sampling and a tab separated bed file with gene names for downstream use with samtools."
  echo "Usage: ./ru_bed.sh -t F8-F9-NDN -n testcase"
  echo ""
  echo "-n : prefix used for output files"
  echo "-t : gene names to target, separated by the '-' character"
  echo "-c : Controls genes to output, separated by '-'. default COL1A and FMR."
  echo "-b : buffer to add to each side of each target. default is 100kb."
  echo "-g : Omit saving a copy of bed file with gene names, default false."
  echo "-e : Ensembl library file (csv). Including one increases performance speed. default included."
  echo "-B : Specify control buffer length, if different from gene buffer. default is 50kb."
  echo "-h: view this help and exit"
  echo ""
  echo ""
}

OPTIONS=""

while getopts "n:t:c:b:ge:B:h" option; do
  case $option in
    n) NAME=$OPTARG ;;
    h) Help
       exit;;
    t) TARGETS=$OPTARG;;
	c) OPTIONS="$OPTIONS -c $OPTARG";;
	b) OPTIONS="$OPTIONS -b $OPTARG";;
    g) OPTIONS="$OPTIONS -g";;
    e) OPTIONS="$OPTIONS -e $OPTARG";;
    B) OPTIONS="$OPTIONS -B $OPTARG";;
  esac
done

if [ -z ${NAME+x} ]
then
    echo ""
    echo "-n is required, please specify a name for output"
    Help
    echo ""
    exit 1
fi

if [ -z ${TARGETS+x} ]
then
    echo ""
    echo "-t is required, please specify at least one target gene"
    Help
    echo ""
    exit 1
fi

COMMAND="Rscript ru_bed.R $OPTIONS $TARGETS $NAME.temp"

echo "running $COMMAND"
echo ""
$COMMAND

cat $NAME.temp.targets.bed | tr -d ' '| tr ',' '\t' > $NAME.targets.bed
rm $NAME.temp.targets.bed
cat $NAME.temp.named.targets.bed | tr -d ' ' | tr ',' '\t' > $NAME.named.targets.bed
rm $NAME.temp.named.targets.bed