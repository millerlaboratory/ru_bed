#!/bin/bash

Help()
{
  echo ""
  echo "makeNamedFile.sh"
  echo ""
  echo "Generates needed analysis files from ONT supplied six column bed files with two targets for each gene"
  echo ""
  echo "Default behavior saves a copy of named and unbuffered files to /n/alignments/bed_targets/"
  echo "disable with flag -s"
  echo ""
  echo "supply path to bed file with -b"
  echo "supply name prefix with -n"
  echo ""
  echo "Usage: ./makeNamedFile.sh -b example.bed -n testcase"
  echo ""
  echo "example generates:"
  echo "testcase.naked.named.targets.bed"
  echo "testcase.named.targets.bed"
  echo ""
}

temp=$(mktemp)-ru
canonlibrary="resources/ensembl.gene.library.canonical.tsv"
synlibrary="resources/ensembl.gene.library.synonyms.tsv"
pseudolibrary="resources/ensembl.gene.library.pseudo.canon.tsv"
pseudosynlibrary="resources/ensembl.gene.library.pseudo.synonym.tsv"
savenetwork=1

fieldstring=6,7,8

while getopts "n:hb:s" option; do
  case $option in
    n) NAME=$OPTARG ;;
    h) Help
       exit;;
    b) BEDFILE=$OPTARG ;;
    s) savenetwork=0;;
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

if [ -z ${BEDFILE+x} ]
then
    echo ""
    echo "-b is required, please specify bed file input"
    Help
    echo ""
    exit 1
fi

# make naked target file

targets=( $(cut -f4 $BEDFILE | tr ',' '\n' | sort | uniq ) )

tempref=$(mktemp -t tmp.ruXXXXXref.tsv)
temperrs=$(mktemp -t tmp.ruXXXXXerrors.tsv)
tempunmodified=$(mktemp -t tmp.ruXXXXXallgenesunmodified.tsv)

for gene in ${targets[@]}
do
    #if first four characters are ENSG look up by column 4, the ENSG identifier, not the gene name
    idfield=1
    if [ ${gene[@]:0:4} == "ENSG" ]
    then
        echo "Warning: trying ENSG identifier to locate $gene"
        echo "Warning: trying ENSG identifier to locate $gene" >> $temperrs
        idfield=4
    fi
    tempgeneout=$(mktemp -t tmp.ruXXXXXgeneout.tsv)  

    #first search library of genes

    #grep $gene $canonlibrary | awk -v gene=$gene -v id=$idfield -F '\t' '{if($(id)==gene){print $6,$7,$8,$(id)}}' | tr ' ' '\t' > $tempgeneout
    grep $gene $canonlibrary | awk -v gene=$gene -v id=$idfield -v fieldstring=$fieldstring -F '\t' 'BEGIN{split(fieldstring,fieldlist,",")}{if($(id)==gene){for(i in fieldlist){printf $fieldlist[i]" "};printf$(id)" \n";}}' | tr ' ' '\t' > $tempgeneout
    
    results=$( wc -l $tempgeneout | tr -s ' ' | cut -d ' ' -f1 )
    if [[ $results -gt 0 ]]
    then
        cat $tempgeneout >> $tempref
        if [[ $results -gt 1 ]]
        then
            echo "Warning: multiple targets match $gene"
            echo "Warning: multiple targets match $gene" >> $temperrs
        fi
    else
        # next search synonyms of genes
        grep $gene $synlibrary | awk -v gene=$gene -v id=$idfield -v fieldstring=$fieldstring -F '\t' 'BEGIN{split(fieldstring,fieldlist,",")}{if($(id)==gene){for(i in fieldlist){printf $fieldlist[i]" "};printf$(id)" \n";}}' | tr ' ' '\t' > $tempgeneout
        results=$( wc -l $tempgeneout | tr -s ' ' | cut -d ' ' -f1 )
        if [[ $results -gt 0 ]]
        then
            cat $tempgeneout >> $tempref
        else
            # next search pseudogenes
            grep $gene $pseudolibrary | awk -v gene=$gene -v id=$idfield -v fieldstring=$fieldstring -F '\t' 'BEGIN{split(fieldstring,fieldlist,",")}{if($(id)==gene){for(i in fieldlist){printf $fieldlist[i]" "};printf$(id)" \n";}}' | tr ' ' '\t' > $tempgeneout
            results=$( wc -l $tempgeneout | tr -s ' ' | cut -d ' ' -f1 )
            if [[ $results -gt 0 ]]
            then
                cat $tempgeneout >> $tempref
            else
                # lastly search synonyms of pseudogenes
                grep $gene $pseudosynlibrary | awk -v gene=$gene -v id=$idfield -v fieldstring=$fieldstring -F '\t' 'BEGIN{split(fieldstring,fieldlist,",")}{if($(id)==gene){for(i in fieldlist){printf $fieldlist[i]" "};printf$(id)" \n";}}' | tr ' ' '\t' > $tempgeneout
                results=$( wc -l $tempgeneout | tr -s ' ' | cut -d ' ' -f1 )
                if [[ $results -gt 0 ]]
                then
                    cat $tempgeneout >> $tempref
                else
                    # log error that gene was not located
                    echo "Error: $gene could not be located"
                    echo "Error: $gene not located" >> $temperrs
                fi
            fi
        fi
    fi
    rm $tempgeneout
done

cat $tempref | sort -k1,1 -k2,2n -k3,3n | uniq > $tempunmodified

nakednamedoutput=$NAME.naked.named.targets.bed

cat $tempunmodified | tr ' ' '\t' > $nakednamedoutput

# make named target file
# the way this is implemented assumes that there is an equivalent '+' target for every '-' target, which appears to be true for ONT's supplied bed files
namedoutput=$NAME.named.targets.bed

sort $BEDFILE -k4,4 -k6,6 | awk 'NR%2==1{start=$2}NR%2==0{print $1,start,$3,$4}' | tr ' ' '\t' | sort -k1,1 -k2,2n -k3,3n > $namedoutput

if [[ $savenetwork -eq 1 ]]
then
    cp $namedoutput /n/alignments/bed_targets/
    cp $nakednamedoutput /n/alignments/bed_targets/
fi

if [ -f $temperrs ]
then
    cp $temperrs $NAME.errors.txt
fi

rm /tmp/tmp.ru*tsv
exit 0