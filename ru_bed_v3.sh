#!/bin/bash

Help()
{
  echo ""
  echo "ru_bed_v3.sh"
  echo "Supply a list of genes separated by ',' and a prefix for the name of your target files."
  echo "By default outputs a bedfile for use with adaptive sampling, a tab separated bed file with gene names \
  for downstream use with samtools, and a bed file containing unbuffered gene targets."
  echo ""
  echo "Default behavior saves a copy of named and unbuffered files to /n/alignments/bed_targets/"
  echo "disable with flag -s"
  echo "Usage: ./ru_bed.sh -t F8,F9,NDN -n testcase"
  echo ""
  echo "-n : prefix used for output files"
  echo "-t : gene names to target, separated by the ',' character"
  echo "-c : Controls genes to output, separated by ','. default COL1A1,FMR1,NDN"
  echo "-b : buffer to add to each side of each target. default is 100kb."
  echo "-B : Specify control buffer length, if different from gene buffer. default is 50kb"
  echo "-N : Naked run. No buffer is used and no controls are used. Overwrites other options"
  echo "-h : view this help and exit"
  echo "-s : do not save files to network, only to supplied output path that is part of prefix"
  echo "-S : output strand in bed file (disabled by default)"
  echo "-m : manual targets to add that can not be found in the ensembl database (e.g. gene clusters) \
             enter each as a coordinate,name pair and separate targets by semicolons. \
             these regions must not contain dashes in the regionnames, please replace with underscores \
             example: -m chrX:12345-23456,random;chrY:23455-23456,random_2"
  echo ""
  echo ""
}

naked=0
controls=( FMR1 COL1A1 NDN )
temp=$(mktemp)-ru
canonlibrary="resources/ensembl.gene.library.canonical.tsv"
synlibrary="resources/ensembl.gene.library.synonyms.tsv"
pseudolibrary="resources/ensembl.gene.library.pseudo.canon.tsv"
pseudosynlibrary="resources/ensembl.gene.library.pseudo.synonym.tsv"
useENSG=0
buffer=100000
controlbuffer=100000
resolution=50000
round=0
savenetwork=1
STRANDED=0

while getopts "n:t:c:b:B:hNRsm:S" option; do
  case $option in
    n) NAME=$OPTARG ;;
    h) Help
       exit;;
    t) TARGETSTRING=$OPTARG;;
	c) controlstring="$OPTARG";;
	b) buffer="$OPTARG";;
    B) controlbuffer="$OPTARG";;
    N) naked=1;;
    R) round=1;;
    s) savenetwork=0;;
    S) STRANDED=1;;
    m) manualtargets="$OPTARG";;
  esac
done

echo ""

if [[ $STRANDED -eq 1 ]]
then
    canonlibrary="resources/ensembl.gene.library.canonical.stranded.tsv"
    synlibrary="resources/ensembl.gene.library.synonyms.stranded.tsv"
    pseudolibrary="resources/ensembl.gene.library.pseudo.canon.stranded.tsv"
    pseudosynlibrary="resources/ensembl.gene.library.pseudo.synonym.stranded.tsv"
    fieldstring=6,7,8,14
else
    fieldstring=6,7,8
fi

if [ -z ${NAME+x} ]
then
    echo ""
    echo "-n is required, please specify a name for output"
    Help
    echo ""
    exit 1
fi

if [ -z ${TARGETSTRING+x} ]
then
    echo ""
    echo "-t is required, please specify at least one target gene"
    Help
    echo ""
    exit 1
fi

if [ -z ${controlstring+x} ]
then
    echo -e "using default controls\n"
else
    controls=( $(echo $controlstring | tr ',' ' ') )
fi

targets=( $( echo $TARGETSTRING | tr ',' ' ' ) )

allgenes=( ${targets[@]} ${controls[@]} )
tempref=$(mktemp -t tmp.ruXXXXXref.tsv)
temperrs=$(mktemp -t tmp.ruXXXXXerrors.tsv)

# search exhaustively for target

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

## add controls
temprefcontrol=$(mktemp -t tmp.ruXXXXXrefcontrol.tsv)

for gene in ${controls[@]}
do
    tempgeneout=$(mktemp -t tmp.ruXXXXXgeneout.tsv)
    grep $gene $canonlibrary | awk -v gene=$gene -v id=1 -v fieldstring=$fieldstring -F '\t' 'BEGIN{split(fieldstring,fieldlist,",")}{if($(id)==gene){for(i in fieldlist){printf $fieldlist[i]" "};printf$(id)" \n";}}' | tr ' ' '\t' > $tempgeneout
    results=$( wc -l $tempgeneout | tr -s ' ' | cut -d ' ' -f1 )
    if [[ $results -gt 0 ]]
    then
        cat $tempgeneout >> $temprefcontrol
    else
        grep $gene $synlibrary | awk -v gene=$gene -v id=1 -v fieldstring=$fieldstring -F '\t' 'BEGIN{split(fieldstring,fieldlist,",")}{if($(id)==gene){for(i in fieldlist){printf $fieldlist[i]" "};printf$(id)" \n";}}' | tr ' ' '\t' > $tempgeneout
        results=$( wc -l $tempgeneout | tr -s ' ' | cut -d ' ' -f1 )
        if [[ $results -gt 0 ]]
        then
            cat $tempgeneout >> $temprefcontrol
        else
            echo "$gene could not be located"
            echo "$gene not located" >> $temperrs
        fi
    fi
    rm $tempgeneout
done

## add manual targets if they exist

if [ -z ${manualtargets+x} ]
then
    #no manual targets exist
    echo ""
else
    echo $manualtargets | tr ':' '\t' | tr '-' '\t' | tr ',' '\t'| tr ';' '\n' >> $tempref
fi

# rounding is disabled by default

tempround=$(mktemp -t tmp.ruXXXXXrefround.tsv)
temproundcontrol=$(mktemp -t tmp.ruXXXXXrefbuffontrol.tsv)
tempbuffer=$(mktemp -t tmp.ruXXXXXrefround.tsv)
tempbuffercontrol=$(mktemp -t tmp.ruXXXXXrefbuffcontrol.tsv)
tempall=$(mktemp -t tmp.ruXXXXXallgenes.tsv)
tempunmodified=$(mktemp -t tmp.ruXXXXXallgenesunmodified.tsv)

if [[ $STRANDED -eq 1 ]]
then
    if [[ $round -eq 1 ]]
    then
        echo "rounding to nearest 50kb"
        cat $tempref | awk -v buffer=$buffer -v res=$resolution '{start=int(($2-buffer)/res)*res;stop=int(($3+buffer)/res)*res;if(start<1){start=1};print $1,start,stop,$4,$5}' > $tempround
        cat $temprefcontrol | awk -v buffer=$buffer -v res=$resolution '{start=int(($2-buffer)/res)*res;stop=int(($3+buffer)/res)*res;if(start<1){start=1};print $1,start,stop,$4,$5}' > $temproundcontrol

        cat $tempround $temproundcontrol| sort -k1,1 -k2,2n -k3,3n | uniq > $tempall
    else
        cat $tempref | awk -v buffer=$buffer '{start=$2-buffer;stop=$3+buffer;if(start<1){start=1};print $1,start,stop,$4,$5}' > $tempbuffer
        cat $temprefcontrol | awk -v buffer=$buffer '{start=$2-buffer;stop=$3+buffer;if(start<1){start=1};print $1,start,stop,$4,$5}' > $tempbuffercontrol
        cat $tempbuffer $tempbuffercontrol | sort -k1,1 -k2,2n -k3,3n | uniq > $tempall
    fi
else
    if [[ $round -eq 1 ]]
    then
        echo "rounding to nearest 50kb"
        cat $tempref | awk -v buffer=$buffer -v res=$resolution '{start=int(($2-buffer)/res)*res;stop=int(($3+buffer)/res)*res;if(start<1){start=1};print $1,start,stop,$4}' > $tempround
        cat $temprefcontrol | awk -v buffer=$buffer -v res=$resolution '{start=int(($2-buffer)/res)*res;stop=int(($3+buffer)/res)*res;if(start<1){start=1};print $1,start,stop,$4}' > $temproundcontrol

        cat $tempround $temproundcontrol| sort -k1,1 -k2,2n -k3,3n | uniq > $tempall
    else
        cat $tempref | awk -v buffer=$buffer '{start=$2-buffer;stop=$3+buffer;if(start<1){start=1};print $1,start,stop,$4}' > $tempbuffer
        cat $temprefcontrol | awk -v buffer=$buffer '{start=$2-buffer;stop=$3+buffer;if(start<1){start=1};print $1,start,stop,$4}' > $tempbuffercontrol
        cat $tempbuffer $tempbuffercontrol | sort -k1,1 -k2,2n -k3,3n | uniq > $tempall
    fi
fi



cat $tempref $temprefcontrol | sort -k1,1 -k2,2n -k3,3n | uniq > $tempunmodified

rm $tempref $temprefcontrol $tempround $temproundcontrol $tempbuffer $tempbuffercontrol

# guarantee that input is sorted by chromosome, then start position, then end position

tempfinal=$(mktemp -t tmp.ruXXXXXfinal.tsv)

fixOverlapRecur () {
    filename=$1
    index=$2
    output=$3

    filelines=( $(wc -l $filename | tr -s ' ' | cut -d ' ' -f1 ))
    let looplimit=$filelines-1
    let i=$index
    while [ $i -lt $looplimit ]
    do
        let j=$i+1
        let k=$j+1
        thisline=( $( sed -n "$j"'p' $filename) )
        nextline=( $( sed -n "$k"'p' $filename) )

        #only find overlap if chromosome is the same
        if [[ ${nextline[0]} == ${thisline[0]} ]]
        then
            if [[ ${nextline[1]} -le ${thisline[2]} ]]
            then
                startpos=${thisline[1]}
                if [[ ${nextline[3]} == ${thisline[3]} ]]
                then
                    genename=${thisline[3]}
                else
                    genename="${thisline[3]}.${nextline[3]}"
                fi
                chrname=${thisline[0]}
                if [[ ${thisline[2]} -gt ${nextline[2]} ]]
                then
                    endpos=${thisline[2]}
                else
                    endpos=${nextline[2]}
                fi
                #cat $filename | sed "$j","$k"d > $temp.for.$index.tsv
                tempfori=$(mktemp -t tmp.ruXXXXXfor"$index".tsv)
                awk -v line=$j 'NR<line{print $0}NR==line{exit}' $filename > $tempfori 
                echo -e "$chrname\t$startpos\t$endpos\t$genename" >> $tempfori
                awk -v line=$k 'NR>line{print $0}' $filename >> $tempfori
                fixOverlapRecur $tempfori $i $output
                return
            else
                ((i++))
                fixOverlapRecur $filename $i $output
                return
            fi
        else
            ((i++))
            fixOverlapRecur $filename $i $output
            return
        fi
    done
    
    cat $filename | tr ' ' '\t' > $output
    return
}

fixOverlapRecurStrandAware () {
    filename=$1
    index=$2
    output=$3

    filelines=( $(wc -l $filename | tr -s ' ' | cut -d ' ' -f1 ))
    let looplimit=$filelines-1
    let i=$index
    while [ $i -lt $looplimit ]
    do
        let j=$i+1
        let k=$j+1
        thisline=( $( sed -n "$j"'p' $filename) )
        nextline=( $( sed -n "$k"'p' $filename) )

        #only find overlap if chromosome is the same
        if [[ ${nextline[0]} == ${thisline[0]} ]]
        then
            if [[ ${nextline[1]} -le ${thisline[2]} ]]
            then
                startpos=${thisline[1]}
                if [[ ${nextline[4]} == ${thisline[4]} ]]
                then
                    genename=${thisline[4]}
                else
                    genename="${thisline[4]}.${nextline[4]}"
                fi
                chrname=${thisline[0]}
                strand=${thisline[3]}
                if [[ ${thisline[2]} -gt ${nextline[2]} ]]
                then
                    endpos=${thisline[2]}
                else
                    endpos=${nextline[2]}
                fi
                #cat $filename | sed "$j","$k"d > $temp.for.$index.tsv
                tempfori=$(mktemp -t tmp.ruXXXXXfor"$index".tsv)
                awk -v line=$j 'NR<line{print $0}NR==line{exit}' $filename > $tempfori 
                echo -e "$chrname\t$startpos\t$endpos\t$strand\t$genename" >> $tempfori
                awk -v line=$k 'NR>line{print $0}' $filename >> $tempfori
                fixOverlapRecurStrandAware $tempfori $i $output
                return
            else
                ((i++))
                fixOverlapRecurStrandAware $filename $i $output
                return
            fi
        else
            ((i++))
            fixOverlapRecurStrandAware $filename $i $output
            return
        fi
    done
    
    cat $filename | tr ' ' '\t' > $output
    return
}


namedoutput=$NAME.named.targets.bed
sequenceroutput=$NAME.targets.bed
nakednamedoutput=$NAME.naked.named.targets.bed

cat $tempunmodified | tr ' ' '\t' > $nakednamedoutput

totb=$(awk '{tot+=($3-$2)}END{print tot}' $nakednamedoutput)
totbM=$(echo "scale=2; $totb / 1000000" | bc )
totbC=$(echo "scale=2; $totb / 31000000" | bc)

echo -e "\ntotal unbuffered target size (Mb): $totbM"
echo "total unbuffered target genome percentage: $totbC %"

if [[ $naked -eq 1 ]]
then
    cat $tempall | tr ' ' '\t' > $namedoutput
    awk '{print $1,$2,$3}' $tempall | tr ' ' '\t' > $sequenceroutput
else
    if [[ $STRANDED -eq 0 ]]
    then
        fixOverlapRecur $tempall 0 $tempfinal
        cp $tempfinal $namedoutput
        awk '{print $1,$2,$3}' $tempfinal | tr ' ' '\t' > $sequenceroutput
    else
        tempposstrand=$(mktemp -t tmp.ruXXXXXpos.tsv)
        tempnegstrand=$(mktemp -t tmp.ruXXXXXneg.tsv)
        tempposresults=$(mktemp -t tmp.ruXXXXXposfinal.tsv)
        tempnegresults=$(mktemp -t tmp.ruXXXXXnegfinal.tsv)
        cat $tempall | awk -v posout=$tempposstrand -v negout=$tempnegstrand '{if($5==1){print $0>$posout}else{print $0>negout}}'
        fixOverlapRecurStrandAware $tempposstrand 0 $tempposresults
        fixOverlapRecurStrandAware $tempnegstrand 0 $tempnegresults
        cat $tempposresults $tempnegresults | sort -k1,1 -k2,2n -k3,3n > $tempfinal
        cp $tempfinal $namedoutput
        awk '{if($4==-1){strand="-"}else{strand="+"};print $1,$2,$3,$5,"0",strand}' $namedoutput | tr ' ' '\t' > $sequenceroutput
    fi

    totbuffered=$(awk '{tot+=($3-$2)}END{print tot}' $tempfinal)
    totbufferedM=$(echo "scale=2; $totbuffered / 1000000" | bc )
    totbufferedC=$(echo "scale=2; $totbuffered / 31000000" | bc)

    echo "total final target size (Mb): $totbufferedM"
    echo "total final target genome percentage: $totbufferedC %"  
fi

if [ -f $temperrs ]
then
    cp $temperrs $NAME.errors.txt
fi

echo ""

if [[ $savenetwork -eq 1 ]]
then
    cp $namedoutput /n/alignments/bed_targets/
    cp $nakednamedoutput /n/alignments/bed_targets/
fi

rm /tmp/tmp.ru*tsv
exit 0