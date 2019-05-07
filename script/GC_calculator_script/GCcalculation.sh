#!/bin/bash

while getopts a: option
do
    case "${option}"
    in
        a) fasta=${OPTARG};;
    esac
done

output=$(basename ${fasta} | sed 's/\..*//g')

awk '{if($1 ~ /^>/){print $1}else{print toupper($0)}}' ${fasta} | \
        sed 's/cds://g' | \
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | \
	paste -s -d"#" | sed 's/#/ /g;s/>/\n/g' | \
        awk -f GCcalculation.awk | tr '#' '\t' | sort -k1,1 -k2,2n > ${output}_basecomp.txt


