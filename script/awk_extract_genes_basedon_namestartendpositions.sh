#####This awk script can extract gene lists, based on the gene name, start and end position (column 1, 2 and 3).
awk 'NR==FNR{ range[$1,$2,$3]; next }
FNR==1
{
    for(x in range) {
        split(x, check, SUBSEP);
        if($1==check[1] && $2>=check[2] && $3<=check[3]) print $0
    }
}
' BOX_centro_A2.txt A2_gene_location_fi.txt
