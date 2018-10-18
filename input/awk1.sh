awk '
NR==FNR{ range[$1,$2,$3]; next }
FNR==3
{
    for(x in range) {
        split(x, check, SUBSEP); 
        if($1==check[1] && $2>=check[2] && $3<=check[3]) print $0
    }
}    
' BOX_NRR.txt A1_gene_location_fi.txt
