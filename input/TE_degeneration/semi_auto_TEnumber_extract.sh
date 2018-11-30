cut -f2 Mvsl_DEnonDE_A2_20kdownstream__TE_number.txt > test.txt

wait

grep -Fv -f test.txt  Mvsl_DEnonDE_A1A2_genelist.txt | cut -f1 > test_no.txt

wait 

awk '{$1 = $1 OFS (NR==1?"TE":"0")} 1' test_no.txt | tr ' ' \\t > test_no_no.txt

wait

cut -f1 test_no_no.txt > test_no1.txt

wait 

cut -f2 test_no_no.txt > test_no2.txt

wait

paste test_no2.txt test_no1.txt > test_no21.txt

wait

sed '1d' test_no21.txt > test_no212.txt

wait

cat Mvsl_DEnonDE_A2_20kdownstream__TE_number.txt test_no212.txt > Mvsl_DEnonDE_A2_20kdownstream_allgenes_TEs.txt
