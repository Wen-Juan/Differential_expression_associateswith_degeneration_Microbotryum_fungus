##These are script to run prank for codon alignment.

module add SequenceAnalysis/MultipleSequenceAlignment/prank/170427

for i in *.fa;
do
        prank -d=$i -o=$i -f=phylips -F -codon

done
