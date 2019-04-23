### This script is to use Kallisto mapping algorithm to systematically quantify a series of RNAseq libraries/samples.

module add UHTS/Analysis/kallisto/0.44.0

for f in *_pairedR1.fq.gz
do
kallisto quant -i ~/Genome/transcriptome.idx -o ./output/${f%%_pairedR1.fq.gz} -b 1000 ${f%%_pairedR1.fq.gz}_pairedR1.fq.gz ${f%%_pairedR1.fq.gz}_pairedR2.fq.gz
done
