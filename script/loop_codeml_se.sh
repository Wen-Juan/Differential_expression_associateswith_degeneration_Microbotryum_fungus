###This script is to simultaneously run thousands of files of paired fasta sequences to run PAML for calculating dN, dS, or dN/dS.
## module
module add Phylogeny/paml/4.9g

for i in *.fa.best.phy; do

out_ctl_file="codeml.test.txt"
echo $out_ctl_file

results_f=`echo $i | sed 's/.fa.best.phy/.results.txt/g'`

echo "seqfile  =  "$i"         * sequence data filename" >> $out_ctl_file
echo "outfile  =  "$results_f"   * main result file name" >> $out_ctl_file
echo "" >> $out_ctl_file
echo "noisy =  9      * 0,1,2,3,9: how much rubbish on the screen" >> $out_ctl_file
echo "verbose =   1      * 1:detailed output" >> $out_ctl_file
echo "runmode =   -2     * -2:pairwise" >> $out_ctl_file
echo "" >> $out_ctl_file
echo "seqtype =   1      * 1:codons" >> $out_ctl_file
echo "CodonFreq =   2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61" >> $out_ctl_file
echo "model =  0      *" >> $out_ctl_file
echo "NSsites =   0      *" >> $out_ctl_file
echo "icode =  0     * 0:universal code" >> $out_ctl_file
echo "" >> $out_ctl_file
echo "fix_kappa =   0      * 1:kappa fixed, 0:kappa to be estimated" >> $out_ctl_file
echo "kappa =  2      * initial or fixed kappa" >> $out_ctl_file
echo "" >> $out_ctl_file
echo "fix_omega =   0     * 1:omega fixed, 0:omega to be estimated" >> $out_ctl_file
echo "omega =  0.001  * 1st fixed omega value [change this]" >> $out_ctl_file
echo "getSE = 1   *" >> $out_ctl_file

# PAML code here

codeml $out_ctl_file

done
