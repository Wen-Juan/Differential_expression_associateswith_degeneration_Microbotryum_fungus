BEGIN{OFS="\t" ; print "gene"OFS"sequenceLength"OFS"GC0"OFS"GC3"} !/^$/ {
	if(length($2)/3 !~ /[0-9]\.[0-9]/){
		SEQSIZE[$1]=length($2)

		for(i=1;i<=length($2);i+=3){
			CODON=substr($2,i,3) ;

			for(j=1;j<=3;j+=1){

				NUCL=substr(CODON,j,1) ;

				if( NUCL ~ /(G|C)/ ) { GC0[$1] += 1 } ;
				if( j == 3 && NUCL ~ /(G|C)/ ){ GC3[$1] += 1  } ;

			}

		}

	}
}END{
	for(j in GC0){
		print j, SEQSIZE[j], GC0[j]/SEQSIZE[j], GC3[j]/(SEQSIZE[j]/3) ;
	}
}
