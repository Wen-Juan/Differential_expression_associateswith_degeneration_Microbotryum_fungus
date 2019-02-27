## requirement bed tools
BIN='/home/hirak/bedtools2/bin'
## Gencode 
## gencode.v29.chr_patch_hapl_scaff.annotation.gtf
GTF_FILE="gencode.v29.chr_patch_hapl_scaff.annotation.gtf"

# extract transcript boundaries 
cat $GTF_FILE | awk 'BEGIN{OFS="\t";} $3=="transcript" {print $1,$4-1,$5,$12}' | tr -d "\"" | tr -d ";" | $BIN/sortBed > gencode_transcript_intervals.bed

# merge exon boundaris 
cat $GTF_FILE | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5,$12}' | tr -d "\"" | tr -d ";" | $BIN/sortBed | $BIN/mergeBed -i - -c 4 -o collapse > gencode_exon_merged.bed

# extract introns per transcript
$BIN/subtractBed -a gencode_transcript_intervals.bed -b gencode_exon_merged.bed -nonamecheck > intron_transcript.bed
