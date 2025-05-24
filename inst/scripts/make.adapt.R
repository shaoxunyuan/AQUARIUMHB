
args = commandArgs(trailingOnly = TRUE)
sampleID = args[1]
options(stringsAsFactors=F)
options(warn=-1)


library(Biostrings)
library(data.table)
sampletable <- fread("sampletable.txt",data.table = F)

AvgSpotLen <- sampletable[sampletable$Run == sampleID,]$AvgSpotLen
adaptlen <- AvgSpotLen/2-1
# ¶ÁÈ¡fastaÎÄ¼þ

input_fasta_file <- paste0(sampleID,"/quant/circRNA_raw.fa")
output_fasta_file <- paste0(sampleID,"/quant/circRNA_final.fa")
file.copy(input_fasta_file,output_fasta_file,overwrite = TRUE)
output_fasta_seq <- readDNAStringSet(output_fasta_file)
for(index in 1:length(output_fasta_seq)){
	seq = output_fasta_seq[index];
	seq_for_adaptors <- DNAStringSet(c(rep(unlist(seq),ceiling(adaptlen/width(seq)))))
	adaptors <- subseq(seq_for_adaptors,start=width(seq_for_adaptors)-adaptlen+1,end=width(seq_for_adaptors));
	new_seqs <- paste0(adaptors, seq)
	output_fasta_seq[index] <- new_seqs
}
writeXStringSet(output_fasta_seq, output_fasta_file, format = "fasta")
