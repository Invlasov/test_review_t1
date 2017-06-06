setwd("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files")
library(XML)
library(PathNet)
source("D:/Ivan/rscripts/exported_functions.r")

full_names_restoration<-function(namesvec,pw_table) {
	replacevec<-vector()
	for (j in 1:length(namesvec)) {
	replacement<-grep(namesvec[j],unique(pairwise_table[,3]),fixed=TRUE,value=TRUE)
		if (length(nchar(replacement))>1) {
		# print(replacement)
		shortest<-(min(nchar(replacement))==nchar(replacement))
		replacement<-replacement[shortest]
		}
	replacevec<-c(replacevec,replacement)
	}
	return(replacevec)
}


pairwise_table<-read.table("pairwise_table_modified_v2.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
# pairwise_table<-read.table("pairwise_table_modified.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

A_custom<-read.table("adjacency_matrix.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="",row.names=1)
A_custom<-data.matrix(A_custom)

root_output_dir<-c("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/output")
dirpath<-c("D:/Microarray_files/mice_steatosis/scan/ALL_fc_files/fc_annotated_entrez_id")
filevector<-list.files(path=dirpath)

selection<-c(1:length(filevector))
# selection<-c(2)
i<-1
for (i in selection) {
filepath<-paste(dirpath,filevector[i],collapse=NULL,sep="/")
evidence<-read.table(filepath, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
threshold_sig<-evidence[de.reg(evidence$Linear.FC,evidence$p.val.padj,linear=TRUE,threshold=1.5),1]
if (length(threshold_sig)<2) {next}
FC<-log(evidence$Linear.FC,base=2)
m<-mean(FC)
s<-sd(FC)
lowtail<-pnorm(log(1/1.5,base=2), mean = m, sd = s, lower.tail = TRUE, log.p = FALSE)
hightail<-(1-pnorm(log(1.5,base=2), mean = m, sd = s, lower.tail = TRUE, log.p = FALSE))
FC_significance_th<-max(lowtail,hightail)

direct_evidence<-vector()
	for (j in 1:nrow(evidence)) {
		if (FC[j]<0) {
		potential<-pnorm(FC[j], mean = m, sd = s, lower.tail = TRUE, log.p = FALSE)
		} else {
		potential<-(1-pnorm(FC[j], mean = m, sd = s, lower.tail = TRUE, log.p = FALSE))
		}
		direct_evidence<-c(direct_evidence,max(potential,((evidence$p.val.padj[j]/0.05)*FC_significance_th)))
	}
evidence$DRE<-direct_evidence
evidence$DRE_sig<-(direct_evidence<FC_significance_th)

evidence$DE_log<-(-log(evidence$DRE,base=10))


# evidence$DE<-(-log(direct_evidence,base=10))
combined_evidence_sig<-evidence[evidence$DRE_sig,1]
fullness<-length(intersect(combined_evidence_sig,threshold_sig))/length(threshold_sig)
overshot<-length(setdiff(combined_evidence_sig,threshold_sig))/length(threshold_sig)
message_1<-paste("Fulness: ",fullness,", Overshot: ",overshot,sep="",collapse=NULL)
print(message_1)

results <- PathNet(Enrichment_Analysis = TRUE,Contextual_Analysis=TRUE,use_sig_pathways=TRUE,
DirectEvidence_info = evidence,
Adjacency = A_custom,
pathway = pairwise_table,
Column_DirectEvidence = ncol(evidence),
n_perm = 5000, threshold = FC_significance_th)

output_enrichment_results<-results$enrichment_results
output_combined_evidence<-results$enrichment_combined_evidence
output_conn_p_value<-results$conn_p_value
output_pathway_overlap<-results$pathway_overlap
# now we need to restore ridiculous short names of the results

colnames(output_conn_p_value)<-full_names_restoration(colnames(output_conn_p_value),pairwise_table)
rownames(output_conn_p_value)<-full_names_restoration(rownames(output_conn_p_value),pairwise_table)

colnames(output_pathway_overlap)<-full_names_restoration(colnames(output_pathway_overlap),pairwise_table)
rownames(output_pathway_overlap)<-full_names_restoration(rownames(output_pathway_overlap),pairwise_table)

output_enrichment_results[,1]<-full_names_restoration(output_enrichment_results[,1],pairwise_table)

new_wd(root_output_dir,"pathnet_enrichment")
outname<-paste(gsub(".txt","",filevector[i]),"pathnet_enrichment.txt",sep="_",collapse=NULL)
write.table(output_enrichment_results,file=outname, quote = FALSE, sep = "\t",row.names=FALSE)

new_wd(root_output_dir,"combined_evidence")
outname<-paste(gsub(".txt","",filevector[i]),"pathnet_combined_evidence.txt",sep="_",collapse=NULL)
write.table(output_combined_evidence,file=outname, quote = FALSE, sep = "\t",row.names=FALSE)

new_wd(root_output_dir,"contextual_association")
outname<-paste(gsub(".txt","",filevector[i]),"pathnet_contextual_association.txt",sep="_",collapse=NULL)
write.table(output_conn_p_value,file=outname, quote = FALSE, sep = "\t",row.names=TRUE)

new_wd(root_output_dir,"pathnet_overlap")
outname<-paste(gsub(".txt","",filevector[i]),"pathnet_pathway_overlap.txt",sep="_",collapse=NULL)
write.table(output_pathway_overlap,file=outname, quote = FALSE, sep = "\t",row.names=TRUE)
}



