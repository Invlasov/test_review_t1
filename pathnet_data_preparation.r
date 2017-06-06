setwd("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files")
library(XML)
library(PathNet)
# load the test dataset

current <- getwd()
setwd(system.file(dir="extdata", package="PathNetData"))
# Begin loading datasets from the text files
brain_regions <- as.matrix(read.table(
file = "brain_regions_data.txt", sep = "\t", header = T))
disease_progression <- as.matrix(read.table(
file = "disease_progression_data.txt", sep = "\t", header = T))
A <- as.matrix(read.table(
file = "adjacency_data.txt", sep = "\t", header = T))
pathway <- read.table(
file = "pathway_data.txt", sep = "\t", header = T)
# Change back to our previous working directory
setwd(current)

filepath<-c("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files/pairwise_table.txt")
pairwise_table<-read.table("pairwise_table.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
# dirpath<-c("D:/Microarray_files/mice_steatosis/scan/ALL_fc_files/fc_annotated_kegg")
filepath<-c("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files/all_annotated_genes_for_kgml/all_common_entrez.txt_kegg_annotated.txt")
pw_filtration_table<-read.table(filepath, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

# first thing that need to be done is reformation of data according to our microarray platform.
# As such, we need to 
# 1)remove pathways from pairwise matrix
# 2)convert names to entrez IDs
# 3)remove all the edges that include missing genes

# pointless idea - there is better way to do that
# remove_path<-unique(c(grep("path:",pairwise_table[,1]),grep("path:",pairwise_table[,2])))
# keep<-setdiff(1:nrow(pairwise_table),remove_path)
# pairwise_table<-pairwise_table[keep,]
keggnames<-unique(pw_filtration_table$KEGG)
convert<-pw_filtration_table[,1]
names(convert)<-pw_filtration_table$KEGG

pairwise_table_v2<-pairwise_table
for (i in 1:nrow(pairwise_table_v2)) {
	if (any(pairwise_table_v2[i,1]==keggnames)) {
	pairwise_table_v2[i,1]<-convert[pairwise_table_v2[i,1]]
	} else {
	pairwise_table_v2[i,1]<-0
	}
	if (any(pairwise_table_v2[i,2]==keggnames)) {
	pairwise_table_v2[i,2]<-convert[pairwise_table_v2[i,2]]
	} else {
	pairwise_table_v2[i,2]<-0
	}
}
merged<-vector()
for (i in 1:nrow(pairwise_table_v2)) {
merged<-c(merged,paste(pairwise_table_v2[i,1],pairwise_table_v2[i,2],pairwise_table_v2[i,3],collapse=NULL,sep="_"))
}
pairwise_table_v2<-pairwise_table_v2[!duplicated(merged),]


# pairwise_table<-pairwise_table[firstrow&secondrow,]

firstrow<-vector()
secondrow<-vector()
for (i in 1:nrow(pairwise_table)) {
	if (any(pairwise_table[i,1]==keggnames)) {
	firstrow<-c(firstrow,TRUE)
	pairwise_table[i,1]<-convert[pairwise_table[i,1]]
	} else {
	firstrow<-c(firstrow,FALSE)
	}
	if (any(pairwise_table[i,2]==keggnames)) {
	secondrow<-c(secondrow,TRUE)
	pairwise_table[i,2]<-convert[pairwise_table[i,2]]
	} else {
	secondrow<-c(secondrow,FALSE)
	}
}
pairwise_table<-pairwise_table[firstrow&secondrow,]


A_side<-as.character(sort(as.numeric(unique(c(pairwise_table[,1],pairwise_table[,2])))))
A_custom<-matrix(rep(0,length(A_side)^2),ncol=length(A_side),nrow=length(A_side))
colnames(A_custom)<-A_side
rownames(A_custom)<-A_side

for (i in 1:nrow(pairwise_table)) {
left<-pairwise_table[i,1]
right<-pairwise_table[i,2]
left<-match(left,A_side)
right<-match(right,A_side)
A_custom[left,right]<-1
}

write.table(pairwise_table, file="pairwise_table.txt", quote = FALSE, sep = "\t",row.names=FALSE,col.names=TRUE)
write.table(pairwise_table_v2, file="pairwise_table_modified_v2.txt", quote = FALSE, sep = "\t",row.names=FALSE,col.names=TRUE)
write.table(A_custom, file="adjacency_matrix.txt", quote = FALSE, sep = "\t",row.names=TRUE,col.names=TRUE)

