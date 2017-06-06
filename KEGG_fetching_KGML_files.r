source("D:/Ivan/rscripts/exported_functions.r")
library(KEGGREST)
# Sys.setenv("http_proxy" = "http://proxy.mh-hannover.de:8080")
annonation_file_path<-c("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files/all_annotated_genes_for_kgml/all_common_entrez.txt_kegg_annotated.txt")
annotation_df<-read.table(annonation_file_path,header=TRUE, sep="\t", stringsAsFactors=FALSE,quote="")
annotation_df[,1]<-as.character(annotation_df[,1])
library(org.Mm.eg.db)
setwd("D:/Microarray_files/TVX_project/Functional_enrcihment/KGML")

annotation_df<-annotation_df[annotation_df$KEGG!="",]
keggenes<-unique(annotation_df$KEGG)

genes_to_pathways<-keggLink("pathway", "hsa")
	
pws<-genes_to_pathways[keggenes]
allpws<-vector()
for (i in 1:length(keggenes)) {
	if (!is.na(genes_to_pathways[keggenes[i]])) {
	allpws<-c(allpws,genes_to_pathways[keggenes[i]])
	}
}
allpws<-unique(unname(allpws))

for (i in seq(1,length(allpws),by=1)) {
# batchy<-allpws[i:min((i+9),length(allpws))]
kegget<-keggGet(allpws[i], option =  "kgml")
outname<-gsub(":","_",paste(allpws[i],".xml",collapse=NULL,sep=""))
write(kegget,file=outname,sep = "\t")
}