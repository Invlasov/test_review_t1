# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("mouse.db0")
source("D:/Ivan/rscripts/exported_functions.r")
library(KEGGREST)
# Sys.setenv("http_proxy" = "http://proxy.mh-hannover.de:8080")
dirpath<-c("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/concordant_files")
library(org.Mm.eg.db)

# directory with genes of interest
# it is implied that first row are ID's, sep is \t
setwd("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files/all_annotated_genes_for_kgml")
i<-1
filevector<-list.files(path=dirpath)
for (l in 1:length(filevector)) {
	filepath<-paste(dirpath,filevector[l],sep="/",collapse=NULL)
	data_df<-read.table(filepath,header=TRUE, sep="\t", stringsAsFactors=FALSE,quote="")
	genenames<-data_df[,1]
	forconv<-paste("ncbi-geneid:",genenames,sep="",collapse=NULL)
	
	outputvec<-vector()
	for (i in 1:length(forconv)) {
	kconv<-keggConv("hsa", forconv[i])
		if (!identical(kconv, character(0))) {
		outputvec<-c(outputvec, keggConv("hsa", forconv[i]))
		}
		else {
		outputvec<-c(outputvec, "")
		}
	}
	data_df$KEGG<-unname(outputvec)
	# res_mf<-formatter(res_mf)
	outname<-paste(filevector[l],"kegg_annotated.txt",sep="_",collapse=NULL)
	write.table(data_df, file=outname, quote = FALSE, sep = "\t",row.names=FALSE)
	# }

}

setwd("D:/Microarray_files/mice_steatosis/scan/functional_enrichment/keggrest")
