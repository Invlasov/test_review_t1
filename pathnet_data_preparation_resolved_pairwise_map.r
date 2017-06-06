setwd("D:/Microarray_files/TVX_project/Functional_enrcihment/pathnet/support_files")
library(XML)
library(PathNet)

current <- getwd()

# load the test dataset
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

dirpath<-c("D:/Microarray_files/TVX_project/Functional_enrcihment/KGML")
filevector<-list.files(path=dirpath)
resolved_relations<-list()
counter<-1
# i<-1
for (i in 1:length(filevector)) {
	filepath<-paste(dirpath,filevector[i],sep="/",collapse=NULL)
	xmlfile <- xmlTreeParse(filepath)
	xml_list <- xmlToList(xmlRoot(xmlfile))
	pwname<-unname(xml_list$.attrs["title"])
	xmllist_nodes<-xml_list[names(xml_list)=="entry"]
	xmllist_edges_relation<-xml_list[names(xml_list)=="relation"]
	if (length(xmllist_edges_relation)==0) {next}
	
	# finally, nodes
	
	formatted_node_list<-list()
	
	for (j in 1:length(xmllist_nodes)) {
	current_node<-xmllist_nodes[[j]]$.attrs
	main<-(names(current_node)!="id" & names(current_node)!="link")
	formatted_node_list[[current_node["id"]]]<-current_node[main]
	
	}
	
	# !reaction
	interest_in_reaction<-FALSE
	# they are typically out of our field of interest. This is low molecular weight compounds
	if (interest_in_reaction){
formatted_reaction_list_id<-list()
formatted_reaction_list_kegg_id<-list()
xmllist_edges_reaction<-xml_list[names(xml_list)=="reaction"]

# I know that this is insane, bear with me
	# this is outer cycle. It goes through all the reactions
	for (j in 1:length(xmllist_edges_reaction)) {
	interactions_in_this_reaction_id<-vector()
	interactions_in_this_reaction_keggid<-vector()

		# this is inner cycle. The structure of every reaction is 1-n inputs + 1-m outputs. every input has an id, and every output has it. Id then kegg id.
		# for the cautiosness sake, i store both
		for (k in seq(1,length(xmllist_edges_reaction[[j]]),by=2)) {
			# now we input all the "normal" id into one list, and all other into another. Name of id in list is the name of reaction. 
			# !not id - they all are just liken ormal numbers, field for errors is enormous here
			if (names(xmllist_edges_reaction[[j]])[k]==".attrs") {
			formatted_reaction_list_id[[xmllist_edges_reaction[[j]][[k]]["name"]]]<-interactions_in_this_reaction_id
			formatted_reaction_list_kegg_id[[xmllist_edges_reaction[[j]][[k]]["name"]]]<-interactions_in_this_reaction_keggid
			} else {
			interactions_in_this_reaction_id<-c(interactions_in_this_reaction_id,xmllist_edges_reaction[[j]][[k]])
			interactions_in_this_reaction_keggid<-c(interactions_in_this_reaction_keggid,xmllist_edges_reaction[[j]][[k+1]])
			}		
		}
	}
	list_of_returns_reactions<-vector()
	for (j in 1:length(formatted_reaction_list_id)) {
	reaction<-formatted_reaction_list_id[[j]]
	for (k in 1:length(reaction)) {
	list_of_returns_reactions<-c(list_of_returns_reactions,formatted_node_list[[reaction[k]]]["type"])
	}
}
}
	
	# relation
	# it is important here to separate gene relations from maplinks.
	formatted_relations<-list()
	for (j in 1:length(xmllist_edges_relation)) {
		if (!is.atomic(xmllist_edges_relation[[j]])) {
		attrs_relation<-xmllist_edges_relation[[j]]$.attrs
		} else {
		attrs_relation<-xmllist_edges_relation[[j]]
		}
		formatted_relations[[j]]<-unlist(attrs_relation[1:2],use.names=FALSE)
		}
			
	
	for (j in 1:length(formatted_relations)) {
	reaction<-formatted_relations[[j]]
	resolved_relations[[counter]]<-c(unname(formatted_node_list[[reaction[1]]]["name"]),unname(formatted_node_list[[reaction[2]]]["name"]),pwname)
	counter<-counter+1
	}
 # now to resolve all that insanity
}
resolved<-as.data.frame(matrix(unlist(resolved_relations), ncol=3,byrow=TRUE),stringsAsFactors=FALSE)

resolved_single_1<-vector()
resolved_single_2<-vector()
resolved_single_3<-vector()
splat<-strsplit(resolved[,1]," ")

for (i in 1:nrow(resolved)) {
	if (length(splat[[i]])==1) {
	resolved_single_1<-c(resolved_single_1,resolved[i,1])
	resolved_single_2<-c(resolved_single_2,resolved[i,2])
	resolved_single_3<-c(resolved_single_3,resolved[i,3])
	} else {
		for (j in 1:length(splat[[i]])) {		
		resolved_single_1<-c(resolved_single_1,splat[[i]][j])
		resolved_single_2<-c(resolved_single_2,resolved[i,2])
		resolved_single_3<-c(resolved_single_3,resolved[i,3])
		}
	}
}


resolved<-data.frame(id1=resolved_single_1,id2=resolved_single_2,title=resolved_single_3,stringsAsFactors=FALSE)

resolved_single_1<-vector()
resolved_single_2<-vector()
resolved_single_3<-vector()

splat<-strsplit(resolved[,2]," ")
for (i in 1:nrow(resolved)) {
	if (length(splat[[i]])==1) {
	resolved_single_1<-c(resolved_single_1,resolved[i,1])
	resolved_single_2<-c(resolved_single_2,resolved[i,2])
	resolved_single_3<-c(resolved_single_3,resolved[i,3])
	} else {
		for (j in 1:length(splat[[i]])) {		
		resolved_single_2<-c(resolved_single_2,splat[[i]][j])
		resolved_single_1<-c(resolved_single_1,resolved[i,1])
		resolved_single_3<-c(resolved_single_3,resolved[i,3])
		}
	}
}
pairwise_with_maps<-data.frame(id1=resolved_single_1,id2=resolved_single_2,title=resolved_single_3,stringsAsFactors=FALSE)
write.table(pairwise_with_maps,"pairwise_table.txt",sep="\t",row.names=FALSE,quote=FALSE)