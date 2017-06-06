This are several scripts that are used to run pathnet package on expression data.
KEGGREST_annotation.r - use this script to define genes of interest and fetch their entrez id's if they are in different id.
requres a list of genes of interest and loaded database for a species. Uses keggrest.

KEGG_fetching_KGML_files.r - fetch the KGML files that describe the nodes topology. uses files from KEGGREST_annotation.r

pathnet_data_preparation_resolved_pairwise_map.r - using KGML files and annotated genes creaates a pair graph data of gene connectivity.

pathnet_data_preparation.r - using the pair graph data from previous example prepares other files pathnet will need to run.

pathnet_run.r - runs pathnet using the supplementary data.