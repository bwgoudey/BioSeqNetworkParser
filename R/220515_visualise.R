

library(data.table)
library(igraph)
library(tidyverse)

constrain_edges_to_nodes <- function(edges, n) {
  edges[edges$to %in% n & edges$from %in% n,]
}

# xxx -------
classify_acc <- function(acc) {
  
  #Regexs are derived from a range of sources
  # https://community.gep.wustl.edu/repository/course_materials_WU/annotation/Genbank_Accessions.pdf
  # https://www.uniprot.org/help/accession_numbers
  # PIR are cooked up. 
  # https://registry.identifiers.org/registry/insdc
  # https://registry.identifiers.org/registry/pdb
  id_regexs=list(
       uniprot_protein="^(([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]))(.\\d+)?$",
       genbank_protein="^[A-Z]{3}(\\d{5}|\\d{7})(\\.\\d+)?$",
       genbank_nucleotide="^(([U]\\d{5})|([A-Z]{2}\\d{6}))(\\.\\d+)?$",
       genbank_contig="^([A-Z]{4}\\d{8}(\\d+)?)|([A-Z]{6}\\d{9}(\\d+)?)$",
       refseq_protein="^((AP|NP|XP|YP|ZP)_\\d+)(\\.\\d+)?$",
       refseq.anon_protein="^((WP)_\\d+)(.\\d+)?$",
       refseq_nucleotide="^(((AC|NC|NG|NM|NR|NT|NW|XM|XR)_\\d+)|(NZ\\_[A-Z]{2,4}\\d+))(\\.\\d+)?$",
       pdb_protein="^[0-9][A-Za-z0-9]{3}(_([ABCLX]?))?$",
       pir_protein="^(([A-Z]{2}\\d{4})|([A-OR-TV-Z]{1}\\d{5}))$"
       )
  
  #Setup a  dataframe containing the return values
  # Base it directly on names of regex bank to ensure no mapping issues
  vals=data.frame(r=names(id_regexs)) %>% separate(r, sep="_", into=c("db", "seq_type"))
  df=data.frame(db=rep(NA,length(acc)), seq_type=rep(NA,length(acc)))
  #Match the regex
  I = which(sapply(id_regexs, grepl, acc), arr.ind=T)
  
  #Check to make sure we don't have more than one match
  if(length(I))
    if(is.integer(I) & length(I)==1) # If only a single acc
      df[1,]=vals[I,]
    else {
      stopifnot(all(!duplicated(I[,'row'])))
      df[I[,'row'],] = vals[I[,'col'], ]
    }
  return(df)
}
# xxx  ------

nfile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_node.csv"
nodes<-fread(nfile, col.names = c("id", "seq_ver", "seq_type", "product", "org", "db", "seq"))
nfile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_append_node.csv"
nodes<-rbind(nodes, fread(nfile, col.names = c("id", "seq_ver", "seq_type", "product", "org", "db", "seq")))
#nodes$id=gsub("\\.[0-9]$", "", nodes$id)
setkey(nodes, 'id')
#Reclassify proteins baed on accessions
nodes[, c('db', 'seq_type')]=classify_acc(nodes$id)

efile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_edge.csv"
edges<-fread(efile, col.names = c("trg", "src", "trg_type", "src_type"))
efile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_append_edge.csv"
edges<-rbind(edges, fread(efile, col.names = c("trg", "src", "trg_type", "src_type")))

#edges$trg = gsub("\\.[0-9]$", "", edges$trg)
#edges$src = gsub("\\.[0-9]$", "", edges$src)

# Reclassify edges based on accessions ----
edges=edges[, c('trg', 'src')] %>% mutate(trg_db=classify_acc(trg)$db, trg_seq_type=classify_acc(trg)$seq_type,src_db=classify_acc(src)$db, src_seq_type=classify_acc(src)$seq_type )

#Remove edges that are not from relevant sources of proteins
is_relevant <- \(seq_type,db) seq_type=="protein" & db %in% c("genbank", "refseq", "regseq.anon", "uniprot")
edges_prot<-edges[is_relevant(edges$trg_seq_type, edges$trg_db) & is_relevant(edges$src_seq_type, edges$src_db), ]

#  xxx -------
g<-graph_from_edgelist(as.matrix(edges_prot)[,c(2,1)], directed = TRUE)
V(g)$db<-nodes[names(V(g)), 'db'][[1]]
V(g)$organism<-nodes[names(V(g)), 'org'][[1]]
V(g)$product<-stringr::str_to_title(nodes[names(V(g)), 'product'][[1]])
V(g)$seq<-stringr::str_sub(cli::hash_md5(nodes[names(V(g)), 'seq'][[1]]), 1, 3)

#Find all connected components that have any sort of path between them
comp=igraph::components(g, mode = "weak")

# Number of sequences
ggplot(nodes, aes(x=seq_type)) + geom_bar() + scale_y_log10()

# Number of connections
ggplot(data.frame(csize=c(comp$csize,unique(comp$csize))), aes(x=csize)) + geom_histogram(bins=100) + scale_y_log10()

# What are the linked records ----------
table(classify_acc(unique(c(edges$rsc, edges$trg))), useNA = 'ifany')

# sdfsdf --------
#Extract the largest connected component
lrg_cc=names(comp$membership[comp$membership %in% tail(order(comp$csize), n = 3)])
lrg_cc_g<-subgraph(g, lrg_cc)

library(RColorBrewer)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

lapply(c("db", "seq", "organism", "product"), function(x) {
  u=unique(vertex_attr(lrg_cc_g, x))
  ggnet2(lrg_cc_g,
       label = TRUE, label.size = 5,
       arrow.size = 12, arrow.gap = 0.025, 
       color=x, palette = setNames(c25[1:length(u)],u))
})


