

library(data.table)
library(igraph)
library(tidyverse)
source('~/Research/BioDbNetQual/BioDbPropagationEval/R/classify_acc.R', echo=TRUE)

constrain_edges_to_nodes <- function(edges, n) {
  edges[edges$to %in% n & edges$from %in% n,]
}

# xxx  ------

nfile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_node.csv"
nodes<-fread(nfile)#, col.names = c("id", "seq_ver", "seq_type", "product", "org", "db", "seq"))
nfile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_append_node.csv"
nodes<-rbind(nodes, fread(nfile))#, col.names = c("id", "seq_ver", "seq_type", "product", "org", "db", "seq")))
#nodes$id=gsub("\\.[0-9]$", "", nodes$id)
setkey(nodes, 'id')
#Reclassify proteins baed on accessions
nodes[, c('db', 'seq_type')]=classify_acc(nodes$id)

#Fix dates
nodes$date_first_upload=lubridate::ymd(nodes$date_first_upload)
nodes$date_last_modified=lubridate::ymd(nodes$date_last_modified)


efile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_edge.csv"
edges<-fread(efile, col.names = c("child", "parent"))#, "parent_type", "child_type"))
efile="/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/tmp_EC_3_4_11_4_partial_sequence_append_edge.csv"
edges<-rbind(edges, fread(efile, col.names =c("child", "parent"))) #c("trg", "src", "parent_type", "child_type")))

#edges$trg = gsub("\\.[0-9]$", "", edges$trg)
#edges$src = gsub("\\.[0-9]$", "", edges$src)
# Reclassify -----
# Reclassify edges based on accessions ----
edges=edges[, c('child', 'parent')] %>% mutate(parent_db=classify_acc(parent)$db, 
                                               parent_seq_type=classify_acc(parent)$seq_type,
                                               child_db=classify_acc(child)$db, 
                                               child_seq_type=classify_acc(child)$seq_type )

#Remove edges that are not from relevant sources of proteins
is_relevant <- \(seq_type,db) seq_type=="protein" & db %in% c("genbank", "refseq", "regseq.anon", "uniprot")
edges_prot<-edges[is_relevant(edges$parent_seq_type, edges$parent_db) & 
                    is_relevant(edges$child_seq_type, edges$child_db), ]


# Descriptive plots ---------------
table(nodes$date_first_upload, useNA = 'ifany')
ggplot(nodes, aes(x=date_first_upload, fill=db)) + geom_density(alpha=0.4)

#Date last modified
ggplot(nodes, aes(x=date_last_modified, fill=db)) + geom_density(alpha=0.4, adjust=2)

#Updates to the sequence
ggplot(nodes, aes(x=seq_version, fill=db)) + geom_bar(position='dodge') + scale_y_log10()


t(t(sort(table(sapply(nodes$taxonomy, \(x) stringr::str_split(x, ",")[[1]][6], USE.NAMES = F)), dec=T)[1:10]))

# Descriptive of the network -----------

all_edges=rbind(edges[,c('child', 'child_db','child_seq_type')], edges[,c('parent', 'parent_db','parent_seq_type')], use.names=F)
all_edges<-all_edges[!duplicated(all_edges$child)]
setkey(all_edges, 'child')

#  xxx -------
g<-graph_from_edgelist(as.matrix(edges_prot)[,c(2,1)], directed = TRUE)
V(g)$db<-all_edges[names(V(g)), 'child_db'][[1]]
V(g)$seq_type<-all_edges[names(V(g)), 'child_seq_type'][[1]]
V(g)$organism<-nodes[names(V(g)), 'organism'][[1]]
V(g)$product<-stringr::str_to_title(nodes[names(V(g)), 'name'][[1]])
V(g)$seq<-stringr::str_sub(cli::hash_md5(nodes[names(V(g)), 'seq'][[1]]), 1, 3)
V(g)$ec<-nodes[names(V(g)), 'ec'][[1]]
V(g)$taxa_id<-nodes[names(V(g)), 'taxa_id'][[1]]
#V(g)$seq_type<-nodes[names(V(g)), 'seq_type'][[1]]
V(g)$date_last_modified<-nodes[names(V(g)), 'date_last_modified'][[1]]
V(g)$date_first_upload<-nodes[names(V(g)), 'date_first_upload'][[1]]
V(g)$date_first_upload<-nodes[names(V(g)), 'date_first_upload'][[1]]

#Find all connected components that have any sort of path between them
comp=igraph::components(g, mode = "weak")



# xxx -------
# Number of sequences
#ggplot(nodes, aes(x=seq_type, fill=seq_type)) + geom_bar() + scale_y_log10()

# Number of connections
ggplot(data.frame(csize=comp$csize), aes(x=csize)) + 
  geom_histogram(bins=100) +
  scale_x_log10() + 
  scale_y_log10()+ 
  labs(title=sprintf("Size of connected components"), x="Connected component size") 


edges %>% 
  group_by(child, child_db) %>% 
  summarise(degree=n()) %>% 
  ungroup() %>%
  ggplot(aes(x=degree)) + geom_histogram()  + scale_y_log10() + 
  labs(title=sprintf("Node degree (one at 113 not shown\n due to log scale)"))

edges %>% 
  group_by(child, child_db) %>% 
  summarise(degree=n()) %>% 
  ungroup() %>%
  ggplot(aes(x=degree, fill=child_db)) + geom_histogram()  + scale_y_log10() + scale_x_log10() + 
  labs(title=sprintf("Node degree (one at 113 not shown\n due to log scale)"))


edges_u=rbind(edges[,c('child', 'child_db','child_seq_type')], edges[,c('parent', 'parent_db','parent_seq_type')], use.names=F)
edges_u %>% 
  group_by(child, child_db) %>% 
  summarise(degree=n()) %>% 
  ungroup() %>%
  ggplot(aes(x=degree, fill=child_db)) + geom_density(alpha=0.5) + scale_x_log10() + scale_y_sqrt(limits=c(0,6)) +
  labs(title=sprintf("Node degree density (without Genbank)")) 

edges_u %>% 
  group_by(child, child_db) %>% 
  summarise(degree=n()) %>% 
  ungroup() %>%
  ggplot(aes(x=degree, fill=child_db)) + geom_histogram(bins=5, position='dodge') + scale_x_log10() + scale_y_log10() 
  labs(title=sprintf("Node degree density (with Genbank)")) 




# What are the linked records ----------
uniq_acc=unique(c(edges$child, edges$parent))
ca=classify_acc(uniq_acc)
table(ca, useNA = 'ifany')

# sdfsdf --------
#Extract the largest connected component
lrg_cc=names(comp$membership[comp$membership %in% tail(order(comp$csize), n = 2)[1]])
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
V(lrg_cc_g)$db_st<-str_c(V(lrg_cc_g)$db, '_', V(lrg_cc_g)$seq_type)
lapply(c("db_st", "db", "seq_type")[1], function(x) {#, "organism", "product"), function(x) {
  u=unique(vertex_attr(lrg_cc_g, x))
  ggnet2(lrg_cc_g, alpha=0.5,
       label = TRUE, label.size = 5,
       arrow.size = 12, arrow.gap = 0.025, 
       color=x, palette = setNames(c25[1:length(u)],u))
})


# How many edges on average (per DB, over time) ----





