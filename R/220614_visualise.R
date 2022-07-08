library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
source('R/classify_acc.R')

coalesce_join <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}


cols<-c("query", "target", "seq_id", "length", "mismatch", "gaps", "q_start", "q_end" ,"t_start", "t_end", "e_value", "bit_score")
# 
#ALL prots -------
all_vs_sp<-data.table::fread('./data/gb_all_220706.res', col.names = cols)
all_vs_sp_split=all_vs_sp %>% 
  group_by(query) %>% 
  slice(which.max(seq_id))  %>% 
  tidyr::separate(query, sep="_", into=c('q_acc', 'q_ec', 'q_type')) %>% 
  tidyr::separate(target, sep="_", into=c('t_acc', 't_ec', 't_type'))
p1<-all_vs_sp_split %>%  ggplot(aes(x=seq_id)) + geom_histogram()+ xlim(0,100) + theme_bw(base_size=14) + scale_y_log10(limits=c(1,1e6))
p1



poor_lab_dat=all_vs_sp_split %>% 
  filter(seq_id<35)#&seq_id<35) %>% 
  #slice(which.max(bit_score))
write.table(file='./notebook/exps/poor_ec_sim.acc', poor_lab_dat %>% select(q_acc), row.names = F, col.names = FALSE, quote = F)


#all_vs_sp %>% group_by(target) %>% summarise(m=max(seq_id)) %>% filter(m<30)  %>% dim()


#all_vs_sp %>% group_by(target) %>% summarise(m=max(bit_score)) %>% ggplot(aes(x=m)) + geom_histogram() + theme_bw() + scale_x_log10(limits = c(10,50))+ scale_y_log10() 

#all_vs_sp %>% group_by(target) %>% summarise(m=max(-log10(e_value))) %>%   ggplot(aes(x=m)) + geom_histogram() + theme_bw() + scale_x_log10()+ scale_y_log10()
#SP prots ------
sp_vs_sp<-read.table('~/Research2/BioDbNetQual/220411_ExamineFunction/data/EC3_4_11_4/diamond/SP_vs_SP.tsv', col.names = cols, row.names = NULL)

sp_vs_sp %>% group_by(query) %>% 
  filter(query!=target) %>%
  summarise(m=max(seq_id)) %>% 
  ggplot(aes(x=m)) + geom_density() + geom_histogram()+ xlim(0,100)
sp_vs_sp %>% group_by(query) %>% summarise(m=max(seq_id)) %>% filter(m<30)  %>% dim()
# Reannot ----------
new_lab_dat<-data.table::fread('./data/poor_ec_sim.fasta2_res', col.names = cols) %>%
  group_by(query) %>%
  slice(which.max(seq_id )) %>% 
  ungroup() %>% 
  #tidyr::separate(query, sep="\\.", into=c('q_acc', 'q_ver')) %>% 
  tidyr::separate(query, sep="_", into=c('q_acc', 'q_ec', 'q_type')) %>% 
  tidyr::separate(target, sep="_", into=c('t_acc', 't_ec', 't_type')) %>% 
  left_join(poor_lab_dat %>% select(q_acc, t_acc, t_ec, seq_id, q_ec, bit_score, e_value),
            by="q_acc", suffix=c(".new", ".old")) 
  #filter(!str_detect(q_ec, t_ec) & !str_detect(t_ec, q_ec) & !str_detect(q_ec, '-') & !str_detect(t_ec, '-'))


#xx<-new_lab_dat %>% 
#  filter(bit_score.new-bit_score.old > 0 & 
#           seq_id.new-seq_id.old > 0) %>% 
#  select(q_acc, t_acc.old, t_acc.new, t_ec.old,t_ec.new,  bit_score.new, bit_score.old, seq_id.new, seq_id.old) %>% 
#  arrange(bit_score.new-bit_score.old )
I=new_lab_dat$bit_score.new-new_lab_dat$bit_score.old <= 10 
xx=new_lab_dat
xx$t_ec.new[I]=xx$t_ec.old[I]
xx$t_acc.new[I]=xx$t_acc.old[I]
xx$seq_id.new[I]=xx$seq_id.old[I]
xx$bit_score.new[I]=xx$bit_score.old[I]
xx$e_value.new[I]=xx$e_value.old[I]



p2<-#inner_join(poor_lab_dat,new_lab_dat, by=c('t_acc', 't_ec', 't_type'), suffix=c(".old", ".new")) %>% 
  xx %>% ungroup() %>% 
  ggplot(aes(x=seq_id.old, y=seq_id.new)) +
  #ggplot(aes(x=bit_score.old, y=bit_score.new-bit_score.old)) +
  #ggplot(aes(x=-log10(e_value.old), y=-log10(e_value.new))) +
  geom_point(alpha=0.5) + 
  theme_bw(base_size=14) +
  labs(x="identity to old label", y="improvement in identity")+ 
  coord_cartesian(xlim=c(10,50),ylim=c(0,100))
p2

x=xx %>% ungroup() %>% 
  mutate(orig_seq_id=cut(seq_id.old, breaks=c(0,27,31, 35)),
         Status=factor(t_ec.new==t_ec.old, c(TRUE, FALSE), c('Same', 'New'))) 
#x$Status=pbsapply(1:nrow(x), function(i) length(intersect(str_split(x[i,]$t_ec.old, ",")[[1]], str_split(x[i,]$t_ec.new, ",")[[1]]))>0)
#x$Status=factor(x$Status, c(TRUE, FALSE), c('Same', 'New'))



p2<-#inner_join(poor_lab_dat,new_lab_dat, by=c('t_acc', 't_ec', 't_type'), suffix=c(".old", ".new")) %>% 
  x %>% ggplot(aes(x=seq_id.new, fill=Status, color=Status)) + 
  #ggplot(aes(x=bit_score.old, y=bit_score.new-bit_score.old)) +
  #ggplot(aes(x=-log10(e_value.old), y=-log10(e_value.new))) +
  geom_histogram(position='dodge') + 
  #geom_density(alpha=0.5) + 
  theme_bw(base_size=14) +
  labs(x="identity to new label", y="count") +# 
 scale_y_log10(limits=c(1,1e6))
  #coord_cartesian(xlim=c(10,50),ylim=c(0,100))
p2

p<-cowplot::plot_grid(plotlist = list(
  p1 + labs(x="Original Sequence Identity", y="# of records") + 
    scale_y_continuous(limits=c(0,1e6)) + 
    geom_vline(xintercept=35, linetype="dashed",color='red'), 
  p2 + labs(x="New Sequence Identity", y="#of records") + 
    scale_y_continuous(limits=c(0,1e6))+ 
    scale_x_continuous(limits=c(0,100))
  ))
ggsave(filename = './notebook/reannot.pdf', plot=p, width = 24, height=8, units="cm")
# ------


#Load in all the data
nodes_all<-data.table::fread('./data/gbbct10.seq_node.csv', 
                      header=T)#, 
                      #col.names=c("trg", "trg_ver", "src","src_ver", "x", "db", "type","x1"), sep = "\t") %>%
edges_all<-read.table('./data/gbbct10.seq_annot_edge.csv', 
                      header=F, 
                      col.names=c("trg", "trg_ver", "src","src_ver", "x", "db", "type","x1"), sep = "\t") %>%
              select(-x1)

edges_all_inference<-read.csv('./data/gbbct10.seq_annot_edge.accs', header=F, col.names='acc')
edges_dead_found<-read.csv('./data/gbbct10.seq_annot_edge.accs_not_removed_but_dead.txt', header=F, col.names='acc')
edges_dead<-read.csv('./data/gbbct10.seq_annot_edge.accs_removed.txt', header=F, col.names='acc')
edges_reannot<-read.csv('./data/gbbct10.seq_annot_edge.accs_reannot.txt', header=F, col.names='acc')
dead_acc=rbind(edges_dead_found, edges_dead)

# Create a list of edges with inference and indicate their status
# Create a list of edges and their 



xx <- nodes_all %>% 
  mutate(inference="none") %>%
  left_join(edges_all %>% select(trg, src), by=c("id"="trg")) %>% 
  mutate(inference=ifelse(id %in% edges_all$trg, "alive", inference)) %>%
  mutate(status=ifelse(id %in% dead_acc$acc, "dead", inference)) %>%
  mutate(status=ifelse(id %in% edges_reannot$acc, "reannot", inference))




edges_all$status='alive'
edges_all_inference$status[edges_all_inference$acc %in% rbind(edges_dead_found, edges_dead)$acc]='dead'
edges_all_inference$status[edges_all_inference$acc %in% edges_reannot$acc]='reannot'
edges_all_inference %>% count(status) %>% 
  mutate(perc = 100* n / nrow(edges_all_inference)) %>% 
  ggplot(aes(x = status, y = perc, fill=status)) + geom_bar(stat = "identity") + theme_bw(base_size = 14) +
  labs(y="Records (%)")

  

