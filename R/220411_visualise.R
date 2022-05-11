library(dplyr)
library(ggplot)
cols<-c("query", "target", "seq_id", "length", "mismatch", "gaps", "q_start", "q_end" ,"t_start", "t_end", "e_value", "bit_score")

#ALL prots
all_vs_sp<-read.table('~/Research2/BioDbNetQual/220411_ExamineFunction/data/EC3_4_11_4/diamond/ALL_vs_SP.tsv', col.names = cols)

all_vs_sp %>% group_by(query) %>% summarise(m=max(seq_id)) %>% ggplot(aes(x=m)) + geom_histogram()+ xlim(0,100)
all_vs_sp %>% group_by(query) %>% summarise(m=max(seq_id)) %>% filter(m<30)  %>% dim()



#SP prots
sp_vs_sp<-read.table('~/Research2/BioDbNetQual/220411_ExamineFunction/data/EC3_4_11_4/diamond/SP_vs_SP.tsv', col.names = cols, row.names = NULL)

sp_vs_sp %>% group_by(query) %>% 
  filter(query!=target) %>%
  summarise(m=max(seq_id)) %>% 
  ggplot(aes(x=m)) + geom_density() + geom_histogram()+ xlim(0,100)
sp_vs_sp %>% group_by(query) %>% summarise(m=max(seq_id)) %>% filter(m<30)  %>% dim()

