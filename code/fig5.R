library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)
library(ggpubr)
test=read.csv('Fig5.csv',header = T)
gene_name=c('TNFRSF14','FSHB','KLHL21','NFKB1','PTPN22','CTLA4','WNT4','MMS22L')
func=c('immune','reproductive','MHC I','immune','immune','immune','female sex development','DNA repair')
p=ggplot(test,aes(SNP,Age))+
  geom_point(aes(size=ifelse(count==0,NA,7)),shape=21,colour='white', fill='#FCD1CC', alpha=0.7)+
  scale_size_identity()+
  theme(panel.grid.major=element_line(linetype=2,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0, face = 'bold', size = 13),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_text(face = 'bold', size = 13))+scale_x_discrete(labels= paste0(unique(test$SNP),'\n',gene_name,'\n',func))

p2=p+theme_tufte()+geom_segment(aes(x = 1 , y = 1, xend = 1, yend = 6), colour='yellow',linetype='dashed')+annotate('text',x=1, y=c(0.7,6.3),label=c('Hay fever','Hypothyroidism'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 2 , y = 3, xend = 2, yend = 6), colour='yellow',linetype='dashed')+annotate('text',x=2, y=c(2.7,5.3,6.3),label=c('Age at first birth','Leiomyoma of uterus','Age at menopause'),colour = c('#1961A2','black','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 3 , y = 3, xend = 3, yend = 6), colour='yellow',linetype='dashed')+annotate('text',x=3, y=c(2.7,6.3),label=c('Schizophrenia','Hypertension'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 4 , y = 1, xend = 4, yend = 6), colour='yellow',linetype='dashed')+annotate('text',x=4, y=c(0.7,6.3),label=c('Hay fever\nDiseases of middle ear and mastoid','Primary biliary cholangitis'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 5 , y = 2, xend = 5, yend = 7), colour='yellow',linetype='dashed')+annotate('text',x=5, y=c(1.7,7.3),label=c('Type 1 diabetes','Basal cell carcinoma'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 6 , y = 2, xend = 6, yend = 7), colour='yellow',linetype='dashed')+annotate('text',x=6, y=c(1.7,7.3),label=c('Type 1 diabetes','Basal cell carcinoma'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 7 , y = 3, xend = 7, yend = 6), colour='yellow',linetype='dashed')+annotate('text',x=7, y=c(2.7,6.3),label=c('Endometriosis','Female genital prolapse'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  geom_segment(aes(x = 8 , y = 3, xend = 8, yend = 7), colour='yellow',linetype='dashed')+annotate('text',x=8, y=c(2.7,7.3),label=c('Age at first birth','Diverticulosis'),colour = c('#1961A2','black'),fontface = "bold", size=6)+
  theme(axis.text = element_text(face="bold", size = 20),axis.title = element_text(face='bold',size = 22))+xlab('')


age=read.table('age.txt',header = F)
age_c=age %>% group_by(V1,V2) %>% tally()
p1=ggplot(age_c, aes(V1, V2))+geom_point(aes(size=n),color=rgb(130,160,193,maxColorValue = 255), alpha=0.7)+
  scale_x_discrete(name ="Age", limits=c("0-10",'10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90'))+
  scale_y_discrete(name ="Age", limits=c("0-10",'10-20','20-30','30-40','40-50','50-60','60-70'))+
  annotate('rect', xmin=2.8,xmax=3.2,ymin=1.8,ymax=3.2,fill = '#D9F5F4',alpha=0.5)+
  annotate('rect',xmin=3.8,xmax=9.2,ymin=0.8,ymax=3.2,fill = '#E9E6F7',alpha=0.5)+
  annotate('rect', xmin=3.8,xmax=9.2,ymin=3.8,ymax=7.2,fill = '#FAF1DE',alpha=0.5)+
  labs(size = 'N')+
  theme(axis.text = element_text(face="bold", size = 20),
        axis.title = element_text(face='bold',size = 22),
        panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = '#DCF5E3', size = 0.5, linetype = 'dotted'),
        plot.background = element_blank())



ggarrange(p1,p2,ncol = 2,nrow = 1,labels = c('a','b'),widths = c(1.3,2))


