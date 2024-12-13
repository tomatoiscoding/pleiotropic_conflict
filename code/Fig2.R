library(ggplot2)
library(stringr)
library(xlsx)
### fig 2a
fig2=read.table('all_trait',header = F)
t=character()
for (i in 1:nrow(fig2)) {
  t=c(t, unlist(str_split(fig2$V1[i],',')))
}
for (i in 1:nrow(fig2)) {
  t=c(t, unlist(str_split(fig2$V2[i],',')))
}
t=unique(t)


cl=xlsx::read.xlsx('class219_immune.xlsx',header = T, sheetIndex = 1)


fig2_matrix=matrix(NA, nrow = nrow(fig2), ncol = nrow(cl))

for (i in 1:nrow(fig2)) {
  t1=unlist(str_split(fig2$V1[i],','))
  a=match(t1,cl$phenocode)
  fig2_matrix[i, a]=1
  t2=unlist(str_split(fig2$V2[i],','))
  b=match(t2,cl$phenocode)
  fig2_matrix[i, b]=2
}

cl1=matrix(NA,length(unique(cl$category))*length(unique(cl$category)),2)
for (i in 1:length(unique(cl$category))) {
  cl1[(length(unique(cl$category))*(i-1)+1):(length(unique(cl$category))*i),1]<-unique(cl$category)[i]
}
for (i in 1:length(unique(cl$category))) {
  cl1[(length(unique(cl$category))*(i-1)+1):(length(unique(cl$category))*i),2]<-unique(cl$category)
}
cl1=as.data.frame(cl1)
cl1$count=0
for (i in 1:nrow(fig2_matrix)) {
  a=unique(cl[which(fig2_matrix[i, ]==1),'category'])
  b=unique(cl[which(fig2_matrix[i, ]==2),'category'])
  for (j in 1:length(a)) {
    for (k in 1:length(b)) {
      cl1[which(cl1$V1==a[j]&cl1$V2==b[k]),'count']=cl1[which(cl1$V1==a[j]&cl1$V2==b[k]),'count']+1
    }
  }
}

for (i in unique(cl$category)) {
  a=which(cl1[, 'V2']==i)
  b=which(cl1[, 'V1']==i)
  c=intersect(a,b)
  a=a[-which(a<=c)]
  b=b[-which(b<=c)]
  cl1[a, 'count']=cl1[a, 'count']+cl1[b, 'count']
  cl1[b, 'count']=cl1[a, 'count']
}

cl1$radius <- sqrt( cl1$count / pi )
colnames(cl1)=c('T1','T2','count','radius')
cl1$T2=as.factor(cl1$T2)
l1=levels(cl1$T2)
a=numeric()
for (i in 1:length(l1)) {
  a=c(a,sum(cl1[which(cl1$T1==l1[i]),'count']))
}

### plot
colfunc <- colorRampPalette(c("#9EC8B9", "#092635"))
co=colfunc(29)[as.numeric(cut(a,breaks = 29))]
l2=as.character(a)
d1=ggplot(cl1,aes(T1,as.numeric(T2)))+
  geom_point(aes(size=ifelse(count==0,NA,radius*7.5)),shape=21,fill=rgb(0, 145, 167, maxColorValue = 255), color=rgb(0, 145, 167, maxColorValue = 255), alpha=0.5)+
  geom_text(aes(label=ifelse(count==0, '', count)),size=4, color='white')+
  scale_size_identity()+
  geom_abline(intercept =0 , slope = 1, colour=rgb(0, 145, 167, maxColorValue = 255), linetype='dotted')+
  theme(
        axis.text.x=element_text(angle=90,hjust=1,vjust=0, face = 'bold', size = 10),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_text(face = 'bold', size = 10),
        axis.text.y.right = element_text(color = co),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = '#DCF5E3', size = 0.5, linetype = 'dotted'),
        plot.background = element_blank())+
  scale_y_continuous(breaks = 1:length(l1),
                     labels = l1,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(l2),
                                         labels = l2))


#### fig2b

library(stringr)
library(tidyr)
library(readxl)

T2=read_xlsx('TableS2.xlsx', skip = 1)
T6=read_xlsx('TableS6.xlsx', skip = 1)

setwd('/Users/uqbbian/pleiotropy/all/clumping/R6')
files=list.files()
files=files[-1]
t1=t
for (f in files) {
  if(startsWith(f, 'finngen')) {
    t=read.table(f,header=T)
    t=t%>%separate(independent_SNP,c('chr','pos','ref','alt'))
    t$chr=str_replace(t$chr,'chr','')
    t$v=paste0(t$chr,':',t$pos)
  } else {
    t <- read.table(f, header = F)
    t$v= paste0(t$V1,':',t$V2)
  }
  t$dis_not=NA
  t$trait_info = f
  snp = character()
  for (i in 1:nrow(T2)) {
    traits <- str_split(T2$trait1[i], ',')[[1]]
    if (f %in% traits) {
      snp <- c(snp, T2$variant[i])
    }
  }
  for (i in 1:nrow(T2)) {
    traits <- str_split(T2$trait2[i], ',')[[1]]
    if (f %in% traits) {
      snp <- c(snp, T2$variant[i])
    }
  }
  if (length(snp) > 0) {
    for (i in 1:length(snp)) {
      a=match(snp[i], t$v)
      if(!is.na(a)) {
        t$dis_not[a] <- 'discordant'
      } else {
        chr <- str_split(snp[i],':')[[1]][1]
        ld <- read.table(paste0('LD/bi_LD_chr',chr,'.tags.list'), header = T)
        ld <- ld%>% separate(SNP, c('chr','pos','ref','alt'))
        ld$v <- paste0(ld$chr, ':', ld$pos)
        candidate <- ld[match(snp[i], ld$v), ]
        tags <- str_split(candidate$TAGS, '\\|')[[1]]
        tags <- as.data.frame(tags)
        tags <- tags%>%separate(tags, c('chr','pos','ref','alt'))
        tags$v <- paste0(tags$chr,':',tags$pos)
        b=intersect(tags$v, t$v)
        if(sum(is.na(b))>0) {
          t$dis_not=t$dis_not
        } else {
          t$dis_not[match(intersect(tags$v, t$v), t$v)] <- 'discordant'
        }
      }
    }
  }
  
  snp = character()
  
  for (i in 1:nrow(T6)) {
    traits <- str_split(T6$trait[i], ',')[[1]]
    if (f %in% traits) {
      snp <- c(snp, T6$variant[i])
    }
  }
  if (length(snp)>0) {
    for (i in 1:length(snp)) {
      a=match(snp[i], t$v)
      if(!is.na(a)) {
        t$dis_not[a] <- 'concordant'
      } else {
        chr <- str_split(snp[i],':')[[1]][1]
        ld <- read.table(paste0('con_LD/con_LD/con_LD_chr',chr,'.tags.list'), header = T)
        ld <- ld%>% separate(SNP, c('chr','pos','ref','alt'))
        ld$v <- paste0(ld$chr, ':', ld$pos)
        candidate <- ld[match(snp[i], ld$v), ]
        tags <- str_split(candidate$TAGS, '\\|')[[1]]
        tags <- as.data.frame(tags)
        tags <- tags%>%separate(tags, c('chr','pos','ref','alt'))
        tags$v <- paste0(tags$chr,':',tags$pos)
        b=intersect(tags$v, t$v)
        if(sum(is.na(b))>0) {
          t$dis_not=t$dis_not
        } else {
          t$dis_not[match(intersect(tags$v, t$v), t$v)] <- 'concordant'
        }
      }
    }
  }
  
  t$dis_not[which(is.na(t$dis_not))] <- 'single_domain'
  t=t[,c('v','dis_not','trait_info')]
  t1=rbind(t1,t)
}

num_domain_bi <- as.data.frame(table(T2$`number of domain`))
num_domain_con <- as.data.frame(table(T6$`number of domain`))
num_domain_con$Var1 <- as.numeric(num_domain_con$Var1)
num_domain_con <- rbind(num_domain_con, list(1, 16904))
num_domain_con$Var1[1:5]=2:6
num_domain_con$Var1 = as.factor(num_domain_con$Var1)

a1=ggplot(num_domain_bi, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity", fill='#bdb5e1', alpha=0.6)+
  labs(x = 'Number of domain', y = 'Number of variants')+
  geom_text(aes(label = Freq), vjust = -0.3, size = 3, fontface = 'bold')+
  ggtitle('Antagonistic variants associated with one or more domains')+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text = element_text(face="bold", size = 12),
        axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold", angle = 90),
        plot.title = element_text(face = 'bold', hjust = 0, size = 12))


b1=ggplot(num_domain_con, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity", fill='#bdb5e1', alpha=0.6)+
  scale_y_break(c(1000, 15000))+
  scale_y_continuous(breaks = c(0, 1000, 16000))+
  scale_x_discrete(limits = factor(1:8))+
  geom_text(aes(label = Freq), vjust = -0.3, size = 3, fontface = 'bold')+
  labs(x = 'Number of domain', y = 'Number of variants')+
  ggtitle('Other genetic variants associated with one or more domains')+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.text = element_text(face="bold", size = 12),
        axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold", angle = 90),
        plot.title = element_text(face = 'bold', hjust = 0, size = 12))

c1=ggarrange(a1, print(b1), align = 'h', ncol = 1)
ggarrange(d1, c1, ncol = 2, nrow = 1, labels = c('a', 'b'))
