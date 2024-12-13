######### PTPN22

auto_hg38=read.table('PTPN22/finngen_R6_AUTOIMMUNE',header = F)
auto_hg37=read.table('PTPN22/auto_hg37.bed',header = F)
auto_hg38$V2=auto_hg37$V2
auto_hg38=auto_hg38[which(auto_hg38$V2>=113541836&auto_hg38$V2<=115044030),]
basal_hg38=read.table('PTPN22/basal_cell',header = F)
basal_hg37=read.table('PTPN22/basal_hg37.bed',header = F)
basal_hg38$V2=basal_hg37$V2
basal_hg38=basal_hg38[which(basal_hg38$V2>=113541836&basal_hg38$V2<=115044030),]
ihs=read.table('PTPN22/ihs',header = F)
ihs=ihs[which(ihs$V4>=113541836&ihs$V4<=115044030),]

pos=numeric()
val=numeric()
while (nrow(ihs)>0) {
  a=ihs[1, 4]+500
  tmp=ihs[which(ihs$V4<a),]
  pos=c(pos, round((tmp[1,4]+tmp[nrow(tmp),4])/2))
  val=c(val, mean(abs(tmp$V5)))
  ihs=ihs[-which(ihs$V4<a),]
}

ihs1=as.data.frame(cbind(pos,val))
ihs=read.table('PTPN22/ihs',header = F)
ihs=ihs[which(ihs$V4>=113541836&ihs$V4<=115044030),]

a=ggplot(basal_hg38, aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#CAC9E7',0.7),size=1.5)+
  labs(x='',y='-log10(P)')+
  annotate("rect", xmin=114000000, xmax=114500000, ymin=0, ymax=Inf, alpha=0.1, fill="#4155A7")+
  ggtitle('Basal cell carcinoma')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))



b=ggplot(auto_hg38, aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#CAC9E7',0.7),size=1.5)+
  labs(x='',y='-log10(P)')+
  annotate("rect", xmin=114000000, xmax=114500000, ymin=0, ymax=Inf, alpha=0.1, fill="#4155A7")+
  ggtitle('Autoimmune diseases')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))



c=ggplot()+ geom_point(data=ihs, aes(x=V4, y=abs(V5)),color=alpha('#CAC9E7',0.7),size=1.5)+
  labs(x='',y='|iHS|')+
  geom_line(data=ihs1, aes(x=pos, y=abs(val)), size=0.6, color=alpha('#4155A7'))+
  annotate("rect", xmin=114000000, xmax=114500000, ymin=0, ymax=Inf, alpha=0.1, fill="#4155A7")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_text(face = 'bold', size = 12),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))


p1=ggarrange(a, b, c, heights = c(2, 2),
             ncol = 1, nrow = 3, align = "v")


p2=annotate_figure(p1, top = text_grob("PTPN22", 
                                       color = "black", face = "bold", size = 14))



####### CD44
ihs=read.table('CD44/ihs',header = F)
ihs=ihs[which(ihs$V1>=34500000&ihs$V1<=35500000), ]
h=read.table('CD44/heel_bone_density',header = F)
lei=read.table('CD44/finngen_R6_CD2_BENIGN_LEIOMYOMA_UTERI',header = F)
h_hg37=read.table('CD44/h_hg37.bed')
lei_hg37=read.table('CD44/lei_hg37.bed')
h$V2=h_hg37$V2
lei$V2=lei_hg37$V2
h=h[which(h$V2>=34500000&h$V2<=35500000),]
lei=lei[which(lei$V2>=34500000&lei$V2<=35500000),]
pos=numeric()
val=numeric()
while (nrow(ihs)>0) {
  a=ihs[1, 1]+500
  tmp=ihs[which(ihs$V1<a),]
  pos=c(pos, round((tmp[1,1]+tmp[nrow(tmp),1])/2))
  val=c(val, mean(abs(tmp$V2)))
  ihs=ihs[-which(ihs$V1<a),]
}

ihs1=as.data.frame(cbind(pos,val))
ihs=read.table('CD44/ihs',header = F)
ihs=ihs[which(ihs$V1>=34500000&ihs$V1<=35500000), ]


a=ggplot(h, aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#CAC9E7',0.7),size=1.5)+
  labs(x='',y='-log10(P)')+
  annotate("rect", xmin=35050000, xmax=35150000, ymin=0, ymax=Inf, alpha=0.1, fill="#4155A7")+
  ggtitle('Heel bone density')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))



b=ggplot(lei, aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#CAC9E7',0.7),size=1.5)+
  labs(x='',y='-log10(P)')+
  annotate("rect", xmin=35050000, xmax=35150000, ymin=0, ymax=Inf, alpha=0.1, fill="#4155A7")+
  ggtitle('Leiomyoma of uterus')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))


c=ggplot()+ geom_point(data=ihs, aes(x=V1, y=abs(V2)),color=alpha('#CAC9E7',0.7),size=1.5)+
  labs(x='',y='|iHS|')+geom_line(data=ihs1, aes(x=pos, y=abs(val)), size=0.6, color=alpha('#4155A7'))+
  annotate("rect", xmin=35050000, xmax=35150000, ymin=0, ymax=Inf, alpha=0.1, fill="#4155A7")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_text(face = 'bold', size = 12),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))


p5=ggarrange(a, b, c, heights = c(2, 2),
             ncol = 1, nrow = 3, align = "v")
p6=annotate_figure(p5, top = text_grob("CD44", 
                                       color = "black", face = "bold", size = 14))




####### NFKB1
hay=read.table('NFKB1/hayfever',header = F)
pbc=read.table('NFKB1/pbc',header = F)
hay_hg37=read.table('NFKB1/hayfever_hg37.bed',header = F)
pbc_hg37=read.table('NFKB1/pbc_hg37.bed',header = F)
hay$V2=hay_hg37$V2
pbc$V2=pbc_hg37$V2
t=read.table('NFKB1/tajima',header = F)
t1=read.table('NFKB1/tajima_10kb',header = F)
beta2=read.table('NFKB1/NFKB1',header = F)


a=ggplot(hay,aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#BCD1E2',0.7),size=1.5) +
  labs(x='',y='-log10(P)')+
  ggtitle('Hayfever')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))+
  annotate("rect", xmin=103420000, xmax=103654000, ymin=0, ymax=Inf, alpha=0.1, fill="#3285B5")


b=ggplot(pbc,aes(x=V2, y=-log10(V3)))+ geom_point(color=alpha('#BCD1E2',0.7),size=1.5) +
  labs(x='',y='-log10(P)')+
  ggtitle('Primary biliary cirrhosis')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))+
  annotate("rect", xmin=103420000, xmax=103654000, ymin=0, ymax=Inf, alpha=0.1, fill="#3285B5")


c=ggplot()+ geom_point(data=t, aes(x=V2, y=V4),color=alpha('#BCD1E2',0.7),size=1.5)+
  labs(x='',y='Tajima\'s D')+
  geom_line(data=t1, aes(x=V2, y=V4), size=0.6, color='#3285B5')+
  annotate("rect", xmin=103500000, xmax=103600000, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3285B5")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))

d=ggplot(beta2,aes(x=V1, y=V2))+ geom_point(color=alpha('#BCD1E2',0.7),size=1.5) +
  labs(x='position',y='Beta2')+
  annotate("rect", xmin=103500000, xmax=103600000, ymin=-2, ymax=Inf, alpha=0.1, fill="#3285B5")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.x=element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))

f2=ggarrange(a, b, c, d, heights = c(2, 2),
             ncol = 1, nrow = 4, align = "v")

f2=annotate_figure(f2, top = text_grob("NFKB1", 
                                       color = "black", face = "bold", size = 14))

### CTLA4
basal=read.table('CTLA4/basal_cell',header = F)
basal_hg37=read.table('CTLA4/basal_cell_hg37.bed',header = F)
hypo=read.table('CTLA4/finngen_R6_E4_HYTHY_AI_STRICT',header=F)
hypo_hg37=read.table('CTLA4/hypo_hg37.bed',header = F)
basal$V2=basal_hg37$V2
basal=basal[which(basal$V2>=204500000&basal$V2<=205000000),]
hypo$V2=hypo_hg37$V2
hypo=hypo[which(hypo$V2>=204500000&hypo$V2<=205000000),]
td=read.table('CTLA4/tajima',header = F)
td=td[-which(is.nan(td$V4)),]
td=td[which(td$V2>=204500000&td$V2<=205000000),]
beta2=read.table('CTLA4/beta2',header = F)
beta2=beta2[which(beta2$V1>=204500000&beta2$V1<=205000000),]
t1=read.table('CTLA4/tajima_10kb',header=F)




a=ggplot(basal, aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#BCD1E2',0.7),size=1.5)+
  labs(x='',y='-log10(P)')+
  annotate("rect", xmin=204680000, xmax=204800000, ymin=0, ymax=Inf, alpha=0.1, fill="#3285B5")+
  ggtitle('Basal cell carcinoma')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))



b=ggplot(hypo, aes(x=V2, y=-log10(V3))) + geom_point(color=alpha('#BCD1E2',0.7),size=1.5)+
  labs(x='',y='-log10(P)')+
  annotate("rect", xmin=204680000, xmax=204800000, ymin=0, ymax=Inf, alpha=0.1, fill="#3285B5")+
  ggtitle('Hypothyroidism')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))




c=ggplot()+geom_point(data=td, aes(x=V2, y=V4),color=alpha('#BCD1E2',0.7),size=1.5)+
  labs(x='',y='Tajima\'s D')+
  geom_line(data=t1, aes(x=V2, y=V4), size=0.6, color='#3285B5')+
  annotate("rect", xmin=204690000, xmax=204710000, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3285B5")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))


d=ggplot()+ geom_point(data=beta2, aes(x=V1, y=V2),color=alpha('#BCD1E2',0.7),size=1.5)+
  geom_raster(stat = 'identity')+
  labs(x='',y='Beta2')+
  annotate("rect", xmin=204690000, xmax=204710000, ymin=0, ymax=Inf, alpha=0.1, fill="#3285B5")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = 'lightgrey', fill = NA, size = 0.5),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = 'bold', size = 12),
        axis.text.x=element_text(face = 'bold', size = 12),
        axis.text.y=element_text(face = 'bold', size = 12))

f1=ggarrange(a, b, c, d, heights = c(2, 2),
             ncol = 1, nrow = 4, align = "v")
f1=annotate_figure(f1, top = text_grob("CTLA4", 
                                       color = "black", face = "bold", size = 14))



ggarrange(p2,p6,f2,f1, ncol = 2, nrow = 2, labels = c('a','b','c','d'))



