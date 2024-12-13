### fig3a

bi_af=read.table('fig3/bi_AF.txt',header = T)
con_af=read.table('/Users/uqbbian/pleiotropy/all/V2/fig3/con_AF.txt',header = T)
bi_af[which(bi_af$bi>0.5), 'bi']=1-bi_af[which(bi_af$bi>0.5), 'bi']
con_af[which(con_af$con>0.5), 'con']=1-con_af[which(con_af$con>0.5), 'con']
a=ggplot(con_af,aes(x=con))+geom_histogram(bins= 30,aes(y=after_stat(density)), fill="light blue", alpha=0.5)+geom_density(alpha=.2, fill="#FF6666", colour='#FF6666')+labs(x='MAF')+ggtitle('Concordant pleiotropic variants affect traits in different domains')+theme(axis.text = element_text(face="bold"),axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"))
b=ggplot(bi_af,aes(x=bi))+geom_histogram(bins= 30,aes(y=after_stat(density)), fill="light blue", alpha=0.5)+geom_density(alpha=.2, fill="#FF6666", colour='#FF6666')+labs(x='MAF')+ggtitle('Antagonistic variants')+theme(axis.text = element_text(face="bold"),axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"))


variable1=c(rep('bi',nrow(bi_frq)),rep('con',nrow(con_frq)))
value1=c(bi_frq$V5,con_frq$V5)
data=as.data.frame(cbind(variable1,value1))

colnames(data)=c('variable','value')
a=ggplot(con_af,aes(x=con))+geom_histogram(bins= 30,aes(y=after_stat(density)), fill="#B3D1E7", alpha=0.5)+
  geom_density(alpha=.3, fill="#D9B9D4", colour='#D9B9D4')+labs(x='MAF')+ggtitle('Concordant pleiotropic variants affect traits in different domains')+
  theme(axis.text = element_text(face="bold", size = 12),axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = 'bold'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))


b=ggplot(bi_af,aes(x=bi))+geom_histogram(bins= 30,aes(y=after_stat(density)), fill="#B3D1E7", alpha=0.5)+
  geom_density(alpha=.3, fill="#D9B9D4", colour='#D9B9D4')+labs(x='MAF')+ggtitle('Antagonistic variants')+
  theme(axis.text = element_text(face="bold", size = 12),axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = 'bold'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

bb=ggarrange(a, b, ncol = 2, nrow = 1)

### fig3b

a=ggplot(median_by_group1, aes(x=factor(trait, levels = c('basal_cell', 'breast_cancer',
                                                          'prostate_cancer',
                                                          'finngen_R6_C3_PROSTATE_EXALLC',
                                                          'finngen_R6_CD2_BENIGN_LEIOMYOMA_UTERI',
                                                          'hayfever',
                                                          'finngen_R6_ASTHMA_MODE',
                                                          'crohn',
                                                          'IBD',
                                                          'UC',
                                                          'T1D',
                                                          'finngen_R6_AUTOIMMUNE',
                                                          'finngen_R6_E4_HYTHY_AI_STRICT',
                                                          'finngen_R6_E4_HYTHYNAS',
                                                          'finngen_R6_HYPOTHYROIDISM',
                                                          'T2D',
                                                          'finngen_R6_E4_DM2_STRICT',
                                                          'finngen_R6_E4_ENDONUTRMET',
                                                          'finngen_R6_E4_LIPOPROT',
                                                          'Disorders_of_lipoid_metabolism',
                                                          'High_cholesterol',
                                                          'finngen_R6_E4_THYROID',
                                                          'finngen_R6_E4_THYTOXGOITDIF',
                                                          'finngen_R6_CHOLELITH_BROAD',
                                                          'finngen_R6_FG_DISVEINLYMPH',
                                                          'Hypertension',
                                                          'finngen_R6_FG_HYPERTENSION',
                                                          'finngen_R6_I9_HYPERTENSION',
                                                          'finngen_R6_I9_ANGIO',
                                                          'finngen_R6_I9_CHD',
                                                          'finngen_R6_I9_CORATHER',
                                                          'stroke',
                                                          'SCZ')), y=median_by_group, color=group))+
  geom_point(size=2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold', size = 10, color=c(rep('#D097BA', 5),
                                                                                            rep('#5DBDE8', 2),
                                                                                            rep('#EE6A33', 8),
                                                                                            rep('#89BCE4', 8),
                                                                                            '#874F8D',
                                                                                            rep('#FBCB1F',8),
                                                                                            '#47AF79')),
        axis.text.y = element_text(face = 'bold', size = 10, color='black'),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size=10, face = 'bold'),
        #panel.background = element_blank(),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        plot.margin = unit(c(.5,6,.5,.5),"lines"))+
  labs(x='', y='Median of genetic variance')+
  labs(color='')+
  scale_x_discrete(labels=c('basal_cell'='Basal cell carcinoma', 'breast_cancer'='Breast cancer',
                            'crohn'='Crohn\'s disease', 'Disorders_of_lipoid_metabolism'='Disorders of lipoid metabolism',
                            'hayfever'='Hayfever',
                            'prostate_cancer'='Prostate cancer',
                            'SCZ'='Schizophrenia',
                            'stroke'='Stroke',
                            'UC'='Ulcerative colitis'))+
  annotate("text", x = Inf, y = 0.004, label = 'Malignant/Benign neoplasm', hjust = -0.25, size = 3, color='#D097BA')+
  annotate("text", x = Inf, y = 0.0035, label = 'Respiratory ', hjust = -0.6, size = 3, color='#5DBDE8')+
  annotate("text", x = Inf, y = 0.003, label = 'Autoimmune ', hjust = -0.55, size = 3, color='#EE6A33')+
  annotate("text", x = Inf, y = 0.0025, label = 'Metabolic ', hjust = -0.7, size = 3, color='#89BCE4')+
  annotate("text", x = Inf, y = 0.002, label = 'Digestive ', hjust = -0.7, size = 3, color='#874F8D')+
  annotate("text", x = Inf, y = 0.0015, label = 'Circulatory ', hjust = -0.6, size = 3, color='#FBCB1F')+
  annotate("text", x = Inf, y = 0.001, label = 'Psychiatric ', hjust = -0.6, size = 3, color='#47AF79')


g = ggplotGrob(a)
g$layout$clip[g$layout$name == "panel"] = "off"
grid.draw(g)


### fig3c

permutation=read.table('/Users/uqbbian/pleiotropy/all/V2/check_permu/result',header = F)

cc=ggplot(permutation, aes(x = V1)) + geom_histogram(fill="#D2E9CB", color="#E8F3E1", alpha=0.6)+ggtitle('Permutation Distribution')+labs(x = 'mean Fst', y = 'Frequency')+
  geom_segment(aes(x=0.134086, xend=0.134086, y=0, yend=1000), colour="#F9CCE0", lwd=1.5, linetype='dotted')+
  annotate('text', x=0.134086, y=1030, label='0.134', color='#F9CCE0', fontface=2)+
  geom_segment(aes(x=0.1158, xend=0.1158, y=0, yend=1000), colour="#D2E9CB", lwd=1.5, linetype='dotted')+
  annotate('text', x=0.1158, y=1030, label='0.116', color="#D2E9CB", fontface=2)+
  theme(axis.text = element_text(face="bold", size = 12),
        axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = 'bold'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))



bb=ggplotGrob(bb)
cc=ggplotGrob(cc)
grid.arrange(bb, g, cc, nrow=3, ncol=1)

