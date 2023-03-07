rm(list=ls())
cat("\014") 
graphics.off()
options(scipen = 999)

library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)

setwd(dirname(envDocument::get_scriptpath()))

feature = read.table('Data For Figure 1A.txt',header=T,stringsAsFactors=F,sep='\t',comment.char=';')
literature = read.table('Data For Figure 1B.txt',header=T,stringsAsFactors=F,sep='\t',comment.char=';')
xmax = 598

p1 = ggplot(feature,aes(xmin=start,xlower=start,xmiddle=0,xupper=end,xmax=end,
                        y=string,group=PlotGroup,fill=Color,color=TextColor,label=Label)) +
  geom_segment(data=feature,aes(x=1,xend=xmax-1,y=string,yend=string),color='black',size=2) +
  geom_boxplot(stat='identity',show.legend=F,position='identity',color='black',size=1,width=1) + ylab('Antigens') +
  geom_text(data=feature %>% filter(Antigen!='ORF1ab'),mapping=aes(x=midpoint,y=string),fontface='bold',size=9,show.legend=F) +
  scale_fill_identity(guide = "legend", aes(labels = lab)) +
  scale_color_identity(guide = "legend", aes(labels = lab)) +
  scale_x_continuous(limits=c(1,xmax),breaks = seq(0,1000,100),expand=c(0,0)) +
  theme(text = element_text(size = 24)) + theme_bw() + ylab('BNT162b4') +
  theme(axis.text = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(angle=0,vjust=0.5,size=30,face='bold'),
        axis.ticks = element_blank(), rect = element_blank(),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

p2 =ggplot(literature,aes(x=position,y=fct_rev(Type),fill=count)) +
  geom_tile() + ylab('Literature\nEpitopes') +
  scale_fill_gradientn(colors=c('white','blue','red'),na.value='white',
                       breaks=seq(0,34,length.out=3),limits=c(0,34)) +
  scale_x_continuous(limits=c(1,xmax),breaks = seq(0,1000,100),expand=c(0,0)) +
  scale_y_discrete(labels=c('Reported Epitope\nDec 2022','Reported Epitope\nSept 2020')) +
  theme(text = element_text(size = 24)) + theme_bw() + labs(fill = '# of Occurences') +
  theme(axis.text.x = element_blank(), axis.title = element_blank(),axis.text.y = element_text(size=20,face='bold',color='black'),
        axis.ticks = element_blank(), rect = element_blank(),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(size=0.5,color='gray70',linetype=3), panel.grid.minor.x = element_blank(),
        legend.position='bottom',legend.justification = c(1,0),legend.key.width=unit(40,'pt'),
        legend.text=element_text(size=12,face='bold'),legend.title=element_text(size=16,face='bold'),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),legend.box.spacing = unit(0,'pt'))

windows(height=600,width=2000)
plot(plot_grid(p1,NULL,p2,align='v',nrow=3,ncol=1,axis = 'lr',
               rel_heights = c(1.2,0.5,2)))

