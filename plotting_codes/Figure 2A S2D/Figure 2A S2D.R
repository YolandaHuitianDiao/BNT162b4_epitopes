rm(list=ls())
cat("\014") 
graphics.off()
options(scipen = 999)

BNTForest = '#3C9673'
BNTJungle = '#005A64'
BNTLime = '#BECD32'
BNTBlue = '#002060'
  
library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)

vgrepl <- Vectorize(grepl, vectorize.args = "pattern")

setwd(dirname(envDocument::get_scriptpath()))

feature = read.table('BNT162b4 Mapping for Figure 2A+S2D.txt',header=T,stringsAsFactors=F,sep='\t',comment.char=';')
hit = read.table('BNT162b4 Epitopes for Figure 2A+S2D.txt',header=T,stringsAsFactors=F,sep='\t',comment.char=';')
hydro_df = read.table('BNT162b4 Hydropathy for Figure S2D.txt',header=T,stringsAsFactors=F,sep='\t',comment.char=';') %>%
  mutate(string='')

xmax=598
x1 = 0.01
x2 = 0
rowmin = min(hit$row)
hit$crow = factor(hit$row,levels=sort(unique(hit$row),decreasing=T))
levels(hit$crow) = c(levels(hit$crow),as.character(rowmin-1:2))

p1 = ggplot(feature,aes(xmin=start,xlower=start,xmiddle=0,xupper=end,xmax=end,
                        y=string,group=PlotGroup,fill=Color,color=TextColor,label=Label)) +
  # plot line background for length of template
  geom_segment(data=feature,aes(x=1,xend=xmax-1,y=string,yend=string),color='black',size=2) +
  geom_segment(data=feature,aes(x=1,xend=xmax-1,y=3,yend=string),color='white',size=0.1) +
  # plot main feautures on the line
  geom_boxplot(stat='identity',show.legend=F,position='identity',color='black',size=1,width=1) + ylab('Antigens') +
  # plot labels on top of antigen line segments (1.75 for hydropathy, 2.1 for epitope)
  geom_text(data=feature,mapping=aes(x=midpoint,y=2.1),fontface='bold',color='black',show.legend=F,size=7) +
  scale_fill_identity(guide = "legend", aes(labels = lab)) +
  scale_color_identity(guide = "legend", aes(labels = lab)) +
  scale_x_continuous(limits=c(1,xmax),breaks = seq(0,1000,100),expand=c(x1,x2)) +
  theme(text = element_text(size = 24)) + theme_bw() + ylab('BNT162b4') +
  theme(axis.text = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(angle=0,vjust=0.2,size=16,face='bold'),
        axis.ticks = element_blank(), rect = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# LABELS BOXES WITH THE EPITOPE ID NUMBER
p2.1 = ggplot(hit,aes(xmin=start,xlower=start,xmiddle=0,xupper=end,xmax=end,
                      y=fct_rev(crow),group=PlotGroup,fill=Color,color=TextColor,label=order)) +
  # plot the range of the observed epitopes
  geom_boxplot(stat='identity',position='identity',width=1,size=0.5,show.legend=F,color='black') + ylab('') +
  geom_text(mapping=aes(x=midpoint,y=fct_rev(crow)),size=5,fontface='bold') +
  scale_fill_identity(guide = "legend", aes(labels = lab)) +
  scale_color_identity(guide = "legend", aes(labels = lab)) +
  scale_x_continuous(limits=c(1,xmax),breaks = seq(0,1000,100),expand=c(x1,x2)) +
  scale_y_discrete(drop=F) +
  theme(text = element_text(size = 24)) + theme_bw() +
  theme(axis.text = element_blank(),axis.title.x = element_blank(), axis.title.y = element_text(angle=0,vjust=0.75,face='bold',size=16),
        axis.ticks = element_blank(),rect = element_blank(),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.minor.x = element_blank(), legend.position ='none')

margin = -0.2

# prints Figure 2A
windows(width=600,height=150)
plot(plot_grid(p1,NULL,p2.1,align='v',nrow=3,ncol=1,axis = 'lr',
               rel_heights = c(2,margin,-rowmin+1)))

# plot hydropathy
h = ggplot(hydro_df,aes(x=position,fill=score,y=string)) + geom_tile() +
  scale_fill_gradient(low='yellow',high='blue',breaks=c(-4.5,4.5),labels=c('Less','More'),
                      guide = guide_colorbar(direction='horizontal',title.position='top',vjust=0.7)) +
  scale_x_continuous(expand=c(x1,x2)) +
  theme_bw() + ylab('Hydropathy') + 
  labs(fill='Hydrophobicity') +
  theme(axis.text = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(angle=0,vjust=0.5,size=16,face='bold'),
        axis.ticks = element_blank(), rect = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = 'top',legend.justification = 'right',legend.title.align=0.5)
windows(width=600,height=100)
plot(h) # plot hydropathy only to generate legend (which will be manually overlaid onto composite plot)

# remove the legend for making the plotgrid
h = h + theme(legend.position='none')
windows(width=600,height=175)

# prints Figure S2D (manually add legend from plot(h))
plot(plot_grid(p1,NULL,h,NULL,p2.1,align='v',nrow=5,ncol=1,axis = 'lr',
               rel_heights = c(2,margin,0.5,margin,-rowmin-0.5)))

