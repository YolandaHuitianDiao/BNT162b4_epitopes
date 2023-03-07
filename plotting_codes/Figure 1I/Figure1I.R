rm(list=ls())
cat("\014") 
graphics.off()
options(scipen = 999)

library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)
library(ggrepel)
library(data.table)

vgrepl <- Vectorize(grepl, vectorize.args = "pattern")

setwd(dirname(envDocument::get_scriptpath()))
fontsize=10


# GENERATE SPIKE PLOTS ----------------------------------------------------
file = 'Spike Mapping for Figure 1I.txt'
vfile = 'Spike Conservation for Figure 1I.csv'

df_variant = fread(vfile,sep=',',header=T,stringsAsFactors=F) %>% as.data.frame() %>%
  mutate(label = paste(variant,mutation,sep='\n'))

### pre-process df_variant to remove consecutive deletions and count a single deletion
del_idx = grepl(pattern='>$',x=df_variant$mutation) # is mutation a deletion?
# is the mutation consecutive with one before it?
pos = df_variant$position 
consecutive_idx = c(FALSE,pos[2:length(pos)]-pos[1:(length(pos)-1)] == 1) 
# is the mutation a deletion AND consecutive?  e.g. is it the same single deletion?
consecutive_del_idx = del_idx & consecutive_idx
df_variant = df_variant[!consecutive_del_idx,] # then remove it except for the first call for plotting

feature = fread(file,header=T,stringsAsFactors=F,sep='\t') %>%
  mutate(Label = gsub(pattern='SP',replacement='S\nP',x=Antigen)) %>%
  mutate(FontSize = fontsize*as.numeric(FontSize))
xmax = 1273

map_plot_spike = ggplot(feature,aes(xmin=start,xlower=start,xmiddle=0,xupper=end,xmax=end,x=midpoint,
                                    y=string,group=PlotGroup,fill=Color,label=Label,color=TextColor)) +
  geom_segment(data=feature,aes(x=1,xend=xmax,y=string,yend=string),color='black',size=2) +
  geom_boxplot(data=feature%>%filter(Antigen!='RBD'),stat='identity',show.legend=F,position='identity',color='black',size=1,width=1) +
  geom_boxplot(data=feature%>%filter(Antigen=='RBD'),stat='identity',show.legend=F,position='identity',color='black',size=1,width=1) +
  scale_fill_identity(guide = "legend", aes(labels = lab)) +
  scale_color_identity(guide = "legend", aes(labels = lab)) +
  geom_text(show.legend=F,size=fontsize,fontface='bold') +
  scale_x_continuous(limits=c(1,xmax),breaks = seq(0,1200,100),expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_bw() + ylab('Spike') +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=35,face='bold',color='black',angle=0,vjust=0.5),
        axis.ticks = element_blank(), rect = element_blank(),axis.text.y=element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position='none')

variant_plot_spike = ggplot(df_variant,aes(x=position,y=fct_rev(variant),fill=variant,label=label,color=variant)) + 
  geom_jitter(shape=25,show.legend=T,size=7,width=0,height=0,color='black') + 
  theme_bw() +
  scale_fill_viridis_d() +
  ylab('Variant') +  scale_x_continuous(limits=c(1,xmax),expand=c(0,0)) +
  scale_y_discrete(expand=c(0.2,0)) + 
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        axis.ticks = element_blank(), rect = element_blank(),axis.text.y=element_text(size=24,face='plain',color='black'),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position='none')


# GENERATE BNT162B4 PLOTS -------------------------------------------------
file = 'BNT162b4 Mapping for Figure 1I.txt'
vfile = 'BNT162b4 Conservation for Figure 1I.csv'
df_variant = fread(vfile,sep=',',header=T,stringsAsFactors=F) %>% as.data.frame() %>%
  mutate(label = paste(variant,mutation,sep='\n'))

variantlist = df_variant %>% select(variant) %>% filter(variant!='Omicron (All)') %>%
  as.matrix()
variantlist = c(variantlist,'Omicron (All)')

scale_fill_variant <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(viridisLite::viridis(n=16)[c(1,3,8,16)],
                      variantlist), 
    ...
  )
} # pegs each variant mutation marker to a color

feature = fread(file,header=T,stringsAsFactors=F,sep='\t') %>%
  mutate(Label = gsub(pattern='ORF1ab',replacement='',x=Antigen)) %>%
  mutate(FontSize = fontsize*as.numeric(FontSize))

map_plot_b4 = ggplot(feature,aes(xmin=start,xlower=start,xmiddle=0,xupper=end,xmax=end,x=midpoint,size=FontSize,
                                 y=string,group=PlotGroup,fill=Color,color=TextColor,label=Label)) +
  geom_segment(data=feature,aes(x=1,xend=598,y=string,yend=string),color='black',size=2) +
  scale_fill_identity(guide = "legend", aes(labels = lab)) +
  scale_color_identity(guide = "legend", aes(labels = lab)) +
  geom_boxplot(stat='identity',show.legend=F,position='identity',color='black',size=1,width=1) + ylab('Antigens') +
  scale_x_continuous(limits=c(1,xmax),breaks = seq(0,1000,100),expand=c(0,0)) +
  theme_bw() + ylab('BNT162b4') +
  theme(axis.text = element_blank(), axis.title.x = element_blank(),axis.title.y=element_text(size=35,face='bold',color='black',angle=0,vjust=0.5),
        axis.ticks = element_blank(), rect = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

variant_plot_b4 = ggplot(df_variant,aes(x=position,y='Mutation',fill=variant,label=label)) + 
  geom_point(shape=25,show.legend=F,size=7.5) + 
  geom_text_repel(nudge_y=0.5,box.padding=1,fontface='plain',size=7.5) +
  theme_bw() +
  scale_fill_variant() +
  scale_x_continuous(limits=c(1,xmax),expand=c(0,0)) +
  ylab('Variant') +
  theme(axis.title = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), rect = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

r = plot_grid(variant_plot_spike,NULL,map_plot_spike,NULL,
              variant_plot_b4,NULL,map_plot_b4,
              nrow=7,ncol=1,align='v',axis='lr',rel_heights=c(6,-0.8,1.7,0.5,3,-1.5,2))
windows(width=1820,height=800)
plot(r) # maximize window







