################################################################################
#NN-MM no middle v.s. pblup
################################################################################
setwd("/Users/tianjing/Library/CloudStorage/Box-Box/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker")

nnmm_no=read.table("nnmm_all_ind_nongenotyped/ALL_accuracy.nnmm.txt")$V2
pblup=read.table("pblup_all_ind_nongenotyped/ALL_accuracy.pblup.txt")$V2
df_no   = data.frame(nnmm_no,pblup)

nnmm_all=read.table("nnmm_all_ind_genotyped/ALL_accuracy.nnmm_allomics.txt")$V2
gblup=read.table("gblup_all_ind_genotyped/ALL_accuracy.gblup.txt")$V2
df_all   = data.frame(nnmm_all,gblup)

library(gridExtra)
library(ggplot2)

myalpha=0.4
p1=ggplot() +
  geom_point(data=df_no,aes(x=pblup,y = nnmm_no)) + 
  xlab("PBLUP") +
  ylab("NN-MM") +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  xlim(0.85,0.895)+
  ylim(0.85,0.895)

p2=ggplot() +
  geom_point(data=df_all,aes(x=gblup,y = nnmm_all)) + 
  xlab("GBLUP") +
  ylab("NN-MM") +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  xlim(0.985,1)+
  ylim(0.985,1)

p_all=grid.arrange(p2, p1, nrow = 1)











