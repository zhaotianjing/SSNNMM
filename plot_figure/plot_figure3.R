################################################################################
#conventional single-step
################################################################################
setwd("/Users/tianjing/Library/CloudStorage/Box-Box/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker")

ptc50_ss=read.table("0.5/accuracy_ss/ALL_accuracy.ss.pct0.5.txt")$V2
a50_ss=mean(ptc50_ss)
sd50_ss=sd(ptc50_ss)

ptc70_ss=read.table("0.7/accuracy_ss/ALL_accuracy.ss.pct0.7.txt")$V2
a70_ss=mean(ptc70_ss)
sd70_ss=sd(ptc70_ss)

ptc90_ss=read.table("0.9/accuracy_ss/ALL_accuracy.ss.pct0.9.txt")$V2
a90_ss=mean(ptc90_ss)
sd90_ss=sd(ptc90_ss)

ptc99_ss=read.table("0.99/accuracy_ss/ALL_accuracy.ss.pct0.99.txt")$V2
a99_ss=mean(ptc99_ss)
sd99_ss=sd(ptc99_ss)


################################################################################
#conventional gblup (genotyped inds only)
################################################################################
ptc50_gb=read.table("0.5/accuracy_gblup/ALL_accuracy.gblup.pct0.5.txt")$V2
a50_gb=mean(ptc50_gb)
sd50_gb=sd(ptc50_gb)

ptc70_gb=read.table("0.7/accuracy_gblup/ALL_accuracy.gblup.pct0.7.txt")$V2
a70_gb=mean(ptc70_gb)
sd70_gb=sd(ptc70_gb)

ptc90_gb=read.table("0.9/accuracy_gblup/ALL_accuracy.gblup.pct0.9.txt")$V2
a90_gb=mean(ptc90_gb)
sd90_gb=sd(ptc90_gb)

ptc99_gb=read.table("0.99/accuracy_gblup/ALL_accuracy.gblup.pct0.99.txt")$V2
a99_gb=mean(ptc99_gb)
sd99_gb=sd(ptc99_gb)

################################################################################
#NN-MM
################################################################################
ptc50_nnmm=read.table("0.5/accuracy_nnmm_pblup_zw1/ALL_accuracy.nnmm_pblup_zw1.pct0.5.txt")$V2
a50_nnmm=mean(ptc50_nnmm)
sd50_nnmm=sd(ptc50_nnmm)

ptc70_nnmm=read.table("0.7/accuracy_nnmm_pblup_zw1/ALL_accuracy.nnmm_pblup_zw1.pct0.7.txt")$V2
a70_nnmm=mean(ptc70_nnmm)
sd70_nnmm=sd(ptc70_nnmm)

ptc90_nnmm=read.table("0.9/accuracy_nnmm_pblup_zw1/ALL_accuracy.nnmm_pblup_zw1.pct0.9.txt")$V2
a90_nnmm=mean(ptc90_nnmm)
sd90_nnmm=sd(ptc90_nnmm)

ptc99_nnmm=read.table("0.99/accuracy_nnmm_pblup_zw1/ALL_accuracy.nnmm_pblup_zw1.pct0.99.txt")$V2
a99_nnmm=mean(ptc99_nnmm)
sd99_nnmm=sd(ptc99_nnmm)


#pblup results (lower bound)
pblup=0.642

#full genotype: gblup results  (upper bound)
full_gblup=0.994

#ss results
a_ss_all=c(a50_ss, a70_ss, a90_ss, a99_ss)
sd_ss_all=c(sd50_ss, sd70_ss, sd90_ss, sd99_ss)/sqrt(10)

#gblup results
a_gb_all=c(a50_gb, a70_gb, a90_gb, a99_gb)
sd_gb_all=c(sd50_gb, sd70_gb, sd90_gb, sd99_gb)/sqrt(10)

#nnmm results
a_nnmm_all=c(a50_nnmm, a70_nnmm, a90_nnmm, a99_nnmm)
sd_nnmm_all=c(sd50_nnmm, sd70_nnmm, sd90_nnmm, sd99_nnmm)/sqrt(10)

pct=c(50, 70, 90,99)
xLabels=c("50","70","90","99")


df_ss   = data.frame(pct,a_ss_all, sd_ss_all)
df_gb   = data.frame(pct,a_gb_all, sd_gb_all)
df_nnmm = data.frame(pct,a_nnmm_all, sd_nnmm_all)

library(ggplot2)
myalpha=0.4
ggplot() +
  geom_line(data=df_ss,aes(x=pct,y = a_ss_all, group = 1,colour ="Conventional single-step"),alpha = myalpha,linetype="solid") +
  geom_point(data=df_ss,aes(x=pct,y = a_ss_all, group = 1,colour ="Conventional single-step"),alpha = myalpha+0.3) +
  geom_errorbar(data=df_ss,aes(x=pct,ymin=a_ss_all-sd_ss_all, ymax=a_ss_all+sd_ss_all,colour ="Conventional single-step"), width=.2,
                alpha =myalpha)+
  geom_line(data=df_gb,aes(x=pct,y = a_gb_all, group = 1,colour ="GBLUP (genotyped individuals only)"),alpha = myalpha) +
  geom_point(data=df_gb,aes(x=pct,y = a_gb_all, group = 1,colour ="GBLUP (genotyped individuals only)"),alpha = myalpha+0.3) +
  geom_errorbar(data=df_gb,aes(x=pct, ymin=a_gb_all-sd_gb_all, ymax=a_gb_all+sd_gb_all,colour ="GBLUP (genotyped individuals only)"), width=.2,
                  alpha = myalpha)+
  geom_line(data=df_gb,aes(x=pct,y = a_nnmm_all, group = 1,colour ="NN-MM"),alpha = myalpha) +
  geom_point(data=df_gb,aes(x=pct,y = a_nnmm_all, group = 1,colour ="NN-MM"),alpha = myalpha+0.3) +
  geom_errorbar(data=df_gb,aes(x=pct, ymin=a_nnmm_all-sd_nnmm_all, ymax=a_nnmm_all+sd_nnmm_all,colour ="NN-MM"), width=.2,
                alpha = myalpha)+
  xlab("% non-genotyped individuals in training dataset") +
  ylab("Prediction accuracy") +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0.6,1.0)+
  scale_x_continuous(breaks = pct,labels = xLabels)+
  scale_colour_manual("", breaks = c("Conventional single-step", "NN-MM","GBLUP (genotyped individuals only)"),values = c("blue","red","black"))+
  geom_hline(yintercept=pblup, linetype="dashed", color = "black", size=0.7)+ #PBLUP
  geom_hline(yintercept=full_gblup, linetype="dashed", color = "black", size=0.7) #PBLUP




