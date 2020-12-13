
#setwd("C:\Users\guowc_000\Documents\Work\HealthCard\epiBMI\Diabetes\epiBMI_sites_addedToAGP2_21")
rm(list=ls())    # clean the work folder.

#dat <- read.csv("GSE87571_62003_76169_EPIC1_120_in2212_2270_epiBMI19of20sites_beta.csv", row.names=1, stringsAsFactors=F)
dat <- read.csv("GSE87571_62003_76169_EPIC1_120_in2212_2270_epiBMI15of20sites_beta.csv", row.names=1, stringsAsFactors=F)
d <- t(dat)
cntrl <- d[1:732,]
t2d <- d[733:790,]
t1d <- d[791:853,]
epic <- d[854:973,]
in2212 <- d[974:1023,]
in2270 <- d[1024:1077,]

sink("mean_sd_GSE87571_15diabetesCpGs_beta.txt")
apply(cntrl,2,mean)
apply(cntrl,2,sd)
apply(t2d,2,mean)
apply(t2d,2,sd)
apply(t1d,2,mean)
apply(t1d,2,sd)
sink()

cntrl_mean <- apply(cntrl,2,mean)
cntrl_sd <- apply(cntrl,2,sd)
t2d_mean <- apply(t2d,2,mean)
t2d_sd <- apply(t2d,2,sd)

# Find and discard CpGs with |effect size| <= 0.2
#diff_cntrl_t2d_means <- t2d_mean - cntrl_mean
#write.csv(diff_cntrl_t2d_means, "diff_cntrl_t2d_means.csv",row.names=F)
#effect_size_cntrlSD <- diff_cntrl_t2d_means / cntrl_sd
#effect_size_t2dSD <- diff_cntrl_t2d_means / t2d_sd
#write.csv(effect_size_cntrlSD, "effect_size_byCntrlSD.csv",row.names = F) 
#write.csv(effect_size_t2dSD, "effect_size_byT2dSD.csv",row.names = F)

get_zscore_site <- function(x, m, s) {
  return ((x - m) / s)
}  

Sites = 15
direction <- c(1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1) # 15 CpGs

#################################################### t2d58     
zscores <- data.frame()
for( i in 1:Sites) {
  for (j in 1:length(t2d[,1])) {
    zscores[j,i] <- get_zscore_site(t2d[j,i], cntrl_mean[i], cntrl_sd[i])
  }
}

  z = zscores
  zscores_people <- vector()
  for(rw in 1:dim(z)[1]) { # person by person
    sum = 0
    for(cl in 1:Sites) { # site by site 
      # if a site is negatively correlated with t2d.
      if(direction[cl] == -1) {
        sum = sum - z[rw,cl]
      # if a site is positively correlated with t2d. 
      } else { 
        sum = sum + z[rw,cl]
      } # end of if
    }   # end of the for loop of 15 sites
    
    zscores_people[rw] = sum
  }     # end of the for loop of all samples

zscores_people_t2d58 <- zscores_people
write.csv(zscores_people_t2d58,"zscores_people_t2d58.csv")
write.csv(zscores,"zscores_sites_t2d58.csv")

##################### t1d63
zscores <- data.frame()
for( i in 1:Sites) {
  for (j in 1:length(t1d[,1])) {
    zscores[j,i] <- get_zscore_site(t1d[j,i], cntrl_mean[i], cntrl_sd[i])
  }
}

z = zscores
zscores_people <- vector()
for(rw in 1:dim(z)[1]) { # person by person
  sum = 0
  for(cl in 1:Sites) { # site by site 
    # if a site is negatively correlated with t2d.
    if(direction[cl] == -1) {
      sum = sum - z[rw,cl]
      # if a site is positively correlated with t2d. 
    } else { 
      sum = sum + z[rw,cl]
    } # end of if
  }   # end of the for loop of 15 sites
  
  zscores_people[rw] = sum
}     # end of the for loop of all samples

zscores_people_t1d63 <- zscores_people
write.csv(format(zscores_people_t1d63,digits=4),"zscores_people_t1d63.csv")
write.csv(format(zscores,digits=4),"zscores_sites_t1d63.csv")

################################################## zymoEPIC1_120
zscores <- data.frame()
for( i in 1:Sites) {
  for (j in 1:length(epic[,1])) {
    zscores[j,i] <- get_zscore_site(epic[j,i], cntrl_mean[i], cntrl_sd[i])
  }
}

z = zscores
zscores_people <- vector()
for(rw in 1:dim(z)[1]) { # person by person
  sum = 0
  for(cl in 1:Sites) { # site by site 
    # if a site is negatively correlated with t2d.
    if(direction[cl] == -1) {
      sum = sum - z[rw,cl]
      # if a site is positively correlated with t2d. 
    } else { 
      sum = sum + z[rw,cl]
    } # end of if
  }   # end of the for loop of 15 sites
  
  zscores_people[rw] = sum
}     # end of the for loop of all samples

zscores_people_epic120 <- zscores_people
write.csv(format(zscores_people_epic120,digits=4),"zscores_people_epic120.csv")
write.csv(format(zscores,digits=4),"zscores_sites_epic120.csv")


############################################### in2212 50 samples
zscores <- data.frame()
for( i in 1:Sites) {
  for (j in 1:length(in2212[,1])) {
    zscores[j,i] <- get_zscore_site(in2212[j,i], cntrl_mean[i], cntrl_sd[i])
  }
}

z = zscores
zscores_people <- vector()
for(rw in 1:dim(z)[1]) { # person by person
  sum = 0
  for(cl in 1:Sites) { # site by site 
    # if a site is negatively correlated with t2d.
    if(direction[cl] == -1) {
      sum = sum - z[rw,cl]
      # if a site is positively correlated with t2d. 
    } else { 
      sum = sum + z[rw,cl]
    } # end of if
  }   # end of the for loop of 15 sites
  
  zscores_people[rw] = sum
}     # end of the for loop of all samples

zscores_people_in2212 <- zscores_people
write.csv(format(zscores_people_in2212,digits=4),"zscores_people_in2212.csv")
write.csv(format(zscores,digits=4),"zscores_sites_in2212.csv")


######################################### in2270 54 samples
zscores <- data.frame()
for( i in 1:Sites) {
  for (j in 1:length(in2270[,1])) {
    zscores[j,i] <- get_zscore_site(in2270[j,i], cntrl_mean[i], cntrl_sd[i])
  }
}

z = zscores
zscores_people <- vector()
for(rw in 1:dim(z)[1]) { # person by person
  sum = 0
  for(cl in 1:Sites) { # site by site 
    # if a site is negatively correlated with t2d.
    if(direction[cl] == -1) {
      sum = sum - z[rw,cl]
      # if a site is positively correlated with t2d. 
    } else { 
      sum = sum + z[rw,cl]
    } # end of if
  }   # end of the for loop of 15 sites
  
  zscores_people[rw] = sum
}     # end of the for loop of all samples

zscores_people_in2270 <- zscores_people
write.csv(format(zscores_people_in2270,digits=4),"zscores_people_in2270.csv")
write.csv(format(zscores,digits=4),"zscores_sites_in2270.csv")

############################################### control GSE87571 732 samples
zscores <- data.frame()
for( i in 1:Sites) {
  for (j in 1:length(cntrl[,1])) {
    zscores[j,i] <- get_zscore_site(cntrl[j,i], cntrl_mean[i], cntrl_sd[i])
  }
}

z = zscores
zscores_people <- vector()
for(rw in 1:dim(z)[1]) { # person by person
  sum = 0
  for(cl in 1:Sites) { # site by site 
    # if a site is negatively correlated with t2d.
    if(direction[cl] == -1) {
      sum = sum - z[rw,cl]
      # if a site is positively correlated with t2d. 
    } else { 
      sum = sum + z[rw,cl]
    } # end of if
  }   # end of the for loop of 15 sites
  
  zscores_people[rw] = sum
}     # end of the for loop of all samples

zscores_people_gse87571 <- zscores_people
write.csv(format(zscores_people_gse87571,digits=4),"zscores_people_gse87571.csv")
write.csv(format(zscores,digits=4),"zscores_sites_gse87571.csv")

cpgs15 <- c("cg06500161","cg11024682","cg00574958","cg03725309","cg06192883","cg06690548","cg12593793","cg14476101","cg17058475","cg16246545","cg18181703","cg25649826","cg19693031","cg22891070","cg07814318")
colnames(zscores) <- cpgs15
cntrl_sites_zscores_summary <- summary(zscores)
write.csv(cntrl_sites_zscores_summary,"cntrl_sites_zscores_summary.csv",row.names = FALSE)

######################################## Plot gen. population, t2d, no diabetes, customers
library(ggpubr)
input_plot = read.csv("zscores_people_forPlot_gse87571_gse62003_gse76169_x_113020.csv")
png("Scatterplot_15-CpG_T2D_Zscore_MRS.png")
sp <- ggscatter(input_plot, x= "age",y="MRS",color="group",
          palette=c("red","blue","gray","orange"),
          xlab="Age. The red line denotes that 95% people without diabetes had MRS < 8.2.", ylab="T2D Methylation Risk Score")
sp + geom_hline(yintercept=8.2, color="red",size=1.2)
dev.off()

###########################################################
#t.test(zscores_people_tests58,zscores_people_t2d58)
#t.test(zscores_people_tests58,zscores_people_t1d63)
#t.test(zscores_people_t2d58,zscores_people_t1d63)

zscores <- read.csv("zscores_sites_t2d58.csv",row.names = 1)
colnames(zscores) <- cpgs15
t2d58_sites_zscores_summary <- summary(zscores)
write.csv(t2d58_sites_zscores_summary,"t2d58_sites_zscores_summary.csv",row.names = FALSE)

get95CI <- function(x) { # normal distribution assumed
  #error <- qnorm(0.975) * sd(x) / sqrt(length(x))
  #left <- mean(x) - error
  left <- quantile(x,probs=c(0.025, 0.975))[1]
  right <- quantile(x,probs=c(0.025, 0.975))[2]
  res <- c(left, right)
  return (res)
}


zscores <- read.csv("zscores_sites_gse87571.csv",row.names = 1)
resl <- apply(zscores, 2, get95CI)
write.csv(format(resl,digits=4),"zscores_sites_95CI_gse87571.csv",row.names = FALSE)

