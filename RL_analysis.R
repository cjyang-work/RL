library(stringr)
library(sommer)
library(mvtnorm)
library(car)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggradar)
library(ggh4x)
library(emmeans)
library(sf)

# read the marker genotype data.
geno <- readRDS("input/geno.RDS")

# read the phenotypic trait data.
pheno <- read.csv("input/pheno.csv", as.is=TRUE)

# check for number of varieties tested per year.
t(t(colSums(table(pheno[, c("Variety", "Year")]) > 0)))
#2002   19
#2003   19
#2004   25
#2005   27
#2006   29
#2007   24
#2008   22
#2009   21
#2010   22
#2011   22
#2012   29
#2013   34
#2014   38
#2015   24
#2016   20
#2017   21
#2018   25
#2019   21
#2020   14
#2021   11
#2022   10

# check for number of varieties first appear in each year.
t(t(table(tapply(pheno$Year, pheno$Variety, min))))
#2002   19
#2003    9
#2004   11
#2005    7
#2006    5
#2007    7
#2008    6
#2009    6
#2010    8
#2011    7
#2012   14
#2013   12
#2014   13
#2015    4
#2016    5
#2017    3
#2018    9
#2019    8

## use only 2002 to 2019 - no new varieties in 2020/2021/2022.
## require 5-year minimum so start the prediction on the 6th year, i.e. 2007.

# set variable.
year <- 2007:2019
temp <- unique(pheno[, c("Site", "Region")])
temp <- temp[order(temp$Region, temp$Site), ]
region <- temp$Region
names(region) <- temp$Site
first <- tapply(pheno$Year, pheno$Variety, min)

## Loop through the years to get the predictions.
out.old <- out.new <- list()
for(i in 1:length(year)){
  
  # data from year[i]-5 to year[i]-1.
  df <- pheno[pheno$Year%in%c((year[i]-5):(year[i]-1)), ]
  df$Variety <- factor(x=df$Variety, levels=sort(unique(df$Variety)))
  df$Year <- factor(x=df$Year, levels=sort(unique(df$Year)))
  df$Site <- factor(x=df$Site, levels=sort(unique(df$Site)))
  df$Region <- factor(x=df$Region, levels=sort(unique(df$Region)))
  
  # identify new variety in year[i].
  variety0 <- names(first)[first==year[i]]
  
  # identify all varieties in year[i]-5 to year[i].
  variety <- c(levels(df$Variety), variety0)
  
  ### simple means (current AHDB methods).
  # trait to predict in year[i].
  target <- acast(data=pheno[pheno$Year==year[i], ], formula=Variety~Site, value.var="Yield", drop=TRUE)

  # individual site, year[i]-1.
  ind1 <- acast(data=pheno[pheno$Year==year[i]-1, ], formula=Variety~Site, value.var="Yield", drop=TRUE)
  
  # individual site, mean from year[i]-5 to year[i]-1.
  ind5 <- acast(data=df, formula=Variety~Site, value.var="Yield", drop=TRUE, fun.aggregate=mean, na.rm=TRUE)
  
  # all sites, year[i].
  all1 <- rowSums(ind1, na.rm=TRUE)/rowSums(!is.na(ind1))
  
  # all sites, mean from year[i]-5 to year[i]-1.
  all5 <- acast(data=df, formula=Variety~Year, value.var="Yield", drop=TRUE, fun.aggregate=mean, na.rm=TRUE)
  all5 <- rowSums(all5, na.rm=TRUE)/rowSums(!is.na(all5))
  
  # regional site, year[i]-1.
  reg1 <- acast(data=pheno[pheno$Year==year[i]-1, ], formula=Variety~Region, value.var="Yield", drop=TRUE, fun.aggregate=mean, na.rm=TRUE)
  
  # regional site, mean from year[i]-5 to year[i]-1.
  reg5 <- acast(data=df, formula=Variety~Region+Year, value.var="Yield", drop=TRUE, fun.aggregate=mean, na.rm=TRUE)
  reg5 <- cbind(East=rowSums(reg5[,1:5], na.rm=TRUE)/rowSums(!is.na(reg5[,1:5])),
                North=rowSums(reg5[,6:10], na.rm=TRUE)/rowSums(!is.na(reg5[,6:10])),
                West=rowSums(reg5[,11:15], na.rm=TRUE)/rowSums(!is.na(reg5[,11:15])))
  
  ### 5-year GBLUP.
  # model 1: individual site with GxY.
  # identify Site with data in year[i], and at least 3 years of data in year[i]-5 to year[i]-1.
  site0 <- table(pheno[pheno$Year%in%c((year[i]-5):(year[i])), c("Site", "Year")])
  site0 <- rownames(site0)[site0[,6] > 0 & rowSums(site0 > 0) > 3]
  
  # loop.
  ind.blup1 <- matrix(NA, nrow=length(variety)-length(variety0), ncol=length(site0), dimnames=list(variety[!(variety%in%variety0)], site0))
  ind.blup2 <- matrix(NA, nrow=length(variety0), ncol=length(site0), dimnames=list(variety0, site0))
  for(j in 1:length(site0)){
    
    # subset the data.
    df0 <- droplevels(df[df$Site==site0[j], ])
    levels(df0$Variety) <- c(levels(df0$Variety), variety0)
    
    # calculate A.
    A0 <- A.mat(X=geno[levels(df0$Variety),]-1)
    
    # create Y.
    Y0 <- diag(length(levels(df0$Year)))
    rownames(Y0) <- colnames(Y0) <- levels(df0$Year)
    
    # create YA.
    YA0 <- kronecker(Y0, A0, make.dimnames=TRUE)
    
    # model.
    mm0 <- mmer(fixed=Yield~Year,
                random=~vsr(Variety, Gu=A0) + vsr(Year:Variety, Gu=YA0),
                rcov=~units,
                data=df0,
                dateWarning=FALSE,
                verbose=FALSE)
    
    # get the BLUP.
    temp <- mm0$U$`u:Variety`$Yield
    temp1 <- temp[!(names(temp)%in%variety0)]
    temp2 <- temp[names(temp)%in%variety0]
    ind.blup1[names(temp1), site0[j]] <- temp1
    ind.blup2[names(temp2), site0[j]] <- temp2
    
  }
  
  # model 2: all sites without GxY or GxC or GxYxC.
  # add levels for new varieties.
  df0 <- df
  levels(df0$Variety) <- variety
  
  # calculate A.
  A0 <- A.mat(X=geno[variety,]-1)
  
  # model.
  mm0 <- mmer(fixed=Yield~Year+Site,
              random=~vsr(Variety, Gu=A0),
              rcov=~units,
              data=df0,
              dateWarning=FALSE,
              verbose=FALSE)
  
  # get the BLUP.
  temp <- mm0$U$`u:Variety`$Yield
  all.blup1 <- temp[!(names(temp)%in%variety0)]
  all.blup2 <- temp[names(temp)%in%variety0]
  
  # identify site that is found in target.
  site.common <- colnames(target)
  site.common <- site.common[site.common%in%colnames(ind1) & site.common%in%colnames(ind.blup1)]
  
  # identify old/new variety that is found in target.
  variety.old <- rownames(target)
  variety.old <- variety.old[variety.old%in%rownames(ind1)]
  variety.new <- variety0
  
  # expand reg1 and reg5 from region to site.
  reg1 <- reg1[, region[site.common]]
  colnames(reg1) <- site.common
  reg5 <- reg5[, region[site.common]]
  colnames(reg5) <- site.common
  
  # combine all results into a list, and retain only site/variety found in target.
  temp.old <- list(target=target[variety.old, site.common],
                   ind1=ind1[variety.old, site.common],
                   ind5=ind5[variety.old, site.common],
                   all1=all1[variety.old],
                   all5=all5[variety.old],
                   reg1=reg1[variety.old, site.common],
                   reg5=reg5[variety.old, site.common],
                   ind.blup=ind.blup1[variety.old, site.common],
                   all.blup=all.blup1[variety.old])
  temp.new <- list(target=target[variety.new, site.common],
                   ind.blup=ind.blup2[variety.new, site.common],
                   all.blup=all.blup2[variety.new])
  
  # collect the results.
  out.old <- c(out.old, list(temp.old))
  out.new <- c(out.new, list(temp.new))
  
  message(i)
  
}

save(out.old, out.new, file="output/out_20230821.RData")

# set up empty matrices to store Pearson correlations.
temp <- lapply(1:length(year), FUN=function(i) colnames(out.old[[i]]$target))
temp <- sort(unique(unlist(temp)))
temp <- matrix(NA, nrow=length(temp), ncol=length(year), dimnames=list(temp, year))
cor.ind1 <- cor.ind5 <- cor.reg1 <- cor.reg5 <- cor.all1 <- cor.all5 <- cor.ind.blup <- cor.all.blup <- temp

# loop to calculate Pearson correlations.
for(i in 1:length(year)){
  
  # correlations between target and ind1.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$ind1, use="pairwise.complete.obs"))
  cor.ind1[names(temp), i] <- temp
  
  # correlations between target and ind5.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$ind5, use="pairwise.complete.obs"))
  cor.ind5[names(temp), i] <- temp
  
  # correlations between target and reg1.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$reg1, use="pairwise.complete.obs"))
  cor.reg1[names(temp), i] <- temp
  
  # correlations between target and reg5.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$reg5, use="pairwise.complete.obs"))
  cor.reg5[names(temp), i] <- temp
  
  # correlations between target and all1.
  temp <- cor(out.old[[i]]$target, out.old[[i]]$all1, use="pairwise.complete.obs")[,1]
  cor.all1[names(temp), i] <- temp
  
  # correlations between target and all5.
  temp <- cor(out.old[[i]]$target, out.old[[i]]$all5, use="pairwise.complete.obs")[,1]
  cor.all5[names(temp), i] <- temp
  
  # correlations between target and ind.blup.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$ind.blup, use="pairwise.complete.obs"))
  cor.ind.blup[names(temp), i] <- temp
  
  # correlations between target and all.blup.
  temp <- cor(out.old[[i]]$target, out.old[[i]]$all.blup, use="pairwise.complete.obs")[,1]
  cor.all.blup[names(temp), i] <- temp
  
}

# count the number of varieties used for correlating.
temp <- lapply(1:length(year), FUN=function(i) colnames(out.old[[i]]$target))
temp <- sort(unique(unlist(temp)))
temp <- matrix(NA, nrow=length(temp), ncol=length(year), dimnames=list(temp, year))
cor.n <- replicate(n=8, expr=list(temp))
names(cor.n) <- c("ind1", "ind5", "all1", "all5", "reg1", "reg5", "ind.blup", "all.blup")

for(i in 1:length(year)){
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$ind1))
  cor.n$ind1[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$ind5))
  cor.n$ind5[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$all1))
  cor.n$all1[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$all5))
  cor.n$all5[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$reg1))
  cor.n$reg1[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$reg5))
  cor.n$reg5[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$ind.blup))
  cor.n$ind.blup[names(temp), i] <- temp
  
  temp <- colSums(!is.na(out.old[[i]]$target) & !is.na(out.old[[i]]$all.blup))
  cor.n$all.blup[names(temp), i] <- temp
  
}

# save the Pearson correlations.
write.csv(cor.ind1, "output/results/cor_ind1.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.ind5, "output/results/cor_ind5.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.all1, "output/results/cor_all1.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.all5, "output/results/cor_all5.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.reg1, "output/results/cor_reg1.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.reg5, "output/results/cor_reg5.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.ind.blup, "output/results/cor_ind.blup.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.all.blup, "output/results/cor_all.blup.csv", row.names=TRUE, quote=FALSE)

write.csv(cor.n$ind1, "output/results/n_ind1.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$ind5, "output/results/n_ind5.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$all1, "output/results/n_all1.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$all5, "output/results/n_all5.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$reg1, "output/results/n_reg1.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$reg5, "output/results/n_reg5.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$ind.blup, "output/results/n_ind.blup.csv", row.names=TRUE, quote=FALSE)
write.csv(cor.n$all.blup, "output/results/n_all.blup.csv", row.names=TRUE, quote=FALSE)

# set up empty matrices to store Spearman (rank-based) correlations.
temp <- lapply(1:length(year), FUN=function(i) colnames(out.old[[i]]$target))
temp <- sort(unique(unlist(temp)))
temp <- matrix(NA, nrow=length(temp), ncol=length(year), dimnames=list(temp, year))
cors.ind1 <- cors.ind5 <- cors.reg1 <- cors.reg5 <- cors.all1 <- cors.all5 <- cors.ind.blup <- cors.all.blup <- temp

# loop to calculate Spearman correlations.
for(i in 1:length(year)){
  
  # correlations between target and ind1.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$ind1, use="pairwise.complete.obs", method="spearman"))
  cors.ind1[names(temp), i] <- temp
  
  # correlations between target and ind5.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$ind5, use="pairwise.complete.obs", method="spearman"))
  cors.ind5[names(temp), i] <- temp
  
  # correlations between target and reg1.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$reg1, use="pairwise.complete.obs", method="spearman"))
  cors.reg1[names(temp), i] <- temp
  
  # correlations between target and reg5.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$reg5, use="pairwise.complete.obs", method="spearman"))
  cors.reg5[names(temp), i] <- temp
  
  # correlations between target and all1.
  temp <- cor(out.old[[i]]$target, out.old[[i]]$all1, use="pairwise.complete.obs", method="spearman")[,1]
  cors.all1[names(temp), i] <- temp
  
  # correlations between target and all5.
  temp <- cor(out.old[[i]]$target, out.old[[i]]$all5, use="pairwise.complete.obs", method="spearman")[,1]
  cors.all5[names(temp), i] <- temp
  
  # correlations between target and ind.blup.
  temp <- diag(cor(out.old[[i]]$target, out.old[[i]]$ind.blup, use="pairwise.complete.obs", method="spearman"))
  cors.ind.blup[names(temp), i] <- temp
  
  # correlations between target and all.blup.
  temp <- cor(out.old[[i]]$target, out.old[[i]]$all.blup, use="pairwise.complete.obs", method="spearman")[,1]
  cors.all.blup[names(temp), i] <- temp
  
}

# save the Spearman correlations.
write.csv(cors.ind1, "output/results/cors_ind1.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.ind5, "output/results/cors_ind5.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.all1, "output/results/cors_all1.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.all5, "output/results/cors_all5.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.reg1, "output/results/cors_reg1.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.reg5, "output/results/cors_reg5.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.ind.blup, "output/results/cors_ind.blup.csv", row.names=TRUE, quote=FALSE)
write.csv(cors.all.blup, "output/results/cors_all.blup.csv", row.names=TRUE, quote=FALSE)

### identify yield of best predicted variety,
### best observed variety, worst observed variety for year-i, site-j.
# method
method <- c("I1", "I5", "I5B", "R1", "R5", "A1", "A5", "A5B")

# loop through year-i and site-j.
bv <- vector()
for(i in 1:length(year)){
  for(j in 1:ncol(out.old[[i]]$target)){
    
    temp <- out.old[[i]]$target[,j]
    temp <- temp[!is.na(temp)]
    
    temp.ind1 <- out.old[[i]]$ind1[names(temp), j]
    temp.ind5 <- out.old[[i]]$ind5[names(temp), j]
    
    temp.ind.blup <- out.old[[i]]$ind.blup[names(temp), j]
    
    temp.reg1 <- out.old[[i]]$reg1[names(temp), j]
    temp.reg5 <- out.old[[i]]$reg5[names(temp), j]
    
    temp.all1 <- out.old[[i]]$all1[names(temp)]
    temp.all5 <- out.old[[i]]$all5[names(temp)]
    
    temp.all.blup <- out.old[[i]]$all.blup[names(temp)]
    
    temp.variety <- c(names(which.max(temp.ind1)),
                      names(which.max(temp.ind5)),
                      names(which.max(temp.ind.blup)),
                      names(which.max(temp.reg1)),
                      names(which.max(temp.reg5)),
                      names(which.max(temp.all1)),
                      names(which.max(temp.all5)),
                      names(which.max(temp.all.blup)))
    
    bv <- rbind(bv,
                data.frame(Year=year[i],
                           Site=colnames(out.old[[i]]$target)[j],
                           Method=method,
                           PredBest.Variety=temp.variety,
                           PredBest.Yield=temp[temp.variety],
                           ObsBest.Variety=names(which.max(temp)),
                           ObsBest.Yield=max(temp),
                           ObsWorst.Variety=names(which.min(temp)),
                           ObsWorst.Yield=min(temp)))
    
  }
}

site <- sort(unique(bv$Site))
bv$Year <- factor(x=bv$Year, levels=year)
bv$Site <- factor(x=bv$Site, levels=site)
bv$Method <- factor(x=bv$Method, levels=method)

# calculate the deficit (best obs - best pred).
bv$Deficit <- bv$ObsBest.Yield - bv$PredBest.Yield

# calculate the % deficit (best obs - best pred)/(best obs - worst obs).
bv$PD <- (bv$ObsBest.Yield - bv$PredBest.Yield)/(bv$ObsBest.Yield - bv$ObsWorst.Yield)

# fit a linear model for Deficit~Year+Site+Method.
lm0 <- lm(data=bv, formula=Deficit~Year+Site+Method)
anova(lm0)
#            Df  Sum Sq Mean Sq F value    Pr(>F)    
#Year        12   9.989 0.83238  7.8318 3.818e-14 ***
#Site        28  10.696 0.38201  3.5943 1.324e-09 ***
#Method       7   1.455 0.20782  1.9553   0.05809 .  
#Residuals 1096 116.486 0.10628

# get the mean Deficit for each Method.
em0 <- emmeans(lm0, "Method")
data.frame(em0)
#  Method    emmean         SE   df  lower.CL  upper.CL
#1     I1 0.4259561 0.02941174 1096 0.3682465 0.4836658
#2     I5 0.4034387 0.02941174 1096 0.3457290 0.4611483
#3    I5B 0.3681939 0.02941174 1096 0.3104842 0.4259036
#4     R1 0.3467254 0.02941174 1096 0.2890157 0.4044351
#5     R5 0.4027394 0.02941174 1096 0.3450297 0.4604490
#6     A1 0.3168652 0.02941174 1096 0.2591556 0.3745749
#7     A5 0.3567254 0.02941174 1096 0.2990157 0.4144351
#8    A5B 0.3337883 0.02941174 1096 0.2760786 0.3914980

### aggregate the data by county (hard to plot site on map).
# load the site-county information.
temp <- read.csv("input/boundary/site_to_county.csv", as.is=TRUE)
site2county <- temp$County
names(site2county) <- temp$Site

# first extract only the important columns from bv.
bv2 <- bv
bv2$County <- site2county[as.character(bv2$Site)]
bv2 <- bv2[, c(12,1,3,10,11)]

# take the average if there are multiple site within each county.
temp1 <- tapply(bv2$Deficit, bv2[, c("County", "Year", "Method")], mean)
temp2 <- tapply(bv2$PD, bv2[, c("County", "Year", "Method")], mean)
bv2 <- cbind(melt(temp1, na.rm=TRUE), melt(temp2, na.rm=TRUE))
bv2 <- bv2[, c(1,2,3,4,8)]
colnames(bv2)[4:5] <- c("Deficit", "PD")

# set the columns to characters (factors don't work well in map plotting).
bv2$County <- as.character(bv2$County)
bv2$Method <- as.character(bv2$Method)

# take the average across year for each county-method.
temp1 <- tapply(bv2$Deficit, bv2[, c("County", "Method")], mean)
temp2 <- tapply(bv2$PD, bv2[, c("County", "Method")], mean)
bv3 <- cbind(melt(temp1, na.rm=TRUE), melt(temp2, na.rm=TRUE))
bv3 <- bv3[, c(1,2,3,6)]
colnames(bv3)[3:4] <- c("Deficit", "PD")

# set the columns to characters (factors don't work well in map plotting).
bv3$County <- as.character(bv3$County)
bv3$Method <- as.character(bv3$Method)


### save the results for plotting.
save(bv, bv2, bv3,
     cor.ind1, cor.ind5, cor.ind.blup, cor.reg1, cor.reg5, cor.all1, cor.all5, cor.all.blup,
     cors.ind1, cors.ind5, cors.ind.blup, cors.reg1, cors.reg5, cors.all1, cors.all5, cors.all.blup,
     geno, pheno, out.old,
     site, method, year, region,
     site2county, first, file="input/results.RData")


# load the results to plot.
load("input/results.RData")

# load the map data.
load("input/boundary/uk.RData")

### Fig 1A. Map of number of trials in each county over the years (2002 to 2019).
# prepare the data.
temp <- table(pheno[pheno$Year < 2020, c("County", "Year")])
dim(temp) # 39 18
temp <- rowSums(temp > 0)
uk0 <- uk2
uk0$count <- temp[uk0$County]
uk0 <- uk0[, c(1,3,2)]

# plot.
p0 <- ggplot() +
  geom_sf(data=uk0, aes(fill=count)) +
  theme(panel.background=element_rect(fill="#99CCFF"), panel.grid=element_blank()) +
  theme(axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(axis.title=element_blank()) +
  theme(axis.ticks.length=unit(0, "pt")) +
  theme(legend.background=element_rect(fill="#99CCFF", color="#000000", linewidth=0.1)) +
  theme(legend.position=c(0.85, 0.85)) +
  theme(legend.key.size=unit(0.7, "lines")) +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=6)) +
  theme(plot.margin=margin(t=0, r=0, b=0, l=0)) +
  scale_fill_gradient(low="#FFFFFF", high="#0000FF", na.value="#DDDDDD", breaks=c(0,6,12,18), limits=c(0,18), guide=guide_colorbar(frame.colour="#000000", ticks.colour="#000000"))

ggsave(p0,
       filename="output/fig/1A_data_count.png",
       height=4,
       width=2.2,
       units="in",
       dpi=600)

# note: don't save maps as svg - too big.

### Fig 1B. Map of relevant weather station data.

# create  an empty map, and add the weather data layer in Inkscape.
p0 <- ggplot() +
  geom_sf(data=uk2) +
  theme(panel.background=element_rect(fill="#99CCFF"), panel.grid=element_blank()) +
  theme(axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(axis.title=element_blank()) +
  theme(axis.ticks.length=unit(0, "pt")) +
  theme(plot.margin=margin(t=0, r=0, b=0, l=0))

ggsave(p0,
       filename="output/fig/1B_weather_blank.png",
       height=4,
       width=2.2,
       units="in",
       dpi=600)

# prepare the data.
df <- read.csv("input/historical_weather.csv", as.is=TRUE)
df <- df[df$Year < 2020, ] # Year 2002 to 2019.
df$Trange_C <- df$Tmax_C - df$Tmin_C
df <- df[, c(1:6, 9, 7:8)]
df2 <- data.frame(Station=unique(df$Station),
                  Tmax=tapply(df$Tmax_C, df$Station, mean, na.rm=TRUE),
                  Tmin=tapply(df$Tmin_C, df$Station, mean, na.rm=TRUE),
                  Trange=tapply(df$Trange_C, df$Station, mean, na.rm=TRUE),
                  Frost=tapply(df$Frost_Day, df$Station, mean, na.rm=TRUE),
                  Rain=tapply(df$Rain_mm, df$Station, mean, na.rm=TRUE))
rownames(df2) <- NULL


# rotate coordinate.
rot <- function(r, t){
  x <- r*cos(90*pi/180)
  y <- r*sin(90*pi/180)
  xp <- x*cos(t*pi/180) + y*sin(t*pi/180)
  yp <- -x*sin(t*pi/180) + y*cos(t*pi/180)
  return(c(xp, yp))
}

# base pentagons for the radar plot.
df00 <- data.frame(rbind(cbind(t(sapply(0:5, FUN=function(i) rot(r=100, t=72*i))), group=1),
                        cbind(t(sapply(0:5, FUN=function(i) rot(r=75, t=72*i))), group=2),
                        cbind(t(sapply(0:5, FUN=function(i) rot(r=50, t=72*i))), group=3),
                        cbind(t(sapply(0:5, FUN=function(i) rot(r=25, t=72*i))), group=4)))
colnames(df00)[1:2] <- c("x", "y")

# base axes for the radar plot.
df01 <- data.frame(x=c(rep(0, 5), df00$x[1:5]),
                   y=c(rep(0, 5), df00$y[1:5]),
                   group=rep(1:5,2))

# label for each vertex on the radar plot.
sapply(2:6, FUN=function(i) range(df2[,i]))
# 13.83889  4.875926 5.450000 0.2222222  44.67685
# 19.15741 10.592593 9.807407 5.4768519 129.13981

df02 <- data.frame(df00[c(1:5), 1:2],
                   lab=c("Tmax\n14-19 C",
                         "Tmin\n5-11 C",
                         "Trange\n6-10 C",
                         "Frost\n0-6 day",
                         "Rain\n43-127 mm"))
df02$x <- c(df02$x[1],
            df02$x[2]+60,
            df02$x[3]+20,
            df02$x[4]-20,
            df02$x[5]-90)
df02$y <- c(df02$y[1]+40,
            df02$y[2]+10,
            df02$y[3]-40,
            df02$y[4]-40,
            df02$y[5]+10)

# prepare the data for radar plot.
df3 <- vector()
for(i in 1:nrow(df2)){
  temp <- rbind(rot(r=100*(df2$Tmax[i]-min(df2$Tmax))/diff(range(df2$Tmax)), t=72*0),
                rot(r=100*(df2$Tmin[i]-min(df2$Tmin))/diff(range(df2$Tmin)), t=72*1),
                rot(r=100*(df2$Trange[i]-min(df2$Trange))/diff(range(df2$Trange)), t=72*2),
                rot(r=100*(df2$Frost[i]-min(df2$Frost))/diff(range(df2$Frost)), t=72*3),
                rot(r=100*(df2$Rain[i]-min(df2$Rain))/diff(range(df2$Rain)), t=72*4),
                rot(r=100*(df2$Tmax[i]-min(df2$Tmax))/diff(range(df2$Tmax)), t=72*5))
  temp <- data.frame(df2$Station[i], temp)
  colnames(temp) <- c("Station", "x", "y")
  df3 <- rbind(df3, temp)
}

p0 <- ggplot() +
  geom_polygon(data=df00, aes(x=x, y=y, group=group), color="#AAAAAA", fill="#FFFFFF", linewidth=0.1) +
  geom_line(data=df01, aes(x=x, y=y, group=group), color="#999999", linewidth=0.2) +
  geom_polygon(data=df3, aes(x=x, y=y), color="#FF5555", fill="#FF5555", alpha=0.5, linewidth=0.5) +
  facet_wrap(vars(Station), ncol=2) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title=element_blank(), axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(size=6)) +
  scale_x_continuous(expand=c(0.10, 0)) +
  scale_y_continuous(expand=c(0.10, 0)) +
  coord_fixed()

ggsave(plot=p0,
       filename="output/fig/1B_weather_radar_plot.svg",
       height=8,
       width=4,
       units="in")

# create an empty radar plot with annotations.
p0 <- ggplot() +
  geom_polygon(data=df00, aes(x=x, y=y, group=group), color="#AAAAAA", fill="#FFFFFF", linewidth=0.2) +
  geom_line(data=df01, aes(x=x, y=y, group=group), color="#999999", linewidth=0.4) +
  geom_text(data=df02, aes(x=x, y=y, label=lab), size=2) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title=element_blank(), axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(size=6)) +
  scale_x_continuous(expand=c(0.50, 0)) +
  scale_y_continuous(expand=c(0.20, 0)) +
  coord_fixed()

ggsave(plot=p0,
       filename="output/fig/1B_weather_radar_empty.svg",
       height=1.1,
       width=2,
       units="in")

### Fig 1C. Scatter plot of predicted vs observed ranks.
# Example with Coaltown of Balgonie, 2007 (best year) and 2016 (worst year).
# get the yield ranks for Coaltown of Balgonie, 2007 (best year).
df1 <- data.frame(Year=2007,
                  Observed=out.old[[1]]$target[, "Coaltown of Balgonie"],
                  I1=out.old[[1]]$ind1[, "Coaltown of Balgonie"],
                  I5=out.old[[1]]$ind5[, "Coaltown of Balgonie"],
                  R1=out.old[[1]]$reg1[, "Coaltown of Balgonie"],
                  R5=out.old[[1]]$reg5[, "Coaltown of Balgonie"],
                  A1=out.old[[1]]$all1,
                  A5=out.old[[1]]$all5)
df1 <- df1[!is.na(df1$Observed), ]
for(i in 2:8) df1[,i] <- rank(df1[,i])

# get the yield ranks for Coaltown of Balgonie, 2016 (worst year).
df2 <- data.frame(Year=2016,
                  Observed=out.old[[10]]$target[, "Coaltown of Balgonie"],
                  I1=out.old[[10]]$ind1[, "Coaltown of Balgonie"],
                  I5=out.old[[10]]$ind5[, "Coaltown of Balgonie"],
                  R1=out.old[[10]]$reg1[, "Coaltown of Balgonie"],
                  R5=out.old[[10]]$reg5[, "Coaltown of Balgonie"],
                  A1=out.old[[10]]$all1,
                  A5=out.old[[10]]$all5)
df2 <- df2[!is.na(df2$Observed), ]
for(i in 2:8) df2[,i] <- rank(df2[,i])

# combine and organize the data for plotting.
df <- rbind(df1, df2)
df <- melt(df, id.vars=c("Year", "Observed"))
colnames(df)[3:4] <- c("Method", "Predicted")

# get the best prediction to annotate.
df0 <- vector()
for(i in unique(df$Year)){
  for(j in unique(df$Method)){
    temp <- df[df$Year==i & df$Method==j, ]
    temp <- temp[which.max(temp$Predicted), ]
    df0 <- rbind(df0, temp)
  }
}

# create a scatter plot of predicted vs observed ranks.
p0 <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill="#FFFFFF", color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=Inf, color="#9999FF", linewidth=0.5) +
  geom_segment(data=df0, aes(x=Observed, xend=Observed, y=-Inf, yend=Predicted), col="#FF9999", linewidth=0.5) +
  geom_point(data=df, aes(x=Observed, y=Predicted), size=0.5) +
  facet_grid(rows=vars(Method), cols=vars(Year)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=6)) +
  theme(strip.text=element_text(size=8, margin=margin(0.1,0.1,0.1,0.1,unit="line"))) +
  theme(axis.ticks=element_line(linewidth=0.2)) +
  scale_x_continuous(breaks=c(0,5,10,15), minor_breaks=0:15, guide="axis_minor") +
  scale_y_continuous(breaks=c(0,5,10,15), minor_breaks=0:15, guide="axis_minor") +
  xlab("observed rank") +
  ylab("predicted rank") +
  coord_fixed(xlim=c(0,15), ylim=c(0,15))

ggsave(p0,
       filename="output/fig/1C_ranks.svg",
       height=4,
       width=2.2,
       units="in")


### Fig 2. Boxplot of all correlations.
# prepare the data for plotting.
df <- rbind(cbind(melt(cors.ind1), Method="I1"),
            cbind(melt(cors.ind5), Method="I5"),
            cbind(melt(cors.ind.blup), Method="I5B"),
            cbind(melt(cors.reg1), Method="R1"),
            cbind(melt(cors.reg5), Method="R5"),
            cbind(melt(cors.all1), Method="A1"),
            cbind(melt(cors.all5), Method="A5"),
            cbind(melt(cors.all.blup), Method="A5B"))
colnames(df) <- c("Site", "Year", "Correlation", "Method")
df$Year <- factor(x=df$Year, levels=year)
df$Method <- factor(x=df$Method, levels=method)

round(tapply(df$Correlation, df$Method, median, na.rm=TRUE), 3)
#   I1    I5   I5B    R1    R5    A1    A5   A5B 
#0.414 0.324 0.518 0.496 0.406 0.518 0.438 0.571

# plot the heatmap.
p0 <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill="#FFFFFF", color="#DDDDDD") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=c(1.5, 3.5, 5.5, 7.5), ymax=c(2.5, 4.5, 6.5, 8.5), fill="#EEEEEE", color="#FFFFFF00") +
  annotate("segment", x=0.571, xend=0.571, y=-Inf, yend=Inf, color="#FF9999", linewidth=1) +
  geom_boxplot(data=df, aes(y=Method, x=Correlation), na.rm=TRUE, width=0.5, outlier.size=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=6)) +
  theme(axis.ticks=element_line(color="#DDDDDD")) +
  scale_y_discrete(limits=rev) +
  xlab("Correlation") +
  ylab("Method") +
  coord_cartesian(xlim=c(-1,1))

ggsave(p0,
       filename="output/fig/2_boxplot_spearman.png",
       height=3,
       width=7,
       units="in",
       dpi=600)

### Fig 3A. Heatmap of Spearman correlations for 6 cities (Y: Method, X: Year, Panel: Site).
# prepare the data for plotting.
df <- rbind(cbind(melt(cors.ind1), Method="I1"),
            cbind(melt(cors.ind5), Method="I5"),
            cbind(melt(cors.ind.blup), Method="I5B"),
            cbind(melt(cors.reg1), Method="R1"),
            cbind(melt(cors.reg5), Method="R5"),
            cbind(melt(cors.all1), Method="A1"),
            cbind(melt(cors.all5), Method="A5"),
            cbind(melt(cors.all.blup), Method="A5B"))
colnames(df) <- c("Site", "Year", "Correlation", "Method")
df$Year <- factor(x=df$Year, levels=year)
df$Method <- factor(x=df$Method, levels=method)

# filter the data to 6 cities with the most data.
df <- df[df$Site%in%c("Coaltown of Balgonie", "Laurencekirk", "St Boswells", "Stockbridge", "Tain", "Wymondham"), ]
df$Site <- droplevels(df$Site)

# plot the heatmap.
p0 <- ggplot() +
  geom_tile(data=df, aes(x=Year, y=Method, fill=Correlation), color="#555555") +
  facet_wrap(vars(Site), ncol=3, dir="h") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.ticks=element_line(color="#555555")) +
  theme(axis.text.x=element_text(angle=40, hjust=1, vjust=1, size=6)) +
  theme(axis.text.y=element_text(size=6)) +
  theme(axis.title=element_text(size=8)) +
  theme(axis.ticks=element_line(color="#555555", linewidth=0.5)) +
  theme(legend.title=element_text(size=8, vjust=2), legend.text=element_text(size=6)) +
  theme(legend.key.width=unit(0.75, "line"), legend.key.height=unit(0.75, "line")) +
  theme(legend.text.align=1) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  scale_fill_gradient2(low="#FF0000", mid="#FFFFFF", high="#0000FF", limits=c(-1,1), breaks=c(-1,-0.5,0,0.5,1), name="Correlation", na.value="#AAAAAA") +
  scale_y_discrete(limits=rev, expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  coord_fixed() +
  guides(fill=guide_colorbar(ticks.colour="#555555", ticks.linewidth=0.2, frame.colour="#555555", frame.linewidth=0.2, title.position="top", title.hjust=0.5))

ggsave(p0,
       filename="output/fig/3A_heatmap_subset.png",
       height=3.5,
       width=7,
       units="in",
       dpi=600)

ggsave(p0,
       filename="output/fig/3A_heatmap_subset.svg",
       height=3.5,
       width=7,
       units="in")

# some numbers for results section.
round(tapply(df$Correlation, df[, c("Site", "Method")], mean, na.rm=TRUE), 3)
#                      Method
#Site                      I1    I5   I5B    R1    R5    A1    A5   A5B
#  Coaltown of Balgonie 0.253 0.269 0.407 0.301 0.299 0.353 0.301 0.421
#  Laurencekirk         0.517 0.323 0.632 0.592 0.460 0.578 0.506 0.643
#  St Boswells          0.292 0.318 0.412 0.441 0.487 0.392 0.417 0.485
#  Stockbridge          0.382 0.205 0.514 0.363 0.318 0.482 0.326 0.563
#  Tain                 0.571 0.314 0.538 0.551 0.345 0.572 0.322 0.524
#  Wymondham            0.300 0.259 0.420 0.410 0.322 0.431 0.375 0.459

round(rowSums(tapply(df$Correlation, df[, c("Site", "Method")], mean, na.rm=TRUE))/6, 3)
#Coaltown of Balgonie         Laurencekirk          St Boswells 
#               0.434                0.708                0.541 
#         Stockbridge                 Tain            Wymondham 
#               0.526                0.623                0.496 

### Fig 3B. Jitter and mean plot of spearman correlations for each method by Site.
# get the relevant data.
df <- rbind(cbind(melt(cors.ind1), Method="I1"),
            cbind(melt(cors.ind5), Method="I5"),
            cbind(melt(cors.ind.blup), Method="I5B"),
            cbind(melt(cors.reg1), Method="R1"),
            cbind(melt(cors.reg5), Method="R5"),
            cbind(melt(cors.all1), Method="A1"),
            cbind(melt(cors.all5), Method="A5"),
            cbind(melt(cors.all.blup), Method="A5B"))
colnames(df) <- c("Site", "Year", "Correlation", "Method")
df$Year <- factor(x=df$Year, levels=year)
df$Method <- factor(x=df$Method, levels=method)

# filter the data to 6 cities with the most data.
df <- df[df$Site%in%c("Coaltown of Balgonie", "Laurencekirk", "St Boswells", "Stockbridge", "Tain", "Wymondham"), ]
df$Site <- droplevels(df$Site)

# t.test to compare the correlation means across methods.
pval <- vector()
for(i in levels(df$Site)){
  temp <- vector()
  for(j in 1:7){
    for(k in (j+1):8){
      temp <- c(temp,
                t.test(x=df$Correlation[df$Site==i & df$Method==method[j]],
                       y=df$Correlation[df$Site==i & df$Method==method[k]],
                       paired=TRUE)$p.value)
    }
  }
  pval <- cbind(pval, temp)
}
write.csv(pval, "output/pval_site.csv", row.names=FALSE, quote=FALSE)

# calculate means.
df2 <- tapply(df$Correlation, df[, c("Site", "Method")], mean, na.rm=TRUE)
df2 <- melt(df2)
colnames(df2) <- c("Site", "Method", "Correlation")

# create the plot.
p0 <- ggplot() +
  annotate("tile", y=factor(x=method[c(1,3,5,7)], levels=method), x=0.15, fill="#EEEEEE", color=NA, width=2, height=1) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#BBBBBB") +
  geom_jitter(data=df, aes(y=Method, x=Correlation), size=0.5, height=0.1, width=0, na.rm=TRUE) +
  geom_point(data=df2, aes(y=Method, x=Correlation), size=1, color="#FF9999") +
  facet_wrap(vars(Site), ncol=3, dir="h") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  theme(axis.title=element_text(size=8)) +
  theme(axis.text=element_text(size=6)) +
  theme(axis.ticks=element_line(color="#BBBBBB")) +
  scale_y_discrete(drop=FALSE, limits=rev) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Method") +
  xlab("Correlation")

ggsave(p0,
       filename="output/fig/3B_site_spearman.png",
       height=3.5,
       width=6,
       units="in",
       dpi=600)

ggsave(p0,
       filename="output/fig/3B_site_spearman.svg",
       height=3.5,
       width=6,
       units="in")


### Fig S3. Heatmap of Spearman correlations for all cities (Y: Method, X: Year, Panel: Site).
# prepare the data for plotting.
df <- rbind(cbind(melt(cors.ind1), Method="I1"),
            cbind(melt(cors.ind5), Method="I5"),
            cbind(melt(cors.ind.blup), Method="I5B"),
            cbind(melt(cors.reg1), Method="R1"),
            cbind(melt(cors.reg5), Method="R5"),
            cbind(melt(cors.all1), Method="A1"),
            cbind(melt(cors.all5), Method="A5"),
            cbind(melt(cors.all.blup), Method="A5B"))
colnames(df) <- c("Site", "Year", "Correlation", "Method")
df$Year <- factor(x=df$Year, levels=year)
df$Method <- factor(x=df$Method, levels=method)

# plot the heatmap.
p0 <- ggplot() +
  geom_tile(data=df, aes(x=Year, y=Method, fill=Correlation), color="#555555") +
  facet_wrap(vars(Site), ncol=5, dir="v") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.ticks=element_line(color="#555555")) +
  theme(axis.text.x=element_text(angle=40, hjust=1, vjust=1, size=5)) +
  theme(axis.text.y=element_text(size=7)) +
  theme(axis.title=element_text(size=8)) +
  theme(legend.title=element_text(size=8, vjust=2), legend.text=element_text(size=6)) +
  theme(legend.key.width=unit(0.75, "line"), legend.key.height=unit(0.75, "line")) +
  theme(legend.text.align=1) +
  theme(legend.margin=margin(0,0,0,0)) +
  theme(legend.position=c(0.925,0.045)) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  scale_fill_gradient2(low="#FF0000", mid="#FFFFFF", high="#0000FF", limits=c(-1,1), breaks=c(-1,-0.5,0,0.5,1), name="Correlation", na.value="#AAAAAA") +
  scale_y_discrete(limits=rev, expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  coord_fixed() +
  guides(fill=guide_colorbar(ticks.colour="#555555", ticks.linewidth=0.5, frame.colour="#555555", frame.linewidth=0.5, title.position="top", title.hjust=0.5))

ggsave(p0,
       filename="output/fig/S2_heatmap_full.png",
       height=7,
       width=7,
       units="in",
       dpi=600)


### Fig 4A. Jitter and mean plot of spearman correlations for each method by Year.
# get the relevant data.
df <- rbind(cbind(melt(cors.ind1), Method="I1"),
            cbind(melt(cors.ind5), Method="I5"),
            cbind(melt(cors.ind.blup), Method="I5B"),
            cbind(melt(cors.reg1), Method="R1"),
            cbind(melt(cors.reg5), Method="R5"),
            cbind(melt(cors.all1), Method="A1"),
            cbind(melt(cors.all5), Method="A5"),
            cbind(melt(cors.all.blup), Method="A5B"))
colnames(df) <- c("Site", "Year", "Correlation", "Method")
df$Year <- factor(x=df$Year, levels=year)
df$Method <- factor(x=df$Method, levels=method)

# t.test to compare the correlation means across methods.
pval <- vector()
for(i in levels(df$Year)){
  temp <- vector()
  for(j in 1:7){
    for(k in (j+1):8){
      temp <- c(temp,
                t.test(x=df$Correlation[df$Year==i & df$Method==method[j]],
                       y=df$Correlation[df$Year==i & df$Method==method[k]],
                       paired=TRUE)$p.value)
    }
  }
  pval <- cbind(pval, temp)
}
colnames(pval) <- levels(df$Year)
write.csv(pval, "output/pval_year.csv", row.names=FALSE, quote=FALSE)

# prepare the data for plotting (calculate mean and sd).
df2 <- tapply(df$Correlation, df[, c("Year", "Method")], mean, na.rm=TRUE)
df2 <- melt(df2)
colnames(df2) <- c("Year", "Method", "Correlation")
df2$Year <- factor(x=df2$Year, levels=year)

# create the plot.
p0 <- ggplot() +
  annotate("tile", y=factor(x=method[c(1,3,5,7)], levels=method), x=0.15, fill="#EEEEEE", color=NA, width=2, height=1) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#BBBBBB") +
  geom_jitter(data=df, aes(y=Method, x=Correlation), size=0.25, height=0.1, width=0, na.rm=TRUE) +
  geom_point(data=df2, aes(y=Method, x=Correlation), size=1, color="#FF9999") +
  facet_wrap(vars(Year), ncol=4, dir="h") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  theme(axis.title=element_text(size=8)) +
  theme(axis.text=element_text(size=6)) +
  theme(axis.ticks=element_line(color="#BBBBBB")) +
  scale_y_discrete(drop=FALSE, limits=rev) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#BBBBBB", "#FFBBBB")) +
  ylab("Method") +
  xlab("Correlation")

ggsave(p0,
       filename="output/fig/4A_year_spearman.png",
       height=6.125,
       width=7,
       units="in",
       dpi=600)

ggsave(p0,
       filename="output/fig/4A_year_spearman.svg",
       height=6.125,
       width=7,
       units="in")


### Fig 4B. Line plot of mean spearman correlations for each method by Year.
# get the relevant data.
df <- rbind(cbind(melt(cors.ind1), Method="I1"),
            cbind(melt(cors.ind5), Method="I5"),
            cbind(melt(cors.ind.blup), Method="I5B"),
            cbind(melt(cors.reg1), Method="R1"),
            cbind(melt(cors.reg5), Method="R5"),
            cbind(melt(cors.all1), Method="A1"),
            cbind(melt(cors.all5), Method="A5"),
            cbind(melt(cors.all.blup), Method="A5B"))
colnames(df) <- c("Site", "Year", "Correlation", "Method")
df$Year <- factor(x=df$Year, levels=year)
df$Method <- factor(x=df$Method, levels=method)

# compute the mean across sites for each year/method.
df <- tapply(df$Correlation, df[, c("Method", "Year")], mean, na.rm=TRUE)
df <- melt(df)
df$Year <- factor(df$Year, levels=year)
colnames(df) <- c("Method", "Year", "Correlation")

df$Method1 <- NA
df$Method2 <- NA
for(i in 1:nrow(df)){
  
  if(df$Method[i]%in%c("I1", "I5", "I5B")){
    df$Method1[i] <- "I"
  } else if(df$Method[i]%in%c("R1", "R5")){
    df$Method1[i] <- "R"
  } else if(df$Method[i]%in%c("A1", "A5", "A5B")){
    df$Method1[i] <- "A"
  }
  
  if(df$Method[i]%in%c("I1", "R1", "A1")){
    df$Method2[i] <- "1"
  } else if(df$Method[i]%in%c("I5", "R5", "A5")){
    df$Method2[i] <- "5"
  } else if(df$Method[i]%in%c("I5B", "A5B")){
    df$Method2[i] <- "5B"
  }
  
}
df$Method1 <- factor(df$Method1, levels=c("I", "R", "A"))
df$Method2 <- factor(df$Method2, levels=c("1", "5", "5B"))

# create the plot.
p0 <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#BBBBBB") +
  geom_line(data=df, aes(x=Year, y=Correlation, color=Method1, linetype=Method2, group=Method), linewidth=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.key=element_blank()) +
  theme(legend.key.height=unit(0.6, "line"), legend.key.width=unit(2, "lines")) +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=6)) +
  theme(axis.title=element_text(size=8)) +
  theme(axis.text=element_text(size=6)) +
  theme(axis.ticks=element_line(color="#BBBBBB")) +
  scale_color_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62")) +
  scale_linetype_manual(values=c(3,2,1))

ggsave(p0,
       filename="output/fig/4B_year_spearman.png",
       height=1.6,
       width=5,
       units="in",
       dpi=600)

ggsave(p0,
       filename="output/fig/4B_year_spearman.svg",
       height=1.6,
       width=5,
       units="in")


### Fig 5. Map showing % deficit in yield for each method (averaged across years).
# loop and plot the % deficit for each method.
# manually combine the 8 sub-figures.
for(j in method){
  
  # get the data for method-j.
  temp <- bv3[bv3$Method==j, c("County", "PD")]
  PD <- temp$PD*100
  names(PD) <- temp$County
  
  # add % deficit to uk2.
  uk0 <- uk2
  uk0$PD <- PD[uk0$County]
  uk0 <- uk0[, c(1,3,2)]
  
  # plot.
  p0 <- ggplot() +
    geom_sf(data=uk0, aes(fill=PD)) +
    theme(panel.background=element_rect(fill="#99CCFF"), panel.grid=element_blank()) +
    theme(axis.text=element_blank(), axis.ticks=element_blank()) +
    theme(axis.title=element_blank()) +
    theme(axis.ticks.length=unit(0, "pt")) +
    theme(legend.position="none") +
    theme(plot.margin=margin(t=0, r=0, b=0, l=0)) +
    scale_fill_gradient(low="#FFFFFF", high="#FF0000", na.value="#DDDDDD", limits=c(0,100))
  
  ggsave(p0,
         filename=paste("output/fig/map/PD_", j, ".png", sep=""),
         height=3,
         width=1.625,
         units="in",
         dpi=600)
  
}

# create a dummy plot to get the legend.
p0 <- ggplot() +
  geom_col(aes(x=c(1,2), y=c(1,2), fill=c(0,100))) +
  scale_fill_gradient(low="#FFFFFF", high="#FF0000", na.value="#DDDDDD", limits=c(0,100), name="PD",
                      guide=guide_colorbar(frame.colour="#000000", ticks.colour="#000000"))
ggsave(p0,
       filename="output/fig/map/PD_legend.svg",
       height=3,
       width=3,
       units="in")


### Fig 6. Cumulative PD
# for each Method-Year, calculate mean in PD (averaged across all trial sites).
df <- tapply(100*bv$PD, bv[, c("Method", "Year")], mean)
df <- melt(df)
colnames(df) <- c("Method", "Year", "PD")
df <- lapply(method, FUN=function(i) df[df$Method==i, ])
for(i in 1:length(df)) df[[i]]$CPD <- cumsum(df[[i]]$PD)
df <- do.call(rbind, df)
df$Year <- factor(df$Year, levels=year)

df$Method1 <- NA
df$Method2 <- NA
for(i in 1:nrow(df)){
  
  if(df$Method[i]%in%c("I1", "I5", "I5B")){
    df$Method1[i] <- "I"
  } else if(df$Method[i]%in%c("R1", "R5")){
    df$Method1[i] <- "R"
  } else if(df$Method[i]%in%c("A1", "A5", "A5B")){
    df$Method1[i] <- "A"
  }
  
  if(df$Method[i]%in%c("I1", "R1", "A1")){
    df$Method2[i] <- "1"
  } else if(df$Method[i]%in%c("I5", "R5", "A5")){
    df$Method2[i] <- "5"
  } else if(df$Method[i]%in%c("I5B", "A5B")){
    df$Method2[i] <- "5B"
  }
  
}
df$Method1 <- factor(df$Method1, levels=c("I", "R", "A"))
df$Method2 <- factor(df$Method2, levels=c("1", "5", "5B"))

# create the plot.
p0 <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#BBBBBB") +
  geom_line(data=df, aes(x=Year, y=CPD, color=Method1, linetype=Method2, group=Method), linewidth=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.key=element_blank()) +
  theme(legend.key.height=unit(0.6, "line"), legend.key.width=unit(2, "lines")) +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=6)) +
  theme(axis.title=element_text(size=8)) +
  theme(axis.text=element_text(size=6)) +
  theme(axis.ticks=element_line(color="#BBBBBB")) +
  scale_color_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62")) +
  scale_linetype_manual(values=c(3,2,1))

ggsave(p0,
       filename="output/fig/6_cumulative_PD.png",
       width=7,
       height=3,
       units="in",
       dpi=600)

ggsave(p0,
       filename="output/fig/6_cumulative_PD.svg",
       width=7,
       height=3,
       units="in",
       dpi=600)







### Fig S1. Map of number of trials in each county over the years (2002 to 2019, similar to 1A).
# prepare the data.
temp <- table(pheno[pheno$Year < 2020, c("County", "Year")])
dim(temp) # 39 18
temp <- rowSums(temp > 0)
uk0 <- uk2
uk0$count <- temp[uk0$County]
uk0 <- uk0[, c(1,3,2)]

# plot.
p0 <- ggplot() +
  geom_sf(data=uk0, aes(fill=count)) +
  theme(panel.background=element_rect(fill="#99CCFF"), panel.grid=element_blank()) +
  theme(axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(axis.title=element_blank()) +
  theme(axis.ticks.length=unit(0, "pt")) +
  theme(legend.background=element_rect(fill="#99CCFF", color="#000000", linewidth=0.1)) +
  theme(legend.position=c(0.10, 0.90)) +
  theme(legend.key.size=unit(1, "lines")) +
  theme(legend.title=element_text(size=10), legend.text=element_text(size=8)) +
  theme(legend.spacing.y=unit(0.8, "lines")) +
  theme(plot.margin=margin(t=0, r=0, b=0, l=0)) +
  scale_fill_gradient(low="#FFFFFF", high="#0000FF", na.value="#DDDDDD", breaks=c(0,6,12,18), limits=c(0,18), guide=guide_colorbar(frame.colour="#000000", ticks.colour="#000000"))

ggsave(p0,
       filename="output/fig/S1_data_count.png",
       height=9,
       width=4.95,
       units="in",
       dpi=600)

# note: don't save maps as svg - too big.


