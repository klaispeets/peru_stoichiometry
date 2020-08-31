setwd("/Users/riina82/work/OTHERS/Kristian/peru")
library(mgcv)
library(vegan)
rm(list=ls())

load("peru_short.RData")
data = peru

#PON misses the values on 8th day, DON & DOP on the 28th day. Use each mesocosms 10 day value for missing PON, and mean of 26 & 30 day for DON & DOP.
mesoc = unique(data$MesoC)
for(i in 1:length(mesoc)){
  idx = which(data$DAY == 8 & data$MesoC == mesoc[i])
  idx2 =  which(data$DAY == 10 & data$MesoC == mesoc[i])
  data$PON[idx] = data$PON[idx2] 
  
  idx = which(data$DAY == 28 & data$MesoC == mesoc[i])
  idx2 = which(data$DAY %in% c(26,30) & data$MesoC == mesoc[i])
  
  data$DON[idx] = mean(data$DON[idx2])
  data$DOP[idx] = mean(data$DOP[idx2])
  
}


#Select the biogeochemical variables
set1 = subset(data, select = c("DIN", "PO4", "Dsi", "Chla", "POC", "PON", "POP", "Bsi", "DON", "DOP"))

#Select the phytoplankton variables
set2 = subset(data, select= c("FL4.group", "Synechococcus", "Cryptophytes","Chains","Picoeukaryotes","Nano","Mikro.I", "Mikro.II"   ))

library(vegan)

mds_env = metaMDS(set1, trymax = 100)
mds_phy = metaMDS(set2, trymax = 100)

env_mds = data.frame(mds_env$points)
phy_mds = data.frame(mds_phy$points)
names(env_mds) = c("Emds1", "Emds2")
names(phy_mds) = c("Pmds1", "Pmds2")

result = data.frame(subset(data, select = c("MesoC", "DAY", "Treatment", "APA", "LAP")), env_mds, phy_mds)

#Analyses
library(mgcv)
#APA
gam_APA= gam(APA ~ s(Emds1, k = 4) + s(Emds2, k = 4) + s(Pmds1, k = 4) + s(Pmds2, k = 4), data = result)
summary(gam_APA)#73

#Figure 1
par(mfrow=c(2,2), mar = c(2.5,2.5,1,1), mgp = c(1.5,0.5, 0))
library(colorRamps)
library(plotrix)
col=matlab.like(16)[cut(as.numeric(result$DAY),16)]
col.text <- seq(min(result$DAY),max(result$DAY),length.out=5)
col.text <- round(col.text,digits=0)
xl = quantile(seq(min(result$Emds1), max(result$Emds1), length.out = 10),0.9)
yb = quantile(seq(min(result$APA), max(result$APA), length.out = 10),0.4)
xr = quantile(seq(min(result$Emds1), max(result$Emds1), length.out = 10),1)
yt = quantile(seq(min(result$APA), max(result$APA), length.out = 10),0.95)

for(i in 1:4){
newdata = result[,c(6:9)]
set = c(1:4)[-i]
for(j in 1:3){newdata[,set[j]] = mean(newdata[,set[j]])}   
pred = predict(gam_APA, newdata = newdata, se =T)
plot(result$APA ~ newdata[,i], ylab = "", xlab = "", main = "", col = col, pch = 16)
if(i==1){color.legend(xl,yb,xr,yt, col.text,matlab.like(16),align="lt",gradient="y", cex=0.7)}
idx = order(newdata[,i])
lines(pred$fit[idx] ~ newdata[idx,i],col = 1)
lines(pred$fit[idx] + 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
lines(pred$fit[idx] - 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
}

#Figure 2
gam_LAP= gam(LAP ~ s(Emds1, k = 4) + s(Emds2, k = 4) + s(Pmds1, k = 4) + s(Pmds2, k = 4), data = result)

yb = quantile(seq(min(result$LAP), max(result$LAP), length.out = 10),0.4)
yt = quantile(seq(min(result$LAP), max(result$LAP), length.out = 10),0.95)
xr = quantile(seq(min(result$Emds2), max(result$Emds2), length.out = 10),1)
xl = quantile(seq(min(result$Emds2), max(result$Emds2), length.out = 10),0.9)

par(mfrow=c(2,2), mar = c(2.5,2.5,1,1), mgp = c(1.5,0.5, 0))
library(colorRamps)
library(plotrix)
for(i in 1:4){
  newdata = result[,c(6:9)]
  set = c(1:4)[-i]
  for(j in 1:3){
    newdata[,set[j]] = mean(newdata[,set[j]])
  }   
  pred = predict(gam_LAP, newdata = newdata, se =T)
  plot(result$LAP ~ newdata[,i], ylab = "", xlab = "", main = "", col = col, pch = 16)
  if(i==2){color.legend(xl,yb,xr,yt, col.text,matlab.like(16),align="lt",gradient="y", cex=0.7)}
   idx = order(newdata[,i])
  lines(pred$fit[idx] ~ newdata[idx,i],col = 1)
  lines(pred$fit[idx] + 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
  lines(pred$fit[idx] - 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
}



summary(gam(APA ~ s(Emds1, k = 4), data = result))
summary(gam(APA ~ s(Emds2, k = 4), data = result))
summary(gam(APA ~ s(Pmds1, k = 4), data = result))
summary(gam(APA ~ s(Pmds2, k = 4), data = result))

summary(gam(LAP ~ s(Emds1, k = 4), data = result))
summary(gam(LAP ~ s(Emds2, k = 4), data = result))
summary(gam(LAP ~ s(Pmds1, k = 4), data = result))
summary(gam(LAP ~ s(Pmds2, k = 4), data = result))

summary(gam(APA ~ s(Emds1, k = 4) + s(Emds2, k = 4), data = result))
summary(gam(APA ~ s(Pmds1, k = 4) + s(Pmds2, k = 4), data = result))
summary(gam(LAP ~ s(Emds1, k = 4) + s(Emds2, k = 4), data = result))
summary(gam(LAP ~ s(Pmds1, k = 4) + s(Pmds2, k = 4), data = result))




summary(gam_LAP) # R2 0.249, mostly linear.

summary(gam(LAP ~ s(Emds1, k = 4), data = result))#R2 adj: 0.175
summary(gam(LAP ~ s(Emds2, k = 4), data = result))#0.134
summary(gam(LAP ~ s(Pmds1, k = 4), data = result))#0.0645
summary(gam(LAP ~ s(Pmds2, k = 4), data = result))#0.08

par(mfrow=c(1,1),mar = c(2.5,2.5,1,1), mgp = c(1.5,0.5, 0), cex.axis = 0.9)
sp = data.frame(mds_env$species)
p = data.frame(mds_env$points)
plot(MDS1~ MDS2, data = p, col = col, pch = 16, xlab = "", ylab = "", main = "")
par(new=TRUE)
plot(MDS1~ MDS2, data = sp, axes = F, col = 2, type = "n", ylab = "", xlab ="")
text(sp$MDS2, sp$MDS1, labels = row.names(sp), cex = 0.8, col =2)
col=matlab.like(16)[cut(as.numeric(result$DAY),16)]
col.text <- seq(min(result$DAY),max(result$DAY),length.out=5)
col.text <- round(col.text,digits=0)
color.legend(-0.1936824,-0.05306023, -0.1816683,0.07814014,   col.text,matlab.like(16),align="lt",gradient="y")


sp = data.frame(mds_phy$species)
p = data.frame(mds_phy$points)
plot(MDS1~ MDS2, data = p, col =col, pch = 16, xlab = "", ylab = "", main = "")
par(new=TRUE)

#Widen xlim
x1 = min(sp$MDS2)
x2 = max(sp$MDS2)
step = (x2-x1)/15
x1 = x1-step
x2 = x2+step

plot(MDS1~ MDS2, data = sp, axes = F, col = 2, type = "n", ylab = "", xlab ="", xlim = c(x1,x2))
text(sp$MDS2, sp$MDS1, labels = row.names(sp), cex = 0.8, col =2)
pos = locator(2)
color.legend(pos$x[1], pos$y[1], pos$x[2], pos$y[2],  col.text,matlab.like(16),align="lt",gradient="y")


load("peru_data_with_bacteria.RData")
bac2 = bac[,c(6:126)]
mds_bac = metaMDS(bac2, trymax = 100)

sp = data.frame(mds_bac$species)
p = data.frame(mds_bac$points)
plot(MDS1~ MDS2, data = p, col =col, pch = 16, xlab = "", ylab = "", main = "")

par(new=TRUE)
x1 = min(sp$MDS2)
x2 = max(sp$MDS2)
step = (x2-x1)/15
x1 = x1-step
x2 = x2+step
plot(MDS1~ MDS2, data = sp, axes = F, col = 2, type = "n", ylab = "", xlab ="", xlim = c(x1,x2))
text(sp$MDS2, sp$MDS1, labels = row.names(sp), cex = 0.6, col =2)

pos = locator(2)
color.legend(pos$x[1], pos$y[1], pos$x[2], pos$y[2],  col.text,matlab.like(16),align="lt",gradient="y")


results = data.frame(bac[,c(1:5)],p)
gam_LAP = gam(LAP ~ s(MDS1, k = 4) + s(MDS2, k = 4), data = results)
gam_APA = gam(APA ~ s(MDS1, k = 4) + s(MDS2, k = 4), data = results)
summary(gam_LAP)

yb = quantile(seq(min(results$APA), max(results$APA), length.out = 10),0.4)
yt = quantile(seq(min(result$APA), max(result$APA), length.out = 10),0.95)
xr = quantile(seq(min(results$MDS2), max(results$MDS2), length.out = 10),1)
xl = quantile(seq(min(results$MDS2), max(results$MDS2), length.out = 10),0.9)


par(mfrow=c(2,2), mar = c(2.5,2.5,2,1), mgp = c(1.5,0.5, 0))
for(i in 1:2){
  newdata = results[,c(6:7)]
  if(i==1){newdata$MDS2 = mean(newdata$MDS2)}else{newdata$MDS1 = mean(newdata$MDS1)}
  pred = predict(gam_APA, newdata = newdata, se =T)
  plot(results$APA ~ newdata[,i], ylab = "", xlab = "", main = "", col = col, pch = 16)
  if(i==2){color.legend(xl,yb,xr,yt, col.text,matlab.like(16),align="lt",gradient="y", cex=0.7)}
  idx = order(newdata[,i])
  lines(pred$fit[idx] ~ newdata[idx,i],col = 1)
  lines(pred$fit[idx] + 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
  lines(pred$fit[idx] - 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
}

for(i in 1:2){
  newdata = results[,c(6:7)]
  if(i==1){newdata$MDS2 = mean(newdata$MDS2)}else{newdata$MDS1 = mean(newdata$MDS1)}
  pred = predict(gam_LAP, newdata = newdata, se =T)
  plot(results$LAP ~ newdata[,i], ylab = "", xlab = "", main = "", col = col, pch = 16)
  idx = order(newdata[,i])
  lines(pred$fit[idx] ~ newdata[idx,i],col = 1)
  lines(pred$fit[idx] + 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
  lines(pred$fit[idx] - 2*pred$se.fit[idx] ~ newdata[idx,i],col = 1, lty = 2)
}

summary(gam(APA ~ s(MDS1, k = 4), data = results))# R=0.47
summary(gam(APA ~ s(MDS2, k = 4), data = results))# R=0.311
summary(gam(LAP ~ s(MDS1, k = 4), data = results))# R=0.02
summary(gam(LAP ~ s(MDS2, k = 4), data = results))# R=0.31

summary(gam(APA ~ s(MDS1, k = 4) + s(MDS2, k = 4), data = results))
summary(gam(LAP ~ s(MDS1, k = 4) + s(MDS2, k = 4), data = results))


plot(gam_APA)
plot(gam_LAP)

library(clipr)
write_clip(results)
