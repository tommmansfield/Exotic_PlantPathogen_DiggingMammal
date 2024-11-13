# WORKING DIRECTORY  ----
setwd("C:/Users/tom_m/Desktop/Uni/PhD_Stats/Dig_veg_data_sheets")
setwd("F:/PhD_Murdoch/PhD_Stats/Dig_veg_data_sheets")

# ~ ----


# Digging density estimates::DISTANCE ----

library(dplyr); library(Distance); library(knitr); library(beepr)


raw.digs.df = read.csv("combined.digs.csv", header = T,sep = ",", stringsAsFactors = F)
raw.digs.df$distance <- as.numeric(raw.digs.df$distance)


cf <- convert_units("meter", "kilometer", "square meter")
hist(raw.digs.df$distance, xlab="Distance (m)", main="Digging line transects")



# Testing for best distance models to use ----
## Testing: Zone ----
rd.zone.df <- subset (raw.digs.df, select = -Site); names(rd.zone.df)[names(rd.zone.df)=="Zone"] <- "Region.Label"; rd.zone.df[1] <- NULL 
split_rd <- split(rd.zone.df, f = rd.zone.df$Region.Label)  
split_rd

rd_zone.df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(rd_zone.df) <- c('Sample.ID','Estimate','se','cv','lcl','ucl','df')  


# for-loop over rows
for(i in 1:length(split_rd)) {       
  row <- split_rd[[i]]
  print(i)
  #print(row)
  print(row[1,4])
  if (nrow(row) > 1) {
    tryCatch({ 
      rd.hn <- ds(data=row, key="hn", adjustment=NULL, convert_units=cf)
      rd.hr <- ds(data=row, key="hr", adjustment=NULL, convert_units=cf)
      mod <- summarize_ds_models(rd.hn, rd.hr)
      print(mod)
      str(mod)
      DeltaAIC <- mod[,c('$\\Delta$AIC')]
      print(DeltaAIC)
      print(gof_ds(rd.hn))
      print(gof_ds(rd.hr))
      indexmin <- which.min(DeltaAIC)
      bestmod <- mod[,c('Key function')][indexmin]
      print(indexmin)
      print(bestmod)
      if (bestmod == "Hazard-rate") {
        print(rd.hr$dht$individuals$D)
        new_row <- rd.hr$dht$individuals$D
        new_row[new_row == 'Total'] <- row[1,4]
        new_row$model <- bestmod
        rd_zone.df <- rbind(rd_zone.df, new_row)
      } else {
        print(rd.hn$dht$individuals$D)
        new_row <- rd.hn$dht$individuals$D
        new_row[new_row == 'Total'] <- row[1,4]
        new_row$model <- bestmod
        rd_zone.df <- rbind(rd_zone.df, new_row)
      }
    }, error=function(e) {
      print('Error')
    })
  }
}; beep()

summary(rd_zone.df)
View(rd_zone.df)

# Summary:
# Infested selected for Half-normal with lower AIC, however QQplot & Cramer-von Mises test resulted in p > 0.01 so we reject and use Hazard-rate
# Both zones select for Hazard-rate model


# Now just checking Site * Zone to see if that is influential


raw.digs.df = read.csv("combined.digs.csv", header = T,sep = ",", stringsAsFactors = F)
raw.digs.df$distance <- as.numeric(raw.digs.df$distance)


cf <- convert_units("meter", "kilometer", "square meter")
hist(raw.digs.df$distance, xlab="Distance (m)", main="Digging line transects")




# Testing for best distance models to use ----
## Testing: Zone * Site ----
rd.zone.site.df <- raw.digs.df
rd.zone.site.df$zone.site <- paste0(rd.zone.site.df$Zone,rd.zone.site.df$Site)
rd.zone.site.df[, c('Zone','Site')] <- list(NULL); names(rd.zone.site.df)[names(rd.zone.site.df)=="zone.site"] <- "Region.Label"; rd.zone.site.df[1] <- NULL 

split_rd.zone.site.df <- split(rd.zone.site.df, f = rd.zone.site.df$Region.Label)  
split_rd.zone.site.df

zone.site.aic.df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(zone.site.aic.df) <- c('Sample.ID','Estimate','se','cv','lcl','ucl','df')  


# for-loop over rows
for(i in 1:length(split_rd.zone.site.df)) {       
  row <- split_rd.zone.site.df[[i]]
  print(i)
  #print(row)
  print(row[1,5])
  if (nrow(row) > 1) {
    tryCatch({ 
      rd.hn <- ds(data=row, key="hn", adjustment=NULL, convert_units=cf)
      rd.hr <- ds(data=row, key="hr", adjustment=NULL, convert_units=cf)
      mod <- summarize_ds_models(rd.hn, rd.hr)
      print(mod)
      str(mod)
      DeltaAIC <- mod[,c('$\\Delta$AIC')]
      print(DeltaAIC)
      print(gof_ds(rd.hn))
      print(gof_ds(rd.hr))
      indexmin <- which.min(DeltaAIC)
      bestmod <- mod[,c('Key function')][indexmin]
      print(indexmin)
      print(bestmod)
      if (bestmod == "Hazard-rate") {
        print(rd.hr$dht$individuals$D)
        new_row <- rd.hr$dht$individuals$D
        new_row[new_row == 'Total'] <- row[1,4]
        new_row$model <- bestmod
        zone.site.aic.df <- rbind(zone.site.aic.df, new_row)
      } else {
        print(rd.hn$dht$individuals$D)
        new_row <- rd.hn$dht$individuals$D
        new_row[new_row == 'Total'] <- row[1,4]
        new_row$model <- bestmod
        zone.site.aic.df <- rbind(zone.site.aic.df, new_row)
      }
    }, error=function(e) {
      print('Error')
    })
  }
}; beep()



summary(zone.site.aic.df)
View(zone.site.aic.df)


# Summary:
# Hazard-rate best for Alps, Black Cockatoo***, Falls, Superblock
# Half-normal selected for quail
# ***Black Cockatoo Infested: AIC selected for Half-normal, but was very close 
#                             Goodness of fit test says that both are appropriate, and better t-stat for Hazard-rate.
#                             Therefore, going with Hazard-rate since can't compare models for zones



# ~ ----
# Generating actual digging density estimates ----

raw.digs.df = read.csv("combined.digs.csv", header = T,sep = ",", stringsAsFactors = F)
raw.digs.df$distance <- as.numeric(raw.digs.df$distance)

cf <- convert_units("meter", "kilometer", "square meter")

nest.model1 <- ds(raw.digs.df, key = "hr", adjustment = NULL, convert_units = 0.001)
plot(nest.model1, nc = 12)
# The above plots the detection distance probability with a hazard-rate function 


rd.full.df <- subset (raw.digs.df, select = -Zone); names(rd.full.df)[names(rd.full.df)=="Site"] <- "Sample.Label"; rd.full.df[2] <- NULL 
rd.full.df$Sample.Label <- paste(rd.full.df$Session, rd.full.df$Sample.Label, sep="_"); rd.full.df[5] <- NULL

split_full_rd <- split(rd.full.df, f = rd.full.df$Sample.Label)  
split_full_rd

rd_dist.df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(rd_dist.df) <- c('Sample.ID','Estimate','se','cv','lcl','ucl','df')  

# Run Hr for all, but Hn for Quail
# Join data together
# check that the models are assigned appropriately
for(i in 1:length(split_full_rd)) {       
  row <- split_full_rd[[i]]
  print(i)
  print(row)
  if (nrow(row) > 1) {
    tryCatch({ 
      rd.hr <- ds(data=row, key="hr", adjustment=NULL, convert_units=cf)
      print(rd.hr$dht$individuals$D)
      new_row <- rd.hr$dht$individuals$D
      new_row[new_row == 'Total'] <- row[1,4]# CHECK THESE ROWS WORK (ABOVE AND BELOW)
      rd_dist.df <- rbind(rd_dist.df, new_row)
    }, error=function(e) {
      print('Error')
    })
  }
}


plot(rd.hr, breaks = cutpoints)
summary(rd.hr)

test <- ds(data=rd.full.df, key="hr", adjustment=NULL, convert_units=cf)
print(test$dht$individuals$D)

summary(test)
cutpoints <- c(0,1,2,3,4,5,6,7,8,9,10)
plot(test, breaks = cutpoints, main = "Hazard rate model, quenda foraging line transects")


rd_dist.df 
write.csv(rd_dist.df, "C:/Users/tom_m/Desktop/Uni/PhD_Stats/Dig_veg_data_sheets/rd_dist.csv")
#copy paste this onto "dig.veg.csv" for stats

# Then get the density estimates using hn for Quail and add separately
q_rd_dist.df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(q_rd_dist.df) <- c('Sample.ID','Estimate','se','cv','lcl','ucl','df')  

q.rd.hn <- ds(data=split_full_rd$Spring21_Quail, key="hn", adjustment=NULL, convert_units=cf)
new_row <- q.rd.hn$dht$individuals$D; new_row[new_row == 'Total'] <- "Spring21_Quail"
q_rd_dist.df <- rbind(q_rd_dist.df, new_row)

q.rd.hn <- ds(data=split_full_rd$Winter21_Quail, key="hn", adjustment=NULL, convert_units=cf)
new_row <- q.rd.hn$dht$individuals$D; new_row[new_row == 'Total'] <- "Winter21_Quail"
q_rd_dist.df <- rbind(q_rd_dist.df, new_row)

q.rd.hn <- ds(data=split_full_rd$Summer22_Quail, key="hn", adjustment=NULL, convert_units=cf)
new_row <- q.rd.hn$dht$individuals$D; new_row[new_row == 'Total'] <- "Summer22_Quail"
q_rd_dist.df <- rbind(q_rd_dist.df, new_row)

q.rd.hn <- ds(data=split_full_rd$Summer22_Quail, key="hn", adjustment=NULL, convert_units=cf)
new_row <- q.rd.hn$dht$individuals$D; new_row[new_row == 'Total'] <- "Autumn22_Quail"
q_rd_dist.df <- rbind(q_rd_dist.df, new_row)

write.csv(q_rd_dist.df, "C:/Users/tom_m/Desktop/Uni/PhD_Stats/Dig_veg_data_sheets/quail_dist.csv")

#copy paste this onto "dig.veg.csv" for stats




# ~ ----

# STATISTICS ----

setwd("C:/Users/tom_m/Desktop/Uni/PhD_Stats/Dig_veg_data_sheets")
setwd("F:/PhD_Murdoch/PhD_Stats/Dig_veg_data_sheets")

library(dplyr);library(lme4);library(lmerTest);library(emmeans);library(tweedie);library(statmod)
library(piecewiseSEM);library(MuMIn);library(ggeffects);library(ggplot2); library(cowplot); library(DHARMa)
library(effects)

## df digging estimates & habitat----
dig.veg.df = read.csv("dig.veg.df.csv", header = T,sep = ",", stringsAsFactors = F)
dig.veg.df$Zone <- recode(dig.veg.df$Zone, Uninfested = 'Non-infested')
dig.veg.df$Zone <- factor(x=dig.veg.df$Zone, levels=c("Non-infested","Infested")) #to set order
dig.veg.df$Session <- factor(x=dig.veg.df$Session, levels=c("Winter21","Spring21","Summer22","Autumn22")) #to set order

## Prep & Check LMM data ----
# Zone * Session
str(dig.veg.df)

dig.veg.df$sample.id<-as.factor(dig.veg.df$sample.id)
dig.veg.df$Site<-as.factor(dig.veg.df$Site)
dig.veg.df$Transect<-as.factor(dig.veg.df$Transect)
dig.veg.df$Zone<-as.factor(dig.veg.df$Zone)
dig.veg.df$Vegetation<-as.factor(dig.veg.df$Vegetation)
dig.veg.df$Fire<-as.factor(dig.veg.df$Fire)
dig.veg.df$Soil<-as.factor(dig.veg.df$Soil)
dig.veg.df$Session<-as.factor(dig.veg.df$Session)

dig.veg.df$GT<-as.numeric(dig.veg.df$GT)
dig.veg.df$HGT<-as.numeric(dig.veg.df$HGT)
dig.veg.df$Gtree.tunnels<-as.numeric(dig.veg.df$Gtree.tunnels)


dig.veg.df$Shrub<-as.numeric(dig.veg.df$Shrub)
dig.veg.df$Subshrub<-as.numeric(dig.veg.df$Subshrub)
dig.veg.df$Herb<-as.numeric(dig.veg.df$Herb)
dig.veg.df$Graminoid<-as.numeric(dig.veg.df$Graminoid)

dig.veg.df$Litter<-as.numeric(dig.veg.df$Litter)
dig.veg.df$Bare<-as.numeric(dig.veg.df$Bare)




dotchart(x = dig.veg.df$GT, xlab = "GT", ylab = "order of the data from file")
dotchart(x = dig.veg.df$HGT, xlab = "HGT", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Used.HGT, xlab = "Used.HGT", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Shrub, xlab = "Shrub", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Subshrub, xlab = "Subshrub", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Herb, xlab = "Herb", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Graminoid, xlab = "Graminoid", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Litter, xlab = "Litter", ylab = "order of the data from file")
dotchart(x = dig.veg.df$Bare, xlab = "Bare", ylab = "order of the data from file")

dotchart(x = dig.veg.df$Density.Estimate, xlab = "Diggings", ylab = "order of the data from file")


# Datasets ----

dig.df <- dig.veg.df[-c(7,8,16,47,48,56,87,88,96,127,128,136),]
veg.df <- dig.veg.df[dig.veg.df$Session == "Winter21", ]



## Pearson's Litter:bare ----


correlation = cor.test(veg.df$Litter, veg.df$Bare, method = 'pearson')

correlation
# t = -1.331, df = 38, p-value = 0.1911 n/s
# 95% confidence interval: -0.4903213  0.1075316
# Correlation (r): -0.2110473

plot(veg.df$Litter ~ veg.df$Bare)


# ~ ----

# LMM Habitat Z*V*S(S) ----

## LMM Lifeform Total Cover
LMShrub <- lmer(Shrub ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
#options(na.action = "na.fail"); d_LMShrub <- dredge(LMShrub); options(na.action = "na.omit"); d_LMShrub
summary(LMShrub)
anova(LMShrub)
rsquared(LMShrub)
plot(resid(LMShrub))
hist(resid(LMShrub))
qqnorm(resid(LMShrub)); qqline(resid(LMShrub))
plot(simulateResiduals(LMShrub)) # good

LMSubshrub <- lmer(Subshrub ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
#options(na.action = "na.fail"); d_LMSubshrub <- dredge(LMSubshrub); options(na.action = "na.omit"); d_LMSubshrub
summary(LMSubshrub)
anova(LMSubshrub)
rsquared(LMSubshrub)
plot(resid(LMSubshrub))
hist(resid(LMSubshrub))
qqnorm(resid(LMSubshrub)); qqline(resid(LMSubshrub))
plot(simulateResiduals(LMSubshrub)) # good

LMHerb <- lmer(Herb ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
#options(na.action = "na.fail"); d_LMHerb <- dredge(LMHerb); options(na.action = "na.omit"); d_LMHerb
summary(LMHerb)
anova(LMHerb)
rsquared(LMHerb)
plot(resid(LMHerb))
hist(resid(LMHerb))
qqnorm(resid(LMHerb)); qqline(resid(LMHerb))
plot(simulateResiduals(LMHerb)) # quantile deviations detected


LMHerbb <- lmer(log(Herb) ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
summary(LMHerbb)
anova(LMHerbb)
rsquared(LMHerbb)
plot(resid(LMHerbb))
hist(resid(LMHerbb))
qqnorm(resid(LMHerbb)); qqline(resid(LMHerbb))
plot(simulateResiduals(LMHerbb)) # no problems anymore


LMGraminoid <- lmer(Graminoid ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
#options(na.action = "na.fail"); d_LMGraminoid <- dredge(LMGraminoid); options(na.action = "na.omit"); d_LMGraminoid
summary(LMGraminoid)
anova(LMGraminoid)
rsquared(LMGraminoid)
plot(resid(LMGraminoid))
hist(resid(LMGraminoid))
qqnorm(resid(LMGraminoid)); qqline(resid(LMGraminoid))
plot(simulateResiduals(LMGraminoid)) # good


## LMM Grasstrees
GT.df <- veg.df[-c(7,8,16),]

LMGT <- lmer(GT ~ Zone + Vegetation + Soil + (1|Site/Transect), data = GT.df)
#options(na.action = "na.fail"); d_LMGT <- dredge(LMGT); options(na.action = "na.omit"); d_LMGT
summary(LMGT)
anova(LMGT)
rsquared(LMGT)
plot(resid(LMGT))
hist(resid(LMGT))
qqnorm(resid(LMGT)); qqline(resid(LMGT))
plot(simulateResiduals(LMGT)) # good

LMHGT <- lmer(HGT ~ Zone + Vegetation + Soil + (1|Site/Transect), data = GT.df)
#options(na.action = "na.fail"); d_LMHGT <- dredge(LMHGT); options(na.action = "na.omit"); d_LMHGT
summary(LMHGT)
anova(LMHGT)
rsquared(LMHGT)
plot(resid(LMHGT))
hist(resid(LMHGT))
qqnorm(resid(LMHGT)); qqline(resid(LMHGT))
plot(simulateResiduals(LMHGT))  # good


# LMM Environmental Coverages
LMLit <- lmer(Litter ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
#options(na.action = "na.fail"); d_LMLit <- dredge(LMLit); options(na.action = "na.omit"); d_LMLit
summary(LMLit)
anova(LMLit)
rsquared(LMLit)
plot(resid(LMLit))
hist(resid(LMLit))
qqnorm(resid(LMLit)); qqline(resid(LMLit))
plot(simulateResiduals(LMLit)) # good

LMBare <- lmer(Bare ~ Zone + Vegetation + Soil + (1|Site/Transect), data = veg.df)
#options(na.action = "na.fail"); d_LMBare <- dredge(LMBare); options(na.action = "na.omit"); d_LMBare
summary(LMBare)
anova(LMBare)
rsquared(LMBare)
plot(resid(LMBare))
hist(resid(LMBare))
qqnorm(resid(LMBare)); qqline(resid(LMBare))
plot(simulateResiduals(LMBare)) # good

# Estimates
ggeffect(LMShrub, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMSubshrub, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMHerb, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMGraminoid, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMGT, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMHGT, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMLit, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMBare, "Zone", ci.lvl = 0.95, back.transform = TRUE)

ggeffect(LMShrub, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMSubshrub, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMHerb, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMGraminoid, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMGT, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMHGT, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMLit, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMBare, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)

ggeffect(LMShrub, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMSubshrub, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMHerb, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMGraminoid, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMGT, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMHGT, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMLit, "Soil", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LMBare, "Soil", ci.lvl = 0.95, back.transform = TRUE)


## GGplot ----
cbPalette <- c("chartreuse4", "tan1")
cbPalette2 <- c("#042E01", "#3E3600")


# GGplot Lifeform Total Cover
TC1 <- ggplot() + 
  geom_boxplot(data=veg.df, aes(x=Zone, y=Shrub, col = Zone, fill = Zone), 
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=veg.df, aes(x=Zone, y=Shrub, color = Zone, fill = Zone), 
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  #annotate(geom = "text", x = 1.5, y = 20, label = "**", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette, labels = c("Non-infested", "Infested")) +
  scale_colour_manual(values = cbPalette2, labels = c("Non-infested", "Infested"))  +
  xlab("") + 
  ylab("Tall Shrub Cover (cm)") +
  ggtitle("(b)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
        legend.key.size = unit(2,"line"), legend.title = element_blank(), 
        legend.text = element_text(size = 12))
TC1

TC2 <- ggplot() + 
  geom_boxplot(data=veg.df, aes(x=Zone, y=Subshrub, col = Zone, fill = Zone), 
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=veg.df, aes(x=Zone, y=Subshrub, color = Zone, fill = Zone), 
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 2500, label = "***", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Low Shrub Cover (cm)") +
  ggtitle("(a)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


TC3 <- ggplot() + 
  geom_boxplot(data=veg.df, aes(x=Zone, y=Herb, col = Zone, fill = Zone), 
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=veg.df, aes(x=Zone, y=Herb, color = Zone, fill = Zone), 
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75),) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Herb Cover (cm)") +
  ggtitle("(c)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


TC4 <- ggplot() + 
  geom_boxplot(data=veg.df, aes(x=Zone, y=Graminoid, col = Zone, fill = Zone), 
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=veg.df, aes(x=Zone, y=Graminoid, color = Zone, fill = Zone), 
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  #annotate(geom = "text", x = 1.5, y = 170, label = "n/s", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Graminoid Cover (cm)") +
  ggtitle("(d)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


# GGplot Grass Tree Densities
GT.df <- veg.df[-c(7,8,16),]
GT1 <- ggplot() + 
  geom_boxplot(data=GT.df, aes(x=Zone, y=HGT, col = Zone, fill = Zone),
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=GT.df, aes(x=Zone, y=HGT, color = Zone, fill = Zone),
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 57, label = "**", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Habitable Xp/600"~m^2) +
  ggtitle("(b)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


GT2 <- ggplot() + 
  geom_boxplot(data=GT.df, aes(x=Zone, y=GT, col = Zone, fill = Zone),
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=GT.df, aes(x=Zone, y=GT, color = Zone, fill = Zone),
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 135, label = "", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Xp/600"~m^2) +
  ggtitle("(a)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

# GGplot Environmental Cover
EC1 <- ggplot() + 
  geom_boxplot(data=veg.df, aes(x=Zone, y=Litter, col = Zone, fill = Zone),
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=veg.df, aes(x=Zone, y=Litter, color = Zone, fill = Zone), 
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 85, label = "", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Leaf Litter Cover (%)") +
  ggtitle("(f)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

EC2 <- ggplot() + 
  geom_boxplot(data=veg.df, aes(x=Zone, y=Bare, col = Zone, fill = Zone),
               outlier.size = 0, outlier.shape = NA)  +
  geom_point(data=veg.df, aes(x=Zone, y=Bare, color = Zone, fill = Zone), 
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 42, label = "*", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette2)  +
  xlab("") + 
  ylab("Bare Ground Cover (%)") +
  ggtitle("(e)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


## Plotting to pdf  ----

# the output will be in the datasheet folder (source)
library(ggpubr)

ggarrange(GT2,GT1,TC1,TC2,TC3,TC4,EC1,EC2, # list of plots
          labels = NA, # labels
          common.legend = T, # COMMON LEGEND
          legend = "bottom", # legend position
          align = "hv", # Align them both, horizontal and vertical
          nrow = 2,ncol = 4) 
# save pdf as 7 x 11.50 inches

 

# ~ ----

# Pearson's correlation tunnels HGT ----

GT.dfb <- veg.df[-c(7,8,16,40),]

correlation = cor.test(GT.dfb$Gtree.tunnels, GT.dfb$HGT, method = 'pearson')

correlation
# t = 6.3186, df = 34, p-value = 3.332e-07 ***
# 95% confidence interval: 0.5357028 0.8566131
# Correlation (r): 0.7348983 

library(ggplot2)

cbPalette <- c("chartreuse4", "tan1")

TunHGT <- ggplot(GT.dfb) +
  geom_smooth(data=GT.dfb, aes(x=HGT, y=Gtree.tunnels), 
              method = lm, linetype = "dashed", color = "black", size = 1) +
  geom_point(data=GT.dfb, aes(x=HGT, y=Gtree.tunnels, fill = Zone, shape = Zone), size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  annotate(geom = "text", x = 8, y = 20, label = "r = 0.73", colour="black", lwd = 5) +
  scale_fill_manual(values = cbPalette, labels = c("Non-infested", "Infested")) +
  ylab("Tunnels/600"~m^2) + 
  xlab("Habitable Xp/600"~m^2) +
  ggtitle("(d)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(legend.position="none")
TunHGT


# ~ ----


# GGplot grass tree tunnels ----
# Chi-squared stats were done in excel

# Make df
Zone <- c("Infested", "Infested","Healthy", "Healthy")
Skirt_Use <- c("Unused", "Tunnel","Unused", "Tunnel")
Skirt_Mean <- c(11.2941176470588,2.82352941176471,21,7.35294117647059)
Skirt_Upper <- c(12.5428766622481,3.67606089369739,22.4218604630395,8.3873359550205)
Skirt_Lower <- c(10.0453586318696,1.97099792983202,19.5781395369605,6.31854639792067)
Skirt_Total <- c(192,48,357,125)

skirt.df <- data.frame(Zone,Skirt_Use,Skirt_Mean,Skirt_Upper,Skirt_Lower)
str(skirt.df)

skirt.df$Zone <- factor(x=skirt.df$Zone, levels=c("Healthy","Infested")) #to set order
skirt.df$Skirt_Use <- factor(x=skirt.df$Skirt_Use, levels=c("Unused","Tunnel")) #to set order

cbPalette <- c("chartreuse4", "tan1")

# Means1
Skirt_Mean <- ggplot(data = skirt.df, aes(x=Skirt_Use, y=Skirt_Mean, fill=Zone)) +
  geom_bar(stat = "identity", colour = "black",position=position_dodge(.9), show.legend = FALSE)+
  geom_errorbar(aes(ymin=Skirt_Lower,ymax=Skirt_Upper), width=.2, position=position_dodge(.9))+
  guides(fill=guide_legend(title="Zone")) +
  scale_fill_manual(values = c("chartreuse4", "tan1"), labels = c("Non-infested", "Infested")) +
  xlab("Xp Skirt Condition") + 
  ylab("Mean Xp/600"~m^2) +
  ggtitle("(e)") +
  theme(axis.title.y = element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.2),
        strip.text = element_text(face = "italic")) +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_text(size = 11))

# ~ ----

# GGplot grass tree heights ----
# Chi-squared stats were done in excel

# get df
gt.height.df = read.csv("Gtree.height.ggplot.csv", header = T,sep = ",", stringsAsFactors = F)

gt.height.df2 <- gt.height.df ## to not overwrite the original DF
gt.height.df2[-1] <- lapply(gt.height.df2[-1], chartr, old = "_", new = "-")

gt.height.df2$Height_Class <- factor(x=gt.height.df2$Height_Class, 
                                    levels=c("0","1-9","10-20","21-30","31-40","41-50",
                                             "51-75","76-100","101-150","151+")) #to set order

str(gt.height.df2)
gt.height.df2$Zone<-as.factor(gt.height.df2$Zone)
gt.height.df2$Mean_Count<-as.numeric(gt.height.df2$Mean_Count)
gt.height.df2$Upper<-as.numeric(gt.height.df2$Upper)
gt.height.df2$Lower<-as.numeric(gt.height.df2$Lower)
gt.height.df2$Total.Count<-as.numeric(gt.height.df2$Total.Count)


cbPalette <- c("chartreuse4", "tan1")

# Means

Gtree_Height_Mean <- ggplot(data = gt.height.df2, aes(x=Height_Class, y=Mean_Count, fill=Zone)) +
  geom_bar(stat = "identity", colour = "black",position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Lower,ymax=Upper), width=.2, position=position_dodge(.9))+
  guides(fill=guide_legend(title="Zone")) +
  annotate(geom = "text", x = 1, y = 19.5, label = "*", colour="black", size = 14) +
  annotate(geom = "text", x = 7, y = 9.5, label = "*", colour="black", size = 14) +
  annotate(geom = "text", x = 8, y = 6, label = "**", colour="black", size = 14) +
  scale_fill_manual(values = c("chartreuse4", "tan1"), labels = c("Non-infested", "Infested")) +
  xlab("Height class (cm)") + 
  ylab("Mean Xp/600"~m^2) +
  ggtitle("(c)") +
  theme(axis.title.y = element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 0.2),
        strip.text = element_text(face = "italic")) +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_text(size = 11), legend.position = c(0.8,0.8), 
        legend.key.size = unit(1.5,'cm'), legend.title = element_text(face = 'bold', size = 18),
        legend.text = element_text(size = 16))

Gtree_Height_Mean



# ~ ----

# Grass tree panel ----

# new versions of GT plots
cbPalette <- c("chartreuse4", "tan1")
cbPalette2 <- c("#042E01", "#3E3600")


# GGplot Grass Tree Densities
GT.df <- veg.df[-c(7,8,16),]
GT1 <- ggplot() + 
  geom_boxplot(data=GT.df, aes(x=Zone, y=HGT, col = Zone, fill = Zone), 
               outlier.size = 0, outlier.shape = NA, show.legend = FALSE)  +
  geom_point(data=GT.df, aes(x=Zone, y=HGT, color = Zone, fill = Zone), show.legend = FALSE,
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 58, label = "**", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette, labels = c("Non-infested", "Infested")) +
  scale_colour_manual(values = cbPalette2, labels = c("Non-infested", "Infested"))  +
  xlab("") + 
  ylab("Habitable Xp/600"~m^2) +
  ggtitle("(b)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


GT2 <- ggplot() + 
  geom_boxplot(data=GT.df, aes(x=Zone, y=GT, col = Zone, fill = Zone),
               outlier.size = 0, outlier.shape = NA, show.legend = FALSE)  +
  geom_point(data=GT.df, aes(x=Zone, y=GT, color = Zone, fill = Zone), show.legend = FALSE,
             shape = 21, size = 3, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) +
  annotate(geom = "text", x = 1.5, y = 135, label = "", colour="black", size = 14) +
  scale_fill_manual(values = cbPalette, labels = c("Non-infested", "Infested")) +
  scale_colour_manual(values = cbPalette2, labels = c("Non-infested", "Infested"))  +
  xlab("") + 
  ylab("Xp/600"~m^2) +
  ggtitle("(a)") +
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())




# the output will be in the datasheet folder (source)

library(cowplot)

cow_gtree_ab <- plot_grid(GT2,GT1,
                            align = "h",
                            ncol = 2,
                            nrow = 1,
                            #labels = c("(a)","(b)"),
                            label_size = 14,
                            label_fontfamily = "sans")

cow_gtree_c <- plot_grid(NULL,Gtree_Height_Mean,
                         nrow = 1,
                         rel_widths = c(0.01,1),
                         #labels = c("(c)"),
                         label_size = 14,
                         label_fontfamily = "sans")

cow_gtree_de <- plot_grid(TunHGT,Skirt_Mean,
                          align = "h",
                          ncol = 2,
                          nrow = 1,
                          greedy = TRUE,
                          #labels = c("(d)","(e)"),
                          label_size = 14,
                          label_fontfamily = "sans")

cow_gtree <- plot_grid(cow_gtree_ab,cow_gtree_c,cow_gtree_de,
                       rel_widths = c(2,2),
                       align = 'v',
                       nrow = 3)
cow_gtree # 15 x 9

# Vegetation panel 2 ----
# Using this instead of the one under the stats:

# the output will be in the datasheet folder (source)
library(ggpubr)

ggarrange(TC2,TC1,TC3,TC4,EC2,EC1, # list of plots
          labels = NA, # labels
          common.legend = T, # COMMON LEGEND
          legend = "bottom", # legend position
          align = "hv", # Align them both, horizontal and vertical
          nrow = 2,ncol = 3) 
# save pdf as 7 x 8.5 inches landscape




# ~ ----

# LMM Digging Estimates----

# Checking if log is important to include - may absorb the extreme outliers
LM1 <- lmer(Density.Estimate ~ Zone + Session + Vegetation + Soil + (1|Site/Transect), data = dig.df)
options(na.action = "na.fail"); d_LM1 <- dredge(LM1); options(na.action = "na.omit"); d_LM1
plot(resid(LM1))
hist(resid(LM1))
qqnorm(resid(LM1)); qqline(resid(LM1))

LM2 <- lmer(log(Density.Estimate + 1) ~ Zone + Session + Vegetation + Soil +(1|Site/Transect), data = dig.df)
options(na.action = "na.fail"); d_LM2 <- dredge(LM2); options(na.action = "na.omit"); d_LM2
plot(resid(LM2))
hist(resid(LM2))
qqnorm(resid(LM2)); qqline(resid(LM2))
# QQPLOT is slightly better

# Dredging this to see which factors are important
LM3 <- lmer(log(Density.Estimate + 1) ~ Zone + Session + Vegetation + Soil + 
              GT + HGT + Litter + Bare + Shrub + Subshrub + Herb + Graminoid + (1|Site/Transect), 
            data = dig.df)
options(na.action = "na.fail"); d_LM3 <- dredge(LM3); options(na.action = "na.omit"); d_LM3
plot(resid(LM3))
hist(resid(LM3))
qqnorm(resid(LM3)); qqline(resid(LM3))

# Zone + Vegetation is the best model

# Checking to see which of these is the best (including Session and Soil as randoms since not main terms now)
LM4 <- lmer(log(Density.Estimate + 1) ~ Zone + Vegetation + (1|Site/Transect), data = dig.df)
LM5 <- lmer(log(Density.Estimate + 1) ~ Zone + Vegetation + (1|Site/Transect) + (1|Session), data = dig.df)
LM6 <- lmer(log(Density.Estimate + 1) ~ Zone + Vegetation + (1|Site/Transect) + (1|Soil), data = dig.df)
LM7 <- lmer(log(Density.Estimate + 1) ~ Zone + Vegetation + (1|Site/Transect) + (1|Session) + (1|Soil), data = dig.df)
anova(LM4, LM5, LM6, LM7)
#LM5 is the lowest AIC

LM5 <- lmer(log(Density.Estimate + 1) ~ Zone + Vegetation + (1|Site/Transect) + (1|Session), data = dig.df)
summary(LM5)
anova(LM5)
rsquared(LM5)
emmeans(LM5, list(pairwise ~ Zone*Vegetation), adjust = "tukey")
plot(resid(LM5))
hist(resid(LM5)) # bell shaped curve
qqnorm(resid(LM5)); qqline(resid(LM5))
plot(simulateResiduals(LM5)) # within-group deviations from uniformity significant (Levene's)
                             # I think this is because of the time-based nature

testDispersion(LM5) # dispersion not significant


# Estimates
ggeffect(LM5, "Zone", ci.lvl = 0.95, back.transform = TRUE)
ggeffect(LM5, "Vegetation", ci.lvl = 0.95, back.transform = TRUE)


## GGplot ----

# Graphing all of these results:
cbPalette <- c("chartreuse4", "tan1")
cbPalette2 <- c("#042E01", "#3E3600")

# Also graphing without the outlier:
ggplot() +
  geom_boxplot(data=dig.df, aes(x=Session, y=Density.Estimate, col = Zone, fill = Zone), 
               outlier.size = 0, outlier.color = NULL, outlier.alpha = 0, outlier.stroke = 0.5)  +
  geom_point(data=dig.df, aes(x=Session, y=Density.Estimate, color = Zone, fill = Zone), 
             shape = 21, size= 2.8, alpha = 0.8, position=position_jitterdodge(dodge.width=0.75)) + 
  scale_fill_manual(values = cbPalette, labels = c("Non-infested", "Infested")) +
  scale_colour_manual(values = cbPalette2, labels = c("Non-infested", "Infested"))  +
  labs(y='Foraging pit density estimates'~(m^2), x='Foraging Period') + 
  labs(fill='Zone*', col = 'Zone*') +  
  theme_classic(base_size = 14, base_line_size = 0.5) +
  theme(axis.text.x = element_text(color="black"),
        axis.title.x = element_text(vjust=-0.5),
        axis.text.y = element_text(color="black"), 
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.18,0.9)) +
  scale_x_discrete(labels = c("Apr May Jun", "Jul Aug Sep", "Oct Nov Dec", "Jan Feb Mar"))
# save as pdf 5.5 x 5.5
