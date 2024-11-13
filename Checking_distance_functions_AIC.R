
setwd("C:/Users/tom_m/Desktop/Uni/PhD_Stats/Dig_veg_data_sheets")

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






