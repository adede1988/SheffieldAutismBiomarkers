



library(lmtest)
library(lsr)
library(tidyverse)
library(car)
library(circular)
library(kableExtra)
library(pracma)
library(glmtoolbox)
source('G:\\My Drive\\GitHub/SheffieldAutismBiomarkers/helperFuncsFinal.R')

## set path to github repo on local machine here: 
path = "G:\\My Drive\\GitHub\\SheffieldAutismBiomarkers\\"
# path = "C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\"
figPath = "G:\\My Drive\\Milne\\pubFigPanels\\"
test <- read.csv(paste(path, "autismBiomarkersAllData3.csv", sep = ''))

test$nbChanOrig[test$nbChanOrig==999] = 124 #due to a data import error, less than 10 participants had their nbChanOrig values missing

dat <- filter(test, group %in% c('CON', 'AD', 'ASD'), eyes == 'open')

#AD subjects above and below 250 months of age
sum(dat$group=='AD' & dat$age > 250) 
sum(dat$group=='AD' & dat$age <= 250)


dat <- filter(dat, age<=250)

dat <- filter(dat, nbChanOrig>20)
dat <- filter(dat, nbChanFinal/nbChanOrig >= .5)

#9 subjects had missing IQ data! these subjects will have dataset/group appropriate means substituted in below
sum(is.na(dat$IQ))

#data description table: 
makeDemoTable(dat)


#age group detection originally done with this code: 
# ageGroups = quantile(filter(dat, group %in% c('AD', 'ASD'))$age , c(0,.33333,.666666,1.0))
# ageGroups[1]= 0
# ageGroups[4]=999
#hard code the age groups 
ageGroups = c(0,97,127,999)


#### randomization code written on 29/11/2022 
#see helperFuncs.R for code specifications 


#### end randomization block ####
# dat = read.csv(paste(path, "testSet.csv", sep= ''))
# 
# dat <- dat %>% select(-X)
varNames = names(dat)[20:1179]


test <- read.csv(paste(path, "autismBiomarkersAllData3.csv", sep = ''))

test$nbChanOrig[test$nbChanOrig==999] = 124 #due to a data import error, less than 10 participants had their nbChanOrig values missing

dat <- filter(test, group %in% c('CON', 'AD', 'ASD'), eyes == 'open')



#AD subjects above and below 250 months of age
sum(dat$group=='AD' & dat$age > 250) 
sum(dat$group=='AD' & dat$age <= 250)


dat <- filter(dat, age<=250)

dat <- filter(dat, nbChanOrig>20)
dat <- filter(dat, nbChanFinal/nbChanOrig >= .5)

# ageGroups = c(0, 999999) #use this to test what happens if all included in one!
# varNames = varNames[randsample(length(varNames),30)]
#run stats on the full dataset
results <- getResults(dat, varNames, ageGroups)
dat <- results[[2]]
results <- results[[1]]

#do analysis a second time but this time do it without the severity split: 
# dat2 = dat
# dat2$group[dat2$group == 'ASD'] = 'AD'
# results2 = getResults(dat2, varNames, ageGroups)
# dat2 <- results2[[2]]
# results2 <- results2[[1]]

repeats = 150
savePath = "R:\\MSS\\Johnson_Lab\\dtf8829\\GitHub\\SheffieldAutismBiomarkers\\bootstrapAnalysis\\"

#what's the stability of the large effect sizes for diagnosis? 
if(!file.exists(paste(savePath, "Diag_stab.csv", sep = ''))){
test <- results %>% filter(Diag > .035)
stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
stabilityMat <- as.data.frame(stabilityMat)
names(stabilityMat) <- test$dependentVariable
lookup <- paste(test$dependentVariable, test$ageGroup)


  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$Diag>.035 & tmpRes2$Diag>.035
  
  }
  stabilityMat_Diag <- stabilityMat
  write.csv(stabilityMat_Diag,
            paste(savePath, "Diag_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "Diag_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_Diag <- stabilityMat

#diag X age
if(!file.exists(paste(savePath, "Diag_Age_stab.csv", sep = ''))){
  test <- results %>% filter(Diag_Age > .035)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$Diag_Age>.035 & tmpRes2$Diag_Age>.035
    
  }
  stabilityMat_Diag_Age <- stabilityMat
  write.csv(stabilityMat_Diag_Age,
            paste(savePath, "Diag_Age_stab.csv", sep = ''))
}else{
  stabilityMat <- read.csv(paste(savePath, "Diag_Age_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_Diag_Age <- stabilityMat

#diag X sex
if(!file.exists(paste(savePath, "Diag_Sex_stab.csv", sep = ''))){
  test <- results %>% filter(Diag_Sex > .035)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$Diag_Sex>.035 & tmpRes2$Diag_Sex>.035
  }
  write.csv(stabilityMat, paste(savePath, "Diag_Sex_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "Diag_Sex_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_Diag_Sex <- stabilityMat

#diag X sex X age
if(!file.exists(paste(savePath, "Diag_Age_Sex_stab.csv", sep = ''))){
  test <- results %>% filter(Diag_Age_Sex > .035)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$Diag_Age_Sex>.035 & tmpRes2$Diag_Age_Sex>.035
  }
  write.csv(stabilityMat, paste(savePath, "Diag_Age_Sex_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "Diag_Age_Sex_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_Diag_Age_Sex <- stabilityMat

#age
if(!file.exists(paste(savePath, "Age_stab.csv", sep = ''))){
  test <- results %>% filter(age > .035 & age<.06)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$age>.035 & tmpRes2$age>.035
  }
  write.csv(stabilityMat, paste(savePath, "Age_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "Age_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_Age <- stabilityMat

#sex
if(!file.exists(paste(savePath, "Sex_stab.csv", sep = ''))){
  test <- results %>% filter(sex > .035 & sex<.06)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$sex>.035 & tmpRes2$sex>.035
  }
  write.csv(stabilityMat, paste(savePath, "Sex_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "Sex_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_Sex <- stabilityMat

#age high
if(!file.exists(paste(savePath, "AgeHigh_stab.csv", sep = ''))){
  test <- results %>% filter(age>.06)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$age>.035 & tmpRes2$age>.035
  }
  write.csv(stabilityMat, paste(savePath, "AgeHigh_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "AgeHigh_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_AgeHigh <- stabilityMat

#sex high
if(!file.exists(paste(savePath, "SexHigh_stab.csv", sep = ''))){
  test <- results %>% filter(sex>.06)
  stabilityMat <- matrix(0, nrow = repeats, ncol = length(test$dependentVariable))
  stabilityMat <- as.data.frame(stabilityMat)
  names(stabilityMat) <- test$dependentVariable
  lookup <- paste(test$dependentVariable, test$ageGroup)
  
  for(ii in 1:repeats){
    print(ii)
    tmpDat = randomSplit(dat, ageGroups)
    tmpRes1 <- getResults(tmpDat[[1]], test$dependentVariable, ageGroups)
    tmpRes1 <- tmpRes1[[1]]
    screen = paste(tmpRes1$dependentVariable, tmpRes1$ageGroup)
    tmpRes1 <- tmpRes1[screen %in% lookup, ]
    tmpRes2 <- getResults(tmpDat[[2]], test$dependentVariable, ageGroups)
    tmpRes2 <- tmpRes2[[1]]
    tmpRes2 <- tmpRes2[screen %in% lookup, ]
    
    stabilityMat[ii,] <- tmpRes1$sex>.035 & tmpRes2$sex>.035
  }
  write.csv(stabilityMat, paste(savePath, "SexHigh_stab.csv", sep = ''))
} else {
  stabilityMat <- read.csv(paste(savePath, "SexHigh_stab.csv", sep = ''))
  stabilityMat <- stabilityMat %>% select(-X)
}
stabilityMat_SexHigh <- stabilityMat



results$Diag_stab <- NA
results$Diag_stab[which(results$Diag>.035)] <- colMeans(stabilityMat_Diag)

results$Diag_Age_stab <- NA
results$Diag_Age_stab[which(results$Diag_Age>.035)] <- colMeans(stabilityMat_Diag_Age)

results$Diag_Sex_stab <- NA
results$Diag_Sex_stab[which(results$Diag_Sex>.035)] <- colMeans(stabilityMat_Diag_Sex)

results$Diag_Age_Sex_stab <- NA
results$Diag_Age_Sex_stab[which(results$Diag_Age_Sex>.035)] <- colMeans(stabilityMat_Diag_Age_Sex)

results$Age_stab <- NA
results$Age_stab[which(results$age>.035 & results$age<.06)] <- colMeans(stabilityMat_Age)
results$AgeHigh_stab <- NA
results$AgeHigh_stab[which(results$age>.06 )] <- colMeans(stabilityMat_AgeHigh)
results$Age_stab[which(results$age>.06 )] <- colMeans(stabilityMat_AgeHigh)

results$Sex_stab <- NA
results$Sex_stab[which(results$sex>.035 & results$sex<.06)] <- colMeans(stabilityMat_Sex)
results$SexHigh_stab <- NA
results$SexHigh_stab[which(results$sex>.06)] <- colMeans(stabilityMat_SexHigh)
results$Sex_stab[which(results$sex>.06)] <- colMeans(stabilityMat_SexHigh)


#### what proportion of EEG variables are predicted by different IVs? ####
#Table 3
results %>% 
  group_by(ageGroup, type) %>% 
  summarize(age = sum(age>.035) / length(age),
            sex = sum(sex>.035) / length(sex), 
            IQ = sum(IQ>.035)/length(IQ), 
            Diag = sum(Diag>.035)/length(Diag), 
            Diag_Age = sum(Diag_Age>.035)/length(Diag_Age), 
            Diag_Sex = sum(Diag_Sex>.035)/length(Diag_Sex), 
            Diag_Age_Sex = sum(Diag_Age_Sex>.035)/length(Diag_Age_Sex))


#### what proportion would have been better fit with a quadratic?  ####
sum(results$nonLin<.05) / length(results$nonLin)
  
 results %>% group_by(type, ageGroup) %>% summarize(nonLinCount = sum(nonLin<.05)/length(nonLin))
    
    
test <- results %>% filter(Diag > .05)

results %>% ggplot(aes(x = Diag, y = DiagBinom, color = ageGroup)) + 
  geom_jitter() + 
  geom_abline(intercept = .05, slope = 0) + 
  geom_vline(xintercept = .05)

results %>% ggplot(aes(x = Diag, y = Diag_resid, color = as.factor(ageGroup))) + 
  geom_jitter() + 
  geom_vline(xintercept = .05)

results %>% ggplot(aes(x = Diag, y = Diag_resid, color = as.factor(ageGroup))) + 
  geom_jitter() + 
  geom_vline(xintercept = .05)


#### general graphics params ####

predNames = c('age', 'sex', 'IQ', 'Diag', 'Diag_Age', 'Diag_Sex', 'Diag_Age_Sex')
stabNames = c('Diag_stab', 'Diag_Age_stab', 'Diag_Sex_stab', 'Diag_Age_Sex_stab',
              'Age_stab', 'AgeHigh_stab', 'Sex_stab', 'SexHigh_stab')
# stabNames = c('Diag_stab', 'Diag_Age_stab', 'Diag_Sex_stab', 'Diag_Age_Sex_stab',
#               'Age_stab',  'Sex_stab', 'SexHigh_stab')
ageLabs = c(paste(as.character(min(filter(dat, age>ageGroups[1] & age<ageGroups[2])$age)), '-',
                  as.character(max(filter(dat, age>ageGroups[1] & age<ageGroups[2])$age)), 
                  ' months', sep = ''), 
            paste(as.character(min(filter(dat, age>ageGroups[2] & age<ageGroups[3])$age)), '-',
                  as.character(max(filter(dat, age>ageGroups[2] & age<ageGroups[3])$age)), 
                  ' months', sep = ''),
            paste(as.character(min(filter(dat, age>ageGroups[3] & age<ageGroups[4])$age)), '-',
                  as.character(max(filter(dat, age>ageGroups[3] & age<ageGroups[4])$age)), 
                  ' months', sep = ''))
names(ageLabs) <- c('1', '2', '3')
myColors = c('#F22F29','#CD4D5D', '#FB9564', '#4ED6B2','#2C64C0', '#1B47A5',    '#7B5C52')
ageColors = c('#551956', '#91a250', '#cc99cc')
names(myColors) <- levels(factor(results$type))
colScale <- scale_color_manual(name = "type",values = myColors)
fillScale <- scale_fill_manual(name = "type",values = myColors)
theme_set(theme_gray(base_size = 30))
results$ageGroup = as.factor(results$ageGroup)
#### Evaluate outliers ####
sum(results$out==0) / length(results$out)
#what was the maximum number of outliers? 
max(results$out)


results %>% group_by(type) %>% summarize(max(out), mean(out), min(out))

outlierCounts <- dat %>%
  group_by(group) %>%
  summarize(across(everything(), ~ sum(is.na(.))))

rowMeans(outlierCounts[,2:1179])


ggplot(results, aes(x=out, color = type, fill = type)) +
  colScale +
  fillScale + 
  geom_histogram(breaks = seq(0,50,1)) +
  coord_cartesian(xlim=c(0,50),ylim = c(0,100))



#### plotting the effect size values ####
thresh = .00


for(vari in 1:7){
  for(agei in 1:3){
    outPlot<-ggplot(filter(results, results[predNames[vari]]>thresh & results$ageGroup == agei),
                    aes_string(x=predNames[vari], color = 'type', fill = 'type')) +
      colScale +
      fillScale + 
      xlab('effect size (\U1D702\U00B2\U209A)') +
      geom_histogram(breaks = seq(0,1,.035)) +
      coord_cartesian(xlim=c(0,.350001),ylim = c(0,100)) +
      theme_classic(base_size = 30)+
      theme(axis.ticks.length=unit(-0.5, "cm"),
            axis.line = element_line(colour = 'black', size = 3),
            axis.ticks = element_line(colour = 'black', size = 3)) + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, .65), breaks = seq(0,.6,.035)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0,100,20)) + 
      geom_vline(xintercept = .035, linetype = 'dashed', size = 5)
    
    png(paste(figPath, 'effect_ALL_', predNames[vari], '_', as.character(agei), '.png',sep=''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }  
}

for(vari in 1:8){
  outPlot<-ggplot(filter(results, !is.na(results[stabNames[vari]])),
                  aes_string(x=stabNames[vari], color = 'type', fill = 'type')) +
    colScale +
    fillScale + 
    xlab('agreement') +
    geom_histogram(breaks = seq(0,1,.08)) +
    theme_classic(base_size = 30)+
    theme(axis.ticks.length=unit(-0.5, "cm"),
          axis.line = element_line(colour = 'black', size = 3),
          axis.ticks = element_line(colour = 'black', size = 3)) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0,1,.08)) +
    scale_y_continuous(expand = c(0, 0))+ 
    geom_vline(xintercept = .64, linetype = 'dashed', size = 5)
  
  png(paste(figPath, 'stability_New_', stabNames[vari], '_', '.png',sep=''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()     
}


effectSummary <- results %>% group_by(ageGroup, type) %>% 
  summarize(Age = round(sum(age>.035) / length(age),2),
            Sex = round(sum(sex>.035) / length(sex),2),
            IQ = round(sum(IQ>.035) / length(IQ), 2),
            Diag = round(sum(Diag>.035) / length(Diag),2), 
            Diag_Age = round(sum(Diag_Age>.035) / length(Diag_Age),2),
            Diag_Sex = round(sum(Diag_Sex>.035) / length(Diag_Sex),2),
            Diag_Age_Sex = round(sum(Diag_Age_Sex>.035) / length(Diag_Age_Sex),2)) 
write.csv(effectSummary, "effectSummary.csv")

results %>% 
  summarize(Diag = sum(Diag>.035) , 
            Diag_Age = sum(Diag_Age>.035) ,
            Diag_Sex = sum(Diag_Sex>.035),
            Diag_Age_Sex = sum(Diag_Age_Sex>.035 ) )

#### plotting relationship between effect size and p-value ####

outPlot <- ggplot(results,
                  aes(x=Diag, y = Diag_p, color = ageGroup)) + 
  geom_jitter(alpha = .5, size = 10) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, .081), breaks = seq(.005,1,.015)) +
  scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
  geom_hline(yintercept = .05, linewidth = 3) + 
  geom_vline(xintercept = .035, linewidth = 3) +
  theme_classic(base_size = 30)+
  ylim(c(0, .25)) +
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  xlab('effect size (\U1D702\U00B2\U209A)') + 
  ylab('p-value') 
  
png(paste(figPath, 'EffectSize_pValRelationship', '.png',sep=''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

#### plotting relationship between different effect sizes ####

diagTypes = c("Diag", "Diag_Age", "Diag_Sex", "Diag_Age_Sex")

for(ii in 1:4){
  for(jj in 1:4){
    outPlot <- ggplot(results,
           aes_string(x=diagTypes[ii], y = diagTypes[jj], color = "ageGroup")) + 
      scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
      geom_jitter(alpha = .5, size = 5) + 
      geom_hline(yintercept = .035) + 
      geom_vline(xintercept = .035) +
      theme_classic(base_size = 30)+
      xlim(c(0, .08)) +
      ylim(c(0, .08)) +
      theme(axis.ticks.length=unit(-0.5, "cm"),
            axis.line = element_line(colour = 'black', size = 3),
            axis.ticks = element_line(colour = 'black', size = 3), 
            legend.position = 'none') + 
      xlab(paste( diagTypes[ii], '(\U1D702\U00B2\U209A)')) + 
      ylab(paste( diagTypes[jj], '(\U1D702\U00B2\U209A)')) 
    
    png(paste(figPath, diagTypes[ii],'___', diagTypes[jj], '.png',sep=''),         # File name
        width=700, height=700)
    print(outPlot)
    dev.off()
    
  }
  
  
}


ggplot(results,
       aes(x=Diag, y = Diag_Age, color = as.factor(ageGroup))) + 
  scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
  geom_jitter(alpha = .5, size = 5) + 
  geom_hline(yintercept = .035) + 
  geom_vline(xintercept = .035) +
  theme_classic(base_size = 30)+
  xlim(c(0, .07)) +
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  xlab('Diag (\U1D702\U00B2\U209A)') + 
  ylab('Diag X Age (\U1D702\U00B2\U209A)') 

ggplot(results,
       aes(x=Diag, y = Diag_Sex, color = as.factor(ageGroup))) + 
  scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
  geom_jitter(alpha = .5, size = 5) + 
  geom_hline(yintercept = .035) + 
  geom_vline(xintercept = .035) +
  theme_classic(base_size = 30)+
  xlim(c(0, .07)) +
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  xlab('Diag (\U1D702\U00B2\U209A)') + 
  ylab('Diag X Sex (\U1D702\U00B2\U209A)') 

ggplot(results,
       aes(x=Diag, y = Diag_Age_Sex, color = as.factor(ageGroup))) + 
  scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
  geom_jitter(alpha = .5, size = 5) + 
  geom_hline(yintercept = .035) + 
  geom_vline(xintercept = .035) +
  theme_classic(base_size = 30)+
  xlim(c(0, .07)) +
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  xlab('Diag (\U1D702\U00B2\U209A)') + 
  ylab('Diag X Age X Sex (\U1D702\U00B2\U209A)') 


#### relationship between stability and effect size ####

sizeVal = 8

outPlot <- ggplot(results, aes(x = age, y = AgeHigh_stab)) + 
  geom_point(size = sizeVal/2) + 
  geom_point(aes(x=age, y = Age_stab),size = sizeVal/2) + 
  geom_point(aes(x=sex, y = Sex_stab),size = sizeVal/2) + 
  geom_point(aes(x=sex, y = SexHigh_stab),size = sizeVal/2) + 
  geom_point(aes(x=Diag, y = Diag_stab), color = 'red', size = sizeVal) +
  geom_point(aes(x=Diag_Age, y = Diag_Age_stab), color = 'red', size = sizeVal) +
  geom_point(aes(x=Diag_Sex, y = Diag_Sex_stab), color = 'red', size = sizeVal) +
  geom_point(aes(x=Diag_Age_Sex, y = Diag_Age_Sex_stab), color = 'red', size = sizeVal) +
  # xlim(c(.034, .1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01), breaks = seq(0,1,.1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(.034, .101), breaks = seq(.035,.1,.015)) +
  geom_hline(yintercept = .64, linewidth = 2, linetype = 'dashed') + 
  theme_classic(base_size = 30)+
  theme(axis.ticks.length=unit(-0.75, "cm"),
        axis.line = element_line(colour = 'black', size = sizeVal),
        axis.ticks = element_line(colour = 'black', size = sizeVal), 
        legend.position = 'none') + 
  xlab('effect size (\U1D702\U00B2\U209A)') + 
  ylab('stability') 

png(paste(figPath, 'EffectSize_stability', '.png',sep=''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

#### relationship between logistic regression vs multiple linear reg. ####

outPlot <- results %>% ggplot(aes(x = Diag, y = logAug - logBase, color = ageGroup)) + 
  geom_jitter(alpha = .5, size = 10) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, .081), breaks = seq(.005,1,.015)) +
  scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
  theme_classic(base_size = 30)+
  # ylim(c(0, .03)) +
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  # geom_smooth(data = results, method = 'lm', aes(x = Diag, y = logAug - logBase)) + 
  xlab('effect size (\U1D702\U00B2\U209A)') + 
  ylab('logistic regression (AUC improvement)') 

png(paste(figPath, 'logisticVlinear', '.png',sep=''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

cor.test(results$Diag, results$logAug-results$logBase)




outPlot <-ggplot(results, aes(x = Diag, y = Diag_residR2, color = ageGroup)) + 
  geom_jitter(alpha = .5, size = 10) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-.001, .081), breaks = seq(.005,1,.015)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-.001, .081), breaks = seq(.005,1,.015)) +
  scale_color_manual(values = setNames(ageColors, levels(results$ageGroup)))+
  theme_classic(base_size = 30)+
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  # geom_smooth(data = results, method = 'lm', aes(x = Diag, y = logAug - logBase)) + 
  xlab('effect size (\U1D702\U00B2\U209A)') + 
  ylab('effect size (residual model; R^2)') 
png(paste(figPath, 'residualVlinear', '.png',sep=''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

cor.test(results$Diag, results$Diag_residR2)

#### plotting specific group differences ####
dat$ageGroup = 1
dat$ageGroup[dat$age>ageGroups[2]] = 2
dat$ageGroup[dat$age>ageGroups[3]] = 3
dat$ageGroup = as.factor(dat$ageGroup)

test <- results %>% filter(Diag>.035, Diag_stab > .64)
for(tt in 1:length(test$Diag)){
 
  targName = test$dependentVariable[tt]

  limVals = quantile(dat[,targName], c(.1,.9), na.rm = T)
  limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
              limVals[2] + (limVals[2]-limVals[1])*.5)

  outPlot <-  dat %>% filter(ageGroup == test$ageGroup[tt]) %>%
    ggplot(aes_string(x = 'ageGroup', y = targName, fill = 'group')) +
    geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.85)) + #width = .6, coef = 0
    geom_jitter(show.legend = F, shape = 21, size = 4, 
                position = position_jitterdodge(dodge.width = .85, jitter.width = .2)) +
    # geom_violin(draw_quantiles = c(.25, .5, .75))+
    scale_fill_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                      labels=c(  'AD', 'ASD', 'CON')) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20),
          legend.position = 'none') +
    ylim(limVals) +
    ylab(targName) +
    ggtitle(paste(targName, '; age group: ', as.character(test$ageGroup[tt]), sep = ''))
  # rect(1, 5, 3, 7, col="white")
  png(paste(figPath, 'bestPredicted_',targName,'_', test$ageGroup[tt], '.png', sep = ''),         # File name
      width=512, height=768)
  print(outPlot)
  dev.off()

}

##diag X sex interactions 

test <- results %>% filter(Diag_Sex>.035, Diag_Sex_stab > .64)
for(tt in 1:length(test$Diag)){

  targName = test$dependentVariable[tt]
  
  limVals = quantile(dat[,targName], c(.1,.9), na.rm = T)
  limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
              limVals[2] + (limVals[2]-limVals[1])*.5)
  
outPlot <-  dat %>% filter(ageGroup == test$ageGroup[tt]) %>%
  ggplot(aes_string(x = 'sex', y = targName, fill = 'group')) +
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.85)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4, 
              position = position_jitterdodge(dodge.width = .85, jitter.width = .2)) +
  # geom_violin(draw_quantiles = c(.25, .5, .75))+
  scale_fill_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                    labels=c(  'AD', 'ASD', 'CON')) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        legend.position = 'none') +
  ylim(limVals) +
  ylab(targName) +
  ggtitle(paste(targName, '; age group: ', as.character(test$ageGroup[tt]), sep = ''))
# geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
# geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
# rect(1, 5, 3, 7, col="white")
png(paste(figPath, 'bestPredictedSEX_',targName,'_', test$ageGroup[tt], '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

}

#diag X age
test <- results %>% filter(Diag_Age>.035, Diag_Age_stab > .64)
for(tt in 1:length(test$Diag)){
  
  targName = test$dependentVariable[tt]
  
  limVals = quantile(dat[,targName], c(.1,.9), na.rm = T)
  limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
              limVals[2] + (limVals[2]-limVals[1])*.5)
  
  outPlot <-  dat %>% filter(ageGroup == test$ageGroup[tt]) %>%
    ggplot(aes_string(x = 'age', y = targName, fill = 'group',
                      color = 'group', group = 'group')) +
    geom_jitter(show.legend = F, shape = 21, size = 4, 
                position = position_jitterdodge(dodge.width = .85, jitter.width = .2)) +
    # geom_violin(draw_quantiles = c(.25, .5, .75))+
    scale_fill_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                      labels=c(  'AD', 'ASD', 'CON')) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20),
          legend.position = 'none') +
    geom_smooth(method = 'lm', se = F, linewidth = 3) +
    scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                      labels=c(  'AD', 'ASD', 'CON')) +
    ylim(limVals) +
    ylab(targName) +
    ggtitle(paste(targName, '; age group: ', as.character(test$ageGroup[tt]), sep = ''))
  # geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
  # geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
  # rect(1, 5, 3, 7, col="white")
  png(paste(figPath, 'bestPredictedAGE_',targName,'_', test$ageGroup[tt], '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
}

#three way interaction plot
test <- results %>% filter(Diag_Age_Sex>.035, Diag_Age_Sex_stab > .64)
for(tt in 1:length(test$Diag)){
  
  targName = test$dependentVariable[tt]
  
  limVals = quantile(dat[,targName], c(.1,.9), na.rm = T)
  limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
              limVals[2] + (limVals[2]-limVals[1])*.5)
  
  outPlot <-  dat %>% filter(ageGroup == test$ageGroup[tt]) %>%
    ggplot(aes_string(x = 'age', y = targName, fill = 'group',
                      color = 'group', group = 'group')) +
    geom_jitter(show.legend = F, shape = 21, size = 4, 
                position = position_jitterdodge(dodge.width = .85, jitter.width = .2)) +
    # geom_violin(draw_quantiles = c(.25, .5, .75))+
    scale_fill_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                      labels=c(  'AD', 'ASD', 'CON')) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20),
          legend.position = 'none') +
    geom_smooth(method = 'lm', se = F, linewidth = 3) +
    scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                       labels=c(  'AD', 'ASD', 'CON')) +
    ylim(limVals) +
    ylab(targName) +
    ggtitle(paste(targName, '; age group: ', as.character(test$ageGroup[tt]), sep = '')) + 
    facet_grid(~sex)
  # geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
  # geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
  # rect(1, 5, 3, 7, col="white")
  png(paste(figPath, 'bestPredictedINT_',targName,'_', test$ageGroup[tt], '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
}

#### bootstrap with smaller sample sizes ##############################################
library(randnet)
allNMI <-data.frame('sampProp' = rep(NA, 100000), 
                    'sampN' = rep(NA, 100000), 
                    'diagNMI' = rep(NA, 100000), 
                    'ageNMI' = rep(NA, 100000), 
                    'sexNMI' = rep(NA, 100000), 
                    'diagProp' = rep(NA, 100000),
                    'ageProp' = rep(NA, 100000),
                    'sexProp' = rep(NA, 100000),
                    'diagAgree' = rep(NA, 100000),
                    'ageAgree' = rep(NA, 100000),
                    'sexAgree' = rep(NA, 100000))
for(ss in 1:7){
  if(file.exists(paste(savePath, "NMIvals", ss, ".csv", sep = ''))){
    NMIvals = read.csv(paste(savePath, "NMIvals", ss, ".csv", sep = ''))
    NMIvals <- NMIvals[,-1]
  } else {
    
    
    options(warn=-1)
    sampSize = seq(.1, .7, .1)
    NMIvals = data.frame('sampProp' = rep(NA, 10000), 
                         'sampN' = rep(NA, 10000), 
                         'diagNMI' = rep(NA, 10000), 
                         'ageNMI' = rep(NA, 10000), 
                         'sexNMI' = rep(NA, 10000), 
                         'diagProp' = rep(NA, 10000),
                         'ageProp' = rep(NA, 10000),
                         'sexProp' = rep(NA, 10000),
                         'diagAgree' = rep(NA, 10000),
                         'ageAgree' = rep(NA, 10000),
                         'sexAgree' = rep(NA, 10000))
    ni = 1
    
    boots = 100
    stabilityMat <- matrix(0, nrow = 10, ncol = length(results$dependentVariable))
    stabilityMat <- as.data.frame(stabilityMat)
    names(stabilityMat) <- results$dependentVariable
    stabilityMatAge = stabilityMat
    stabilityMatSex = stabilityMat
    for(bb in 1:boots){
      curDat = randomSplit(dat, ageGroups, sampSize[ss])[[2]]
      curResults <- getResults(curDat, varNames, ageGroups)
      curResults <- curResults[[1]]
      diag = curResults$Diag
      stabilityMat[bb,] = diag>.035
      stabilityMatAge[bb,] = curResults$age >.035
      stabilityMatSex[bb,] = curResults$sex >.035
      
    }
    #get the NMI values and other summary stats
    for(ii in 1:boots){
      for(jj in 1:boots){
        if(ii>jj){ #there's no directionality
          NMIvals$diagNMI[ni] = NMI(as.numeric(as.vector(stabilityMat[ii,])), as.numeric(as.vector(stabilityMat[jj,])))  
          NMIvals$ageNMI[ni] = NMI(as.numeric(as.vector(stabilityMatAge[ii,])), as.numeric(as.vector(stabilityMatAge[jj,])))  
          NMIvals$sexNMI[ni] = NMI(as.numeric(as.vector(stabilityMatSex[ii,])), as.numeric(as.vector(stabilityMatSex[jj,])))
          NMIvals$diagAgree[ni] = sum(as.numeric(as.vector(stabilityMat[ii,]))==as.numeric(as.vector(stabilityMat[jj,])))/ 2184  
          NMIvals$ageAgree[ni] = sum(as.numeric(as.vector(stabilityMatAge[ii,]))==as.numeric(as.vector(stabilityMatAge[jj,])))/ 2184  
          NMIvals$sexAgree[ni] = sum(as.numeric(as.vector(stabilityMatSex[ii,]))==as.numeric(as.vector(stabilityMatSex[jj,])))/ 2184  
          NMIvals$sampProp[ni] = sampSize[ss]
          NMIvals$sampN[ni] = round(mean(curResults$n))
          NMIvals$diagProp[ni] = mean(rowSums(stabilityMat[c(ii,jj),]) / 2184)
          NMIvals$ageProp[ni] = mean(rowSums(stabilityMatAge[c(ii,jj),]) / 2184)
          NMIvals$sexProp[ni] = mean(rowSums(stabilityMatSex[c(ii,jj),]) / 2184)
          
          ni = ni+1
          
          
          
        }
      }
    }
    NMIvals <- NMIvals %>% filter(!is.na(diagProp))
    write.csv(NMIvals,
              paste(savePath, "NMIvals", ss, ".csv", sep = ''))
  }
  allNMI[((ss-1)*4950+1):(ss*4950),] <- NMIvals
}
allNMI <- allNMI %>% filter(!is.na(diagProp))

#### does anything fundamentally change with balanced data? ####
datBal <- getBalanced(dat)

resultsBal <- getResults(datBal, varNames, ageGroups)
resultsBal <- resultsBal[[1]]




outPlot <-allNMI %>% ggplot(aes(x = sampN, y = diagNMI)) + 
  geom_jitter(alpha = .1, width= 13,color = '#E1C271') + 
  geom_jitter(alpha = .1, width= 13,aes(x = sampN, y = ageNMI), color = '#2f0404') +
  geom_jitter(alpha = .1, width= 13,aes(x = sampN, y = sexNMI), color = '#BF40BF') +
  theme_classic(base_size = 30)+
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  # geom_smooth(data = results, method = 'lm', aes(x = Diag, y = logAug - logBase)) + 
  xlab('sample size') + 
  ylab('replicability (normalized mutual information)') 
png(paste(figPath, 'NMIplot', '.png',sep=''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

test = matrix(seq(1,34650,4950), nrow = 100, ncol = 7, byrow = T)
test2 = matrix(1:100, nrow = 100, ncol = 7)
test3 = test + test2 +300
allNMIsamp = allNMI[c(test3),]
outPlot <- allNMIsamp %>% ggplot(aes(x = sampN, y = diagProp)) + 
  geom_jitter(alpha = .5, width= 8, size = 5, color = '#E1C271') + 
  geom_jitter(alpha = .5, size = 5, aes(x = sampN, y = ageProp), color = '#2f0404', width= 8) +
  geom_jitter(alpha = .5, size = 5, aes(x = sampN, y = sexProp), color = '#BF40BF', width= 8)+ 
  geom_point(data=resultsBal, aes(x = mean(n),
          y = sum(Diag>.035)/length(Diag)), color = 'black', size = 10,
          shape = 21, stroke = 3, fill = 'transparent') +
  geom_point(data=resultsBal, aes(x = mean(n), 
          y = sum(Diag>.035)/length(Diag)), color = '#9dce95', size = 10, alpha = .02) + 
  theme_classic(base_size = 30)+
  theme(axis.ticks.length=unit(-0.5, "cm"),
        axis.line = element_line(colour = 'black', size = 3),
        axis.ticks = element_line(colour = 'black', size = 3), 
        legend.position = 'none') + 
  # geom_smooth(data = results, method = 'lm', aes(x = Diag, y = logAug - logBase)) + 
  xlab('sample size') + 
  ylab('proportion \U1D702\U00B2\U209A > .035') 
png(paste(figPath, 'propSigPlot', '.png',sep=''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()


#### get out a table of the final stable values ####


test <- results %>% select(dependentVariable, ageGroup, Diag, Diag_Age, Diag_Sex, Diag_Age_Sex, Diag_stab, Diag_Age_stab, Diag_Sex_stab, Diag_Age_Sex_stab)

df_long <- test %>% select(-c(Diag_stab, Diag_Age_stab, Diag_Sex_stab, Diag_Age_Sex_stab)) %>%
  pivot_longer(cols = c(Diag, Diag_Age, Diag_Sex, Diag_Age_Sex),
               names_to = "IV",
               values_to = "etaSquared") 

# Pivot longer for stability values
df_stability <- test %>% select(-c(Diag, Diag_Age, Diag_Sex, Diag_Age_Sex)) %>%
  pivot_longer(cols = c(Diag_stab, Diag_Age_stab, Diag_Sex_stab, Diag_Age_Sex_stab),
               names_to = "IV_stability",
               values_to = "stability")

# Extract the base IV names from the stability variables for merging
df_stability <- df_stability %>%
  mutate(IV = sub("_stab", "", IV_stability)) %>%
  select(-IV_stability)

# Merge the etaSquared and stability dataframes
df_final <- df_long %>%
  left_join(df_stability %>% select(dependentVariable, IV, stability), 
            by = c("dependentVariable", "IV")) %>%
  filter(stability>.64, etaSquared >.035)

# Select the desired columns
df_final <- df_final %>%
  select(dependentVariable, etaSquared, IV, stability)

write.csv(df_final, paste(figPath, 'stableVars', '.csv',sep=''))

#### summarize the variability of measures ####

dat$ageGroup = 1
dat$ageGroup[dat$age>=ageGroups[2] & dat$age<ageGroups[3]] = 2
dat$ageGroup[dat$age>=ageGroups[3] & dat$age<ageGroups[4]] = 3
dat$ageGroup <- as.factor(dat$ageGroup)

nonPhaseVarNames = varNames[unlist(lapply(varNames, 
                            function(x) !grepl('PACmi', x ) &
                              !grepl('phase', x )))]
#z-score all of the data within age group
standardDat <- dat[c('group', 'ageGroup', 'sex', nonPhaseVarNames)]
for(agei in 1:3) {
  standardDat[,4:length(names(standardDat))] %>%
    filter(dat$ageGroup == agei) %>%
    apply(2, { function(x)
                    (x - mean(x, na.rm = T))^2 }
          ) ->
    standardDat[dat$ageGroup==agei,4:length(names(standardDat))]
}

AD = standardDat %>% filter(group == "AD" | group == "ASD")
CON = standardDat %>% filter(group == "CON")
ASD = standardDat %>% filter(group == 'ASD')
M = standardDat %>% filter(sex == 'M')
Fe = standardDat %>% filter(sex == 'F')

varRes = matrix(0, nrow = length(results$dependentVariable), ncol = 4)
varRes <- as.data.frame(varRes)
varRes[,1] <- results$dependentVariable
names(varRes) <- c('variable', 'var', 'ageGroup', 'sexVar')



for(agei in 1:3) {
  ageAdjust <- ((agei-1)*((length(names(AD))-3))+1)
  L = length(names(AD)) - 3
  
  AD[,4:length(names(AD))] %>% 
    filter(AD$ageGroup == agei) %>%
    apply(2, function(x) var(x, na.rm = T) ) -> 
    ADvar
  M[,4:length(names(AD))] %>% 
    filter(M$ageGroup == agei) %>%
    apply(2, function(x) var(x, na.rm = T) ) -> 
    Mvar
  Fe[,4:length(names(AD))] %>% 
    filter(Fe$ageGroup == agei) %>%
    apply(2, function(x) var(x, na.rm = T) ) -> 
    Fevar
  
  CON[,4:length(names(CON))] %>% 
    filter(CON$ageGroup == agei) %>%
    apply(2, function(x) var(x, na.rm = T) ) -> 
    conVar
  
  varRes$sexVar[ageAdjust:(ageAdjust+L-1)] <- (Mvar - Fevar) / (Mvar + Fevar)
  varRes$var[ageAdjust:(ageAdjust+L-1)] <- (ADvar - conVar) / (ADvar + conVar)
  varRes$ageGroup[ageAdjust:(ageAdjust+L-1)] <- agei
  # hist((ADvar - conVar) / (ADvar + conVar))
}

ggplot(varRes, aes(y = sexVar, color = as.factor(ageGroup))) + geom_boxplot()
ggplot(varRes, aes(y = var, color = as.factor(ageGroup))) + geom_boxplot()

varRes %>% group_by(ageGroup) %>% 
  summarize(M = mean(var), S = sd(var), d = mean(var) - sd(var))
varRes %>% group_by(ageGroup) %>% 
  summarize(M = mean(sexVar), S = sd(sexVar), d = mean(sexVar) - sd(sexVar))


# 
# for(agei in 1:3) {
#   ageAdjust <- ((agei-1)*(3*(length(names(AD))-2))+1)
#   L = length(names(AD)) - 2
#   varRes$ageGroup[ageAdjust:(ageAdjust+L*3-1)] <- agei
#   
#   AD[,3:length(names(AD))] %>% 
#     filter(AD$ageGroup == agei) %>% 
#     apply(2, function(x) var(x, na.rm = T) ) -> 
#     varRes$var[ageAdjust:(ageAdjust+L-1)]
#   'AD' -> varRes$diagGroup[ageAdjust:(ageAdjust+L-1)]
#   
#   ASD[,3:length(names(AD))] %>% 
#     filter(ASD$ageGroup == agei) %>% 
#     apply(2, function(x) var(x, na.rm = T) ) -> 
#     varRes$var[(ageAdjust+L):(ageAdjust+L*2-1)]
#   'ASD' -> varRes$diagGroup[(ageAdjust+L):(ageAdjust+L*2-1)]
#   
#   CON[,3:length(names(CON))] %>% 
#     filter(CON$ageGroup == agei) %>% 
#     apply(2, function(x) var(x, na.rm = T) ) -> 
#     varRes$var[(ageAdjust+L*2):(ageAdjust+L*3-1)]
#   'CON' -> varRes$diagGroup[(ageAdjust+L*2):(ageAdjust+L*3-1)]
# }



#### scratch ####

#run stats on the training set 
dat = read.csv(paste(path, "trainSet.csv", sep= ''))
# dat = read.csv("C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
dat <- dat %>% select(-X)
resultsTmp <- getResults(dat, varNames, ageGroups)
datTrain <- resultsTmp[[2]]
resultsTRAIN <- resultsTmp[[1]]
#run stats on the test dataset 
dat = read.csv(paste(path, "testSet.csv", sep=''))
# dat = read.csv("C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\testSet.csv")
dat <- dat %>% select(-X)
resultsTmp <- getResults(dat, varNames, ageGroups)
datTest <- resultsTmp[[2]]
resultsTEST <- resultsTmp[[1]]





