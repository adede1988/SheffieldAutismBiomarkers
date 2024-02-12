





library(lsr)
library(tidyverse)
library(car)
library(circular)
library(kableExtra)
library(pracma)
source('R:/MSS/Johnson_Lab/dtf8829/GitHub/SheffieldAutismBiomarkers/helperFuncs.R')

## set path to github repo on local machine here: 
path = "R:\\MSS\\Johnson_Lab\\dtf8829\\GitHub\\SheffieldAutismBiomarkers\\"
# path = "C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\"
figPath = "G:\\My Drive\\Milne\\pubFigPanels\\"
test <- read.csv(paste(path, "autismBiomarkersAllData2.csv", sep = ''))

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


#### randomization code written on 29/11/2022 and run one time
#see helperFuncs.R for code specifications 


#### end randomization block ####
dat = read.csv(paste(path, "trainSet.csv", sep= ''))
# dat = read.csv("C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
dat <- dat %>% select(-X)
varNames = names(dat)[20:1179]

results = data.frame('dependentVariable' = rep('A', 4*length(varNames)),
                     'age' = rep(0, 4*length(varNames)),
                     'sex' = rep(0, 4*length(varNames)),
                     'IQ' = rep(0, 4*length(varNames)),
                     'Diag' = rep(0, 4*length(varNames)),
                     'age_b' = rep(0, 4*length(varNames)),
                     'sex_b' = rep(0, 4*length(varNames)),
                     'IQ_b' = rep(0, 4*length(varNames)),
                     'diag_b1' = rep(0, 4*length(varNames)),
                     'diag_b2' = rep(0, 4*length(varNames)),
                     'n' = rep(0, 4*length(varNames)),
                     'out' = rep(0, 4*length(varNames)),
                     'outlierFirst' = rep(0, 4*length(varNames)),
                     'ageGroup' = rep(0, 4*length(varNames)))



for(ii in 2:length(ageGroups)){
  #down select for age
  curDat = dat[dat$age>ageGroups[ii-1] & dat$age<=ageGroups[ii], ]
  #down select for group ID
  # curDat = curDat[curDat$group == 'CON' | curDat$group == 'AD' |curDat$group == 'ASD',]
  groupID = rep('CON',length(curDat$group))
  
  groupID[curDat$group=='AD'] = 'AD'
  groupID[curDat$group=='ASD'] = 'ASD'
  
  outliers = matrix(0,length(varNames), 3)

  for(tt in 1:length(varNames)){
    if(!grepl('PACmi', varNames[tt] )) { #skipping the non z-scored raw PAC values
    print(tt)
    temp = curDat
    #### outlier removal ####
    #after beginning analysis, it was discovered that some variables contained
    #extreme outliers such that one or two observations were far above the mean
    #and all other variables were below the mean.
    #in these cases, simply eliminating outliers based on a standard z-score
    #would have eliminated all variables since the mean ended up in a position
    #far from any data points. To get around this, data were z scored using only
    #the middle 80% of the data. Then, based on this z-scoring, outliers were removed
    tmpOut <- getModDat(dat, varNames, tt, temp, ageGroups, ii)
    modDat <- tmpOut[[1]]
    dat <- tmpOut[[2]]
   

    #optional plot for this variable
    # modelPlot(modDat, varNames[tt], ii)



    if(!grepl('phase', varNames[tt])){
    
    modDat$diag <- factor(modDat$diag, levels = c("CON", "ASD", "AD"))
    modDat$sex <- factor(modDat$sex, levels = c( "F","M"))  
    curLM = lm(dv ~ age + sex * diag  +IQ , data = modDat)
    aovTab = Anova(curLM, type = 3)
    curLM$coefficients
    ai = (ii-2)*length(varNames)
    
    results[tt+ai,1] = varNames[tt]
    ei = which(row.names(aovTab) == 'Residuals')
    ri = which(row.names(aovTab) == 'diag')
    results$Diag[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
    ri = which(row.names(aovTab) == 'sex:diag')
    results$sexDiag[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
    ri = which(row.names(aovTab) == 'sex')
    results$sex[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
    ri = which(row.names(aovTab) == 'IQ')
    results$IQ[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
    ri = which(row.names(aovTab) == 'age')
    results$age[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
    
    ri = which(names(curLM$coefficients) == 'age')
    results$age_b[tt+ai] = curLM$coefficients[[ri]]
    ri = which(names(curLM$coefficients) == 'sexM')
    results$sex_b[tt+ai] = curLM$coefficients[[ri]]
    ri = which(names(curLM$coefficients) == 'IQ')
    results$IQ_b[tt+ai] = curLM$coefficients[[ri]]
    ri = which(names(curLM$coefficients) == 'diagASD')
    results$diag_b1[tt+ai] = curLM$coefficients[[ri]]
    ri = which(names(curLM$coefficients) == 'diagAD')
    results$diag_b2[tt+ai] = curLM$coefficients[[ri]]
    
    
    results$n[tt+ai] = length(temp[,1])
    results$out[tt+ai] = length(temp[,1])-length(modDat$dv)
    results$ageGroup[tt+ai] = ii-1
    results$dvMean[tt+ai] = modDat$dv_mean[1]
    results$dvSD[tt+ai] = modDat$dv_sd[1]
    
    }else { #phase needs to be treated differently because it's circular
    circdv = circular(modDat$dv, unit = 'degrees')
    dv_mean =  mean(circdv)[[1]]
    dv_sd = sd(circdv)[[1]]
    modDat$dv2 = cbind(cos(modDat$dv*(pi/180)), sin(modDat$dv*(pi/180)))
    curLM = lm(dv2 ~ age + sex * diag  + IQ, data = modDat )
    aovTab = Anova(curLM, type = 3)
    outtests <- car:::print.Anova.mlm
    body(outtests)[[16]] <- quote(invisible(tests))
    body(outtests)[[15]] <- NULL
    tab <- lapply(c("Pillai"), 
                  function(i)  outtests(Anova(curLM, test.statistic=i)))
    aovTab <- do.call(rbind, tab)
    #working out effect size for diag
    ai = (ii-2)*length(varNames)
    ri =  which(row.names(aovTab) == 'diag')
    results[tt+ai,1] = varNames[tt]
    results$Diag[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
                          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
    ri =  which(row.names(aovTab) == 'sex')
    results$sex[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
                          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
    ri =  which(row.names(aovTab) == 'IQ')
    results$IQ[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
                          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
    ri =  which(row.names(aovTab) == 'age')
    results$age[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
                          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
    ri = which(row.names(aovTab) == 'sex:diag')
    results$sexDiag[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
      (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
    
    results$n[tt+ai] = length(temp[,1])
    results$out[tt+ai] = length(temp[,1])-length(modDat$dv)
    results$ageGroup[tt+ai] = ii-1
    results$dvMean[tt+ai] = dv_mean
    results$dvSD[tt+ai] = dv_sd
    }
    

  }
  
 
  }
}

resultsTRAIN = results
datTrain <- dat
#### now get the same set of results but from the test set! 
dat = read.csv(paste(path, "testSet.csv", sep=''))
# dat = read.csv("C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\testSet.csv")
dat <- dat %>% select(-X)
varNames = names(dat)[20:1179]
results = data.frame('dependentVariable' = rep('A', 4*length(varNames)),
                      'age' = rep(0, 4*length(varNames)),
                      'sex' = rep(0, 4*length(varNames)),
                      'IQ' = rep(0, 4*length(varNames)),
                      'Diag' = rep(0, 4*length(varNames)),
                      'age_b' = rep(0, 4*length(varNames)),
                      'sex_b' = rep(0, 4*length(varNames)),
                      'IQ_b' = rep(0, 4*length(varNames)),
                      'diag_b1' = rep(0, 4*length(varNames)),
                      'diag_b2' = rep(0, 4*length(varNames)),
                      'n' = rep(0, 4*length(varNames)),
                      'out' = rep(0, 4*length(varNames)),
                      'outlierFirst' = rep(0, 4*length(varNames)),
                      'ageGroup' = rep(0, 4*length(varNames)))


for(ii in 2:length(ageGroups)){
  #down select for age
  curDat = dat[dat$age>ageGroups[ii-1] & dat$age<=ageGroups[ii], ]
  #down select for group ID
  # curDat = curDat[curDat$group == 'CON' | curDat$group == 'AD' |curDat$group == 'ASD',]
  groupID = rep('CON',length(curDat$group))
  
  groupID[curDat$group=='AD'] = 'AD'
  groupID[curDat$group=='ASD'] = 'ASD'
  
  outliers = matrix(0,length(varNames), 3)
  
  for(tt in 1:length(varNames)){
    if(!grepl('PACmi', varNames[tt] )) { #skipping the non z-scored raw PAC values
      print(tt)
      temp = curDat
      #### outlier removal ####
      #after beginning analysis, it was discovered that some variables contained
      #extreme outliers such that one or two observations were far above the mean
      #and all other variables were below the mean.
      #in these cases, simply eliminating outliers based on a standard z-score
      #would have eliminated all variables since the mean ended up in a position
      #far from any data points. To get around this, data were z scored using only
      #the middle 80% of the data. Then, based on this z-scoring, outliers were removed
      tmpOut <- getModDat(dat, varNames, tt, temp, ageGroups, ii)
      modDat <- tmpOut[[1]]
      dat <- tmpOut[[2]]
      
      
      #optional plot for this variable
      # modelPlot(modDat, varNames[tt], ii)
      
      
      
      if(!grepl('phase', varNames[tt])){
        
        
        modDat$diag <- factor(modDat$diag, levels = c("CON", "ASD", "AD"))
        modDat$sex <- factor(modDat$sex, levels = c( "F","M"))  
        curLM = lm(dv ~ age + sex * diag  +IQ , data = modDat)
        aovTab = Anova(curLM, type = 3)
        curLM$coefficients
        ai = (ii-2)*length(varNames)
        
        results[tt+ai,1] = varNames[tt]
        ei = which(row.names(aovTab) == 'Residuals')
        ri = which(row.names(aovTab) == 'diag')
        results$Diag[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
        ri = which(row.names(aovTab) == 'sex:diag')
        results$sexDiag[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
        ri = which(row.names(aovTab) == 'sex')
        results$sex[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
        ri = which(row.names(aovTab) == 'IQ')
        results$IQ[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
        ri = which(row.names(aovTab) == 'age')
        results$age[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
        
        ri = which(names(curLM$coefficients) == 'age')
        results$age_b[tt+ai] = curLM$coefficients[[ri]]
        ri = which(names(curLM$coefficients) == 'sexM')
        results$sex_b[tt+ai] = curLM$coefficients[[ri]]
        ri = which(names(curLM$coefficients) == 'IQ')
        results$IQ_b[tt+ai] = curLM$coefficients[[ri]]
        ri = which(names(curLM$coefficients) == 'diagASD')
        results$diag_b1[tt+ai] = curLM$coefficients[[ri]]
        ri = which(names(curLM$coefficients) == 'diagAD')
        results$diag_b2[tt+ai] = curLM$coefficients[[ri]]
        
        
        results$n[tt+ai] = length(temp[,1])
        results$out[tt+ai] = length(temp[,1])-length(modDat$dv)
        results$ageGroup[tt+ai] = ii-1
        results$dvMean[tt+ai] = modDat$dv_mean[1]
        results$dvSD[tt+ai] = modDat$dv_sd[1]
        
      }else { #phase needs to be treated differently because it's circular
        circdv = circular(modDat$dv, unit = 'degrees')
        dv_mean =  mean(circdv)[[1]]
        dv_sd = sd(circdv)[[1]]
        modDat$dv2 = cbind(cos(modDat$dv*(pi/180)), sin(modDat$dv*(pi/180)))
        curLM = lm(dv2 ~ age + sex * diag  + IQ, data = modDat )
        aovTab = Anova(curLM, type = 3)
        outtests <- car:::print.Anova.mlm
        body(outtests)[[16]] <- quote(invisible(tests))
        body(outtests)[[15]] <- NULL
        tab <- lapply(c("Pillai"), 
                      function(i)  outtests(Anova(curLM, test.statistic=i)))
        aovTab <- do.call(rbind, tab)
        #working out effect size for diag
        ai = (ii-2)*length(varNames)
        ri =  which(row.names(aovTab) == 'diag')
        results[tt+ai,1] = varNames[tt]
        results$Diag[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
        ri =  which(row.names(aovTab) == 'sex')
        results$sex[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
        ri =  which(row.names(aovTab) == 'IQ')
        results$IQ[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
        ri =  which(row.names(aovTab) == 'age')
        results$age[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
        ri = which(row.names(aovTab) == 'sex:diag')
        results$sexDiag[tt+ai] = (aovTab$Df[1] * aovTab$`approx F`[ri]) / 
          (aovTab$Df[ri] * aovTab$`approx F`[ri] + aovTab$`den Df`[ri])
        
        results$n[tt+ai] = length(temp[,1])
        results$out[tt+ai] = length(temp[,1])-length(modDat$dv)
        results$ageGroup[tt+ai] = ii-1
        results$dvMean[tt+ai] = dv_mean
        results$dvSD[tt+ai] = dv_sd
      }
      
      
    }
    
    
  }
}

resultsTEST = results
#### result difference scores #### 
resultsDif = data.frame('dependentVariable' = resultsTEST$dependentVariable, 
                       'age' = 1 - (abs(resultsTEST$age - resultsTRAIN$age) / (resultsTEST$age + resultsTRAIN$age)),
                       'sex' = 1 - (abs(resultsTEST$sex - resultsTRAIN$sex) / (resultsTEST$sex + resultsTRAIN$sex)),
                       'IQ' = 1 - (abs(resultsTEST$IQ - resultsTRAIN$IQ) / (resultsTEST$IQ + resultsTRAIN$IQ)),
                       'Diag' = 1 - (abs(resultsTEST$Diag - resultsTRAIN$Diag) / (resultsTEST$Diag + resultsTRAIN$Diag)),
                       'sexDiag' = 1 - (abs(resultsTEST$sexDiag - resultsTRAIN$sexDiag) / 
                                          (resultsTEST$sexDiag + resultsTRAIN$sexDiag)),
                       'ageGroup' = resultsTEST$ageGroup )
resultsDif = filter(resultsDif, resultsTEST$dependentVariable != 'A')
resultsTEST = filter(resultsTEST, dependentVariable != 'A')
resultsTRAIN = filter(resultsTRAIN, dependentVariable != 'A')

resultsDif$type = 'A'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('pow', x, ignore.case = T))] = 'power'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('slope', x, ignore.case = T))] = 'slope'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('relalpha', x, ignore.case = T) | 
                        grepl('logalpha', x, ignore.case = T))] = 'alpha'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('sampEnt', x, ignore.case = T))] = 'sampEnt'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('ispc', x, ignore.case = T))] = 'ispc'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('PAC', x, ignore.case = T))] = 'PAC'
resultsDif$type[apply(as.matrix(resultsDif$dependentVariable), 
                      1, function(x) grepl('phase', x, ignore.case = T))] = 'phase'

resultsTRAIN$type = 'A'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('pow', x, ignore.case = T))] = 'power'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('slope', x, ignore.case = T))] = 'slope'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('relalpha', x, ignore.case = T) |
                        grepl('logalpha', x, ignore.case = T))] = 'alpha'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('sampEnt', x, ignore.case = T))] = 'sampEnt'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('ispc', x, ignore.case = T))] = 'ispc'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('PAC', x, ignore.case = T))] = 'PAC'
resultsTRAIN$type[apply(as.matrix(resultsTRAIN$dependentVariable),
                      1, function(x) grepl('phase', x, ignore.case = T))] = 'phase'





#### general graphics params ####

predNames = c('age', 'sex', 'IQ', 'Diag', 'sexDiag')
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
names(myColors) <- levels(factor(resultsTRAIN$type))
colScale <- scale_color_manual(name = "type",values = myColors)
fillScale <- scale_fill_manual(name = "type",values = myColors)
theme_set(theme_gray(base_size = 30))

#### evaluate outliers during analysis #### 

sum(resultsTRAIN$out==0) / length(resultsTRAIN$out)
ggplot(resultsTRAIN, aes(x=out, color = type, fill = type)) +
  colScale +
  fillScale + 
  geom_histogram(breaks = seq(0,22,1)) +
  coord_cartesian(xlim=c(0,22),ylim = c(0,100))


#### plotting the effect size values ####
thresh = .00


for(vari in 1:5){
  for(agei in 1:3){
    outPlot<-ggplot(filter(resultsTRAIN, resultsTRAIN[predNames[vari]]>thresh & resultsTRAIN$ageGroup == agei),
                    aes_string(x=predNames[vari], color = 'type', fill = 'type')) +
      colScale +
      fillScale + 
      xlab('effect size (\U1D702\U00B2\U209A)') +
      geom_histogram(breaks = seq(0,1,.05)) +
      coord_cartesian(xlim=c(0,.6500001),ylim = c(0,100)) +
      theme_classic(base_size = 30)+
      theme(axis.ticks.length=unit(-0.5, "cm"),
            axis.line = element_line(colour = 'black', size = 3),
            axis.ticks = element_line(colour = 'black', size = 3)) + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, .65), breaks = seq(0,.6,.1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0,100,20)) + 
      geom_vline(xintercept = .05, linetype = 'dashed', size = 5)
      
    png(paste(figPath, 'interaction_effect_', predNames[vari], '_', as.character(agei), '.png',sep=''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }  
}


#### what's the relationship between effect size and effect stability? ####
n = length(resultsDif$age)
allEffects = data.frame('agreement' = c(resultsDif$age, resultsDif$sex, resultsDif$IQ, resultsDif$sexDiag),
                        'effectSize' = c(resultsTRAIN$age, resultsTRAIN$sex, resultsTRAIN$IQ, resultsTRAIN$sexDiag),
                        'predictor' = c(rep('age', n), rep('sex', n), rep('IQ', n), rep('sex X diagnosis', n)))

  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'age')[,1],filter(allEffects, predictor == 'age')[,2])
  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'sex')[,1],filter(allEffects, predictor == 'sex')[,2])
  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'IQ')[,1],filter(allEffects, predictor == 'IQ')[,2])
  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'sex X diagnosis')[,1],filter(allEffects, predictor == 'sex X diagnosis')[,2])

outPlot <-ggplot(allEffects, aes(x=agreement, y=effectSize, color = predictor)) + 
                 geom_jitter(size = 5, alpha = .5)+
                 theme_classic(base_size = 30)+
                 theme(axis.ticks.length=unit(-0.5, "cm"),
                 axis.line = element_line(colour = 'black', size = 3),
                 axis.ticks = element_line(colour = 'black', size = 3),
                 legend.position = c(0.37, 0.8)) + 
                 guides(color = guide_legend(override.aes = list(size=10),
                                             nrow=4)) + 
                 ylab('effect size (\U1D702\U00B2\U209A)') + 
                 scale_x_continuous(expand = c(0, 0), limits = c(0, 1.01), breaks = seq(0,1,.25)) +
                 scale_y_continuous(expand = c(0, 0), limits = c(0, .65), breaks = seq(0,.65,.1))
png(paste(figPath, 'interaction_effect_agree_relationship.png', sep = ''),         # File name
    width=400, height=768)
print(outPlot)
dev.off()         

# ggplot(filter(allEffects, predictor == 'sex'), aes(x=effect, y=agreement)) + geom_jitter()


#### agreement between test and train data sets in model fit effect size ####
thresh = .05

for(vari in 1:5){
  # for(agei in 1:3){
    outPlot<-ggplot(filter(resultsDif, resultsTRAIN[predNames[vari]]>thresh),
                    aes_string(x=predNames[vari], color = 'type', fill = 'type')) +
      colScale +
      fillScale + 
      xlab('agreement') +
      geom_histogram(breaks = seq(0,1,.05)) +
      # coord_cartesian(xlim=c(0,1),ylim = c(0,150)) +
      theme_classic(base_size = 30)+
      theme(axis.ticks.length=unit(-0.5, "cm"),
            axis.line = element_line(colour = 'black', size = 3),
            axis.ticks = element_line(colour = 'black', size = 3)) + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0,1,.1)) +
      scale_y_continuous(expand = c(0, 0))+ 
      geom_vline(xintercept = .8, linetype = 'dashed', size = 5)
    
    png(paste(figPath, 'interaction_stability_', predNames[vari], '_', '.png',sep=''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()       
  # }  
}

