





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
    curLM = lm(dv ~ age + sex + diag  +IQ , data = modDat)
    aovTab = Anova(curLM, type = 3)
    curLM$coefficients
    ai = (ii-2)*length(varNames)
    
    results[tt+ai,1] = varNames[tt]
    ei = which(row.names(aovTab) == 'Residuals')
    ri = which(row.names(aovTab) == 'diag')
    results$Diag[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
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
    test = cbind(cos(modDat$dv*(pi/180)), sin(modDat$dv*(pi/180)))
    curLM = lm(test ~ modDat$age + modDat$sex + modDat$diag  + modDat$IQ )
    aovTab = Anova(curLM, type = 3)
    outtests <- car:::print.Anova.mlm
    body(outtests)[[16]] <- quote(invisible(tests))
    body(outtests)[[15]] <- NULL
    tab <- lapply(c("Pillai"), 
                  function(i)  outtests(Anova(curLM, test.statistic=i)))
    tab <- do.call(rbind, tab)
    #working out effect size for diag
    ai = (ii-2)*length(varNames)
    ri = which(grepl('diag', row.names(tab)))
    results[tt+ai,1] = varNames[tt]
    results$Diag[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
                          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
    ri = which(grepl('sex', row.names(tab)))
    results$sex[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
                          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
    ri = which(grepl('IQ', row.names(tab)))
    results$IQ[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
                          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
    ri = which(grepl('age', row.names(tab)))
    results$age[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
                          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
    
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
        curLM = lm(dv ~ age + sex + diag  +IQ , data = modDat)
        aovTab = Anova(curLM, type = 3)
        curLM$coefficients
        ai = (ii-2)*length(varNames)
        
        results[tt+ai,1] = varNames[tt]
        ei = which(row.names(aovTab) == 'Residuals')
        ri = which(row.names(aovTab) == 'diag')
        results$Diag[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
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
        test = cbind(cos(modDat$dv*(pi/180)), sin(modDat$dv*(pi/180)))
        curLM = lm(test ~ modDat$age + modDat$sex + modDat$diag  + modDat$IQ )
        aovTab = Anova(curLM, type = 3)
        outtests <- car:::print.Anova.mlm
        body(outtests)[[16]] <- quote(invisible(tests))
        body(outtests)[[15]] <- NULL
        tab <- lapply(c("Pillai"), 
                      function(i)  outtests(Anova(curLM, test.statistic=i)))
        tab <- do.call(rbind, tab)
        #working out effect size for diag
        ai = (ii-2)*length(varNames)
        ri = which(grepl('diag', row.names(tab)))
        results[tt+ai,1] = varNames[tt]
        results$Diag[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
        ri = which(grepl('sex', row.names(tab)))
        results$sex[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
        ri = which(grepl('IQ', row.names(tab)))
        results$IQ[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
        ri = which(grepl('age', row.names(tab)))
        results$age[tt+ai] = (tab$Df[1] * tab$`approx F`[ri]) / 
          (tab$Df[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
        
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

predNames = c('age', 'sex', 'IQ', 'Diag')
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


for(vari in 1:4){
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
      
    png(paste(figPath, 'effect_', predNames[vari], '_', as.character(agei), '.png',sep=''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }  
}


#### what's the relationship between effect size and effect stability? ####
n = length(resultsDif$age)
allEffects = data.frame('agreement' = c(resultsDif$age, resultsDif$sex, resultsDif$IQ, resultsDif$Diag),
                        'effectSize' = c(resultsTRAIN$age, resultsTRAIN$sex, resultsTRAIN$IQ, resultsTRAIN$Diag),
                        'predictor' = c(rep('age', n), rep('sex', n), rep('IQ', n), rep('diagnosis', n)))

  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'age')[,1],filter(allEffects, predictor == 'age')[,2])
  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'sex')[,1],filter(allEffects, predictor == 'sex')[,2])
  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'IQ')[,1],filter(allEffects, predictor == 'IQ')[,2])
  #how well does stability predict effect size when age is the predictor? 
  cor.test(filter(allEffects, predictor == 'diagnosis')[,1],filter(allEffects, predictor == 'diagnosis')[,2])

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
png(paste(figPath, 'effect_agree_relationship.png', sep = ''),         # File name
    width=400, height=768)
print(outPlot)
dev.off()         

# ggplot(filter(allEffects, predictor == 'sex'), aes(x=effect, y=agreement)) + geom_jitter()


#### agreement between test and train data sets in model fit effect size ####
thresh = .05

for(vari in 1:4){
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
    
    png(paste(figPath, 'stability_', predNames[vari], '_', '.png',sep=''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()       
  # }  
}



#### table 2: what's the number of 'good' predictors as a function of age and predictor? #### 

goodPreds = data.frame('age'= c(0,0,0), 
                       'sex'= c(0,0,0), 
                       'IQ'= c(0,0,0), 
                       'diagnosis'= c(0,0,0),
                       row.names = ageLabs)

#this data frame is used to make supplemental table 1
allVarsPredByDiag = data.frame('dependentVariable' = c(0), 
                               'Diag' = c(0),
                               'diag_b1' = c(0),
                               'diag_b2' = c(0),
                               'age' = c(0))
alli = 1

for(vari in 1:4) { 
  for(agei in 1:3){
    temp = resultsDif[resultsDif$ageGroup == agei,]
    tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
    candidates = which(temp[predNames[vari]] > .8 & tempRes[predNames[vari]] > thresh)
    if(vari==4){
      allVarsPredByDiag[alli:(length(candidates)+alli-1),1:4] = tempRes[candidates,c(1,5,9,10)]
      allVarsPredByDiag[alli:(length(candidates)+alli-1),5] = agei 
      alli = alli + length(candidates)
    }
    goodPreds[agei,vari] = length(candidates)
  }
  
}

goodPreds %>% 
  kbl(align = 'c') %>% 
  kable_classic(full_width = F, 
                font_size = 20) %>%
  column_spec(1, border_right = T)%>%
  footnote(general = "Performance of different predictors in each age group. 
           Number of independent variables predicted with \U1D702\U00B2\U209A > .05 and stability > .80",
           general_title = "Table 2: ",
           footnote_as_chunk = T, title_format = c("italic", "underline")
  )  

chisq.test(goodPreds)


#### choosing the best plot for each 
# datTrain = read.csv("C:\\Users\\Adam Dede\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
# datTrain = read.csv(paste(path, "trainSet.csv", sep = ''))
# datTrain <- datTrain %>% select(-X)
comboDat = rbind(datTrain, dat)
comboDat$ageGroup = 1
comboDat$ageGroup[comboDat$age>ageGroups[2]] = 2
comboDat$ageGroup[comboDat$age>ageGroups[3]] = 3
comboDat$ageGroup = as.factor(comboDat$ageGroup)
#quick check of how many participants had missing data for alpha peak calculation
test = comboDat[,(grepl('logAlpha', names(comboDat)) | grepl('relAlpha', names(comboDat)))]
sum(apply(as.matrix(test), 1, function(x) sum(is.na(x)))> 0) #how many participants had at least one miss? 
sum(apply(as.matrix(test), 1, function(x) sum(is.na(x)))> 5) #how many participants had more than 5 misses? 
1-sum(apply(as.matrix(test), 1, function(x) sum(is.na(x))))/(776*36) #what proportion were well fit? 


## PLOTS FOR AGE! 
vari = 1
  for(agei in 1:3){
    temp = resultsDif[resultsDif$ageGroup == agei,]
    tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
    candidates = which(temp[predNames[vari]] > .8 & tempRes[predNames[vari]] > thresh)
    target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
    targName = temp$dependentVariable[target]
  if(length(target)>0){
    limVals = quantile(comboDat[,targName], c(.1,.9), na.rm = T)
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.2,
                limVals[2] + (limVals[2]-limVals[1])*.2)
  outPlot <- ggplot(comboDat, aes_string(x = 'age', y = targName, color = 'group', shape = 'sex', size = 'IQ')) +
    geom_jitter(alpha = .75)+
    scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                         labels=c(  'AD', 'ASD', 'CON')) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) +
    ylim(limVals) +
    ylab(targName) +
    ggtitle(paste(targName, '; age group: ', as.character(agei), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10)),
           shape = guide_legend(override.aes = list(size=10))) + 
    geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
    geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
  # rect(1, 5, 3, 7, col="white")
  png(paste(figPath, 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  }
  }

## PLOTS FOR SEX! 
vari = 2
for(agei in 1:3){
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  candidates = which(temp[predNames[vari]] > .8 & tempRes[predNames[vari]] > thresh)
  target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
  targName = temp$dependentVariable[target]
  if(length(target)>0){
    limVals = quantile(comboDat[,targName], c(.1,.9), na.rm = T)
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.2,
                limVals[2] + (limVals[2]-limVals[1])*.2)
    outPlot <-  comboDat %>% 
      ggplot(aes_string(x = 'ageGroup', y = targName, fill = 'sex')) +
      # geom_violin(draw_quantiles = c(.25, .5, .75), alpha = .5) + 
      geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.9)) + #width = .6, coef = 0
      geom_jitter(show.legend = F, shape = 21, size = 4,
                  position = position_jitterdodge(dodge.width = .9)) +
      # stat_summary(fun = median, geom = 'crossbar' ) +
      # geom_dotplot(binaxis = 'y', show.legend = F, stackdir = 'center',
      #              position=position_dodge(1), dotsize = .30)+
      theme_classic() +
      theme(axis.line = element_line(color = 'black', size = 3),
            axis.ticks = element_line(colour = "black", size = 2),
            axis.ticks.length=unit(-.25, "cm"),
            text = element_text(size = 20)) +
      ylim(limVals) +
      ylab(targName) +
      ggtitle(paste(targName, '; age group: ', as.character(agei), sep = ''))+
      guides(color = guide_legend(override.aes = list(size=10)),
             shape = guide_legend(override.aes = list(size=10)))  
      # geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
      # geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
    # rect(1, 5, 3, 7, col="white")
    png(paste(figPath, 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }
}

## PLOTS FOR IQ! 
# vari = 3
# for(agei in 1:3){
#   temp = resultsDif[resultsDif$ageGroup == agei,]
#   tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
#   candidates = which(temp[predNames[vari]] < .2 & tempRes[predNames[vari]] > thresh)
#   target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
#   if(length(target)>0){
#     limVals = quantile(comboDat[,varNames[target]], c(.1,.9))
#     limVals = c(limVals[1] - (limVals[2]-limVals[1])*.2,
#                 limVals[2] + (limVals[2]-limVals[1])*.2)
#     outPlot <- ggplot(comboDat, aes_string(x = 'IQ', y = varNames[target], color = 'age', shape = 'group')) +
#       geom_jitter(alpha = .75, size = 7)+
#       # scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
#       #                    labels=c(  'AD', 'ASD', 'CON')) +
#       theme_classic() +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) +
#       ylim(limVals) +
#       ylab(varNames[target]) +
#       ggtitle(paste(varNames[target], '; age group: ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10)),
#              shape = guide_legend(override.aes = list(size=10)))  
#       # geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
#       # geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
#     # rect(1, 5, 3, 7, col="white")
#     png(paste( 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#   }
# }

## PLOTS FOR Diag! 
vari = 4
for(agei in 1:3){
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  candidates = which(temp[predNames[vari]] > .8 & tempRes[predNames[vari]] > thresh)
  # if(agei==1){
  #   target = candidates[6]
  # }else{
    target = candidates[which(tempRes[candidates,'Diag']==max(tempRes[candidates,'Diag']))]
  # }
  targName = temp$dependentVariable[target]
 if(length(target)>0){
    limVals = quantile(comboDat[,targName], c(.1,.9), na.rm = T)
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
                limVals[2] + (limVals[2]-limVals[1])*.5)
    # plotDat <- comboDat %>% filter( group %in% c('CON', 'AD')) %>% arrange(ageGroup, group)
    # n = length(plotDat$age)
    # plotDat[(n+1):(n+3), varNames[target]] = -1000000
    # plotDat[(n+1):(n+3), 'ageGroup'] = c(1,2,3)
    # plotDat[(n+1):(n+3), 'group'] = c('CON', 'CON', 'CON')
    # plotDat$ageGroup = as.factor(plotDat$ageGroup)
    outPlot <-  comboDat %>% #age> ageGroups[agei], age<= ageGroups[agei+1],
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
            text = element_text(size = 20)) +
      ylim(limVals) +
      ylab(targName) +
      ggtitle(paste(targName, '; age group: ', as.character(agei), sep = ''))+
      guides(color = guide_legend(override.aes = list(size=10))) 
    # rect(1, 5, 3, 7, col="white")
    png(paste(figPath, 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }
}

########################################################################
#### signal detection theory based analysis of potential biomarkers ####
########################################################################

SDTTable = data.frame('ageGroup' = c('3-96 mon', '98-126 mon', '128-248 mon'),
                      'AUC'      = c(0,0,0), 
                      'accuracy' = c(0,0,0), 
                      'sensitivity'=c(0,0,0),
                      'specificity'=c(0,0,0),
                      'IG_.0148' = c(0,0,0),
                      'IG_.187' = c(0,0,0),
                      'FPi_.0148' = c(0,0,0),
                      'FPi_.187' = c(0,0,0))

plotValues = data.frame('hitRate' = rep(0,300),
                        'FArate' = rep(0,300), 
                        'precision'=rep(0,300),
                        'infoGain' = rep(0,300),
                        'FPi' = rep(0,300),
                        'ageGroup' = rep(0,300))
for(agei in 1:3){
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  tempTrain = datTrain %>% filter(age>ageGroups[agei], age <=ageGroups[agei+1])
  tempTest = dat %>% filter(age>ageGroups[agei], age <= ageGroups[agei+1])
  # tempTest = dat %>% filter(age>ageGroups[agei], age <=ageGroups[agei+1])
  candidates = which(temp$Diag > .8 & tempRes$Diag > thresh)

  target = candidates[which(tempRes[candidates,'Diag']==max(tempRes[candidates,'Diag']))]
  
  targName = temp$dependentVariable[target]
  limVals = quantile(tempTrain[,targName], c(.025,.975), na.rm = T)
  # limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
  #             limVals[2] + (limVals[2]-limVals[1])*.5)
  

  # make edits here for AD vs. CON OR AD+ASD vs. CON
  groupID = rep(1,length(tempTrain$group)) # 2 = 'control'
  groupID[tempTrain$group != 'CON'] = 2 # 1 = 'patient'
  # groupID[tempTrain$group == 'ASD'] = NA
  
  #ALL AD AND ASD SUBS CODED AS 2
  biomark = tempTrain[,targName]
  missingDat = which(is.na(biomark))
  if(length(missingDat)>0){
    groupID = groupID[-missingDat]
    biomark = biomark[-missingDat]
  }
  missingDat = which(is.na(groupID))
  if(length(missingDat)>0) {
    groupID = groupID[-missingDat]
    biomark = biomark[-missingDat]
  }
  
  #get the variables for the test set: 
  groupID2 = rep(1,length(tempTest$group)) # 2 = 'control'
  groupID2[tempTest$group != 'CON'] = 2 # 2 = 'patient'
  # groupID2[tempTest$group == 'ASD'] = NA
  
  #ALL AD AND ASD SUBS CODED AS 2
  biomark2 = tempTest[,targName]
  missingDat2 = which(is.na(biomark2))
  if(length(missingDat2)>0){
    groupID2 = groupID2[-missingDat2]
    biomark2 = biomark2[-missingDat2]
  }
  missingDat2 = which(is.na(groupID2))
  if(length(missingDat2)>0) {
    groupID2 = groupID2[-missingDat2]
    biomark2 = biomark2[-missingDat2]
  }
  
  criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
  criteria = criteria[1:100]
  # #is the AD group higher or lower than the CON group? 
  # #default is autistic > control
  if(mean(biomark[groupID==1]) > mean(biomark[groupID==2])){
    print("POTENTIAL PROBLEM")
  } 
  
  #inverted criterion logic
  ROC = getROC(biomark, groupID, criteria)
  ROC2 = getROC(biomark2, groupID2, criteria)
  
  # ROC2 = getROC(biomark, groupID, criteria)
  ci = which(ROC$acc == max(ROC$acc))
  ci = ci[1]
  SDTTable$AUC[agei] = paste(round(ROC$AUC[1],2), ' (', round(ROC2$AUC[1],2), ')', sep = '') 
  SDTTable$accuracy[agei] = paste(round(ROC$acc[ci],2), ' (', round(ROC2$acc[ci],2), ')', sep = '')
  SDTTable$sensitivity[agei] = paste(round(ROC$TPR[ci],2), ' (', round(ROC2$TPR[ci],2), ')', sep = '')
  SDTTable$specificity[agei] = paste(round(ROC$TNR[ci],2), ' (', round(ROC2$TNR[ci],2), ')', sep = '')


  
  
  ROC[101:102,] = c(100,100)
  ROC$acc[101:102] = c(0,1)
  # outPlot <- ggplot(ROC, aes(x=FP, y=TP, size = log10(acc))) + 
  #   geom_point(show.legend = T, alpha = .75, color = '#228B22') +
  #   scale_size_continuous(range = c(.25, 45)) + 
  #   geom_abline(slope = 1, intercept = 0, size = 2, linetype = 'dashed') + 
  #   ylim(c(0,1)) + 
  #   xlim(c(0,1)) +
  #   theme_classic() +
  #   theme(axis.line = element_line(color = 'black', size = 3),
  #         axis.ticks = element_line(colour = "black", size = 2),
  #         axis.ticks.length=unit(-.25, "cm"),
  #         text = element_text(size = 20),
  #         panel.grid.major = element_line(color = "grey",
  #                                         size = 2,
  #                                         linetype = 2)) +
  #   ylab('hit rate') +
  #   xlab('false alarm rate') + 
  #   ggtitle(paste('ROC; age group: ', as.character(agei), sep = ''))
  # png(paste( 'bestPredicted_ROC',predNames[vari],'_', agei, '.png', sep = ''),         # File name
  #     width=1024, height=768)
  # print(outPlot)
  # dev.off()
  
 
  
  #simulating various base rates
  brDF = data.frame('base' = c(seq(.001, 1, .01),.0148,.187),
                    'infoGain' = c(seq(.001, 1, .01),.0148,.187),
                    'FPi' = c(seq(.001, 1, .01),.0148,.187),
                    'PPV' = c(seq(.001, 1, .01),.0148,.187),
                    'post' = c(seq(.001, 1, .01),.0148,.187))
 
  for(bi in 1:length(brDF$base)){
    print(bi)
    baseRate = brDF$base[bi]
    ADGroup = sample(biomark2[groupID2==2], round(baseRate*100000), replace = T )
    CONGroup = sample(biomark2[groupID2==1], round((1-baseRate)*100000), replace = T)
    simBiomark = c(ADGroup, CONGroup)
    simID = rep(1,length(simBiomark))
    simID[1:length(ADGroup)] = 2
    simROC = getROC(simBiomark, simID, criteria)
    post = simROC$TP[ci]*baseRate / (simROC$TP[ci]*baseRate + simROC$FP[ci]*(1-baseRate))
   
    brDF$infoGain[bi] = post - baseRate
    
    guess = rep(1,length(simBiomark))
    guess[simBiomark>criteria[ci]] = 2
    
    
    brDF$FPi[bi] = sum(guess==2 & simID==1) / sum(guess==2 & simID==2)
    
    brDF$PPV[bi] = simROC$hits[ci] / (simROC$hits[ci] + simROC$FAs[ci])
    brDF$post[bi] = post
    
  }
  
  SDTTable$IG_.0148[agei] = brDF$infoGain[101]
  SDTTable$IG_.187[agei] = brDF$infoGain[102]
  SDTTable$FPi_.0148[agei] = brDF$FPi[101]
  SDTTable$FPi_.187[agei] = brDF$FPi[102]
  
  plotValues$hitRate[((agei-1)*100+1):((agei*100))] = ROC$TP[1:100]
  plotValues$FArate[((agei-1)*100+1):(agei*100)] = ROC$FP[1:100]
  plotValues$precision[((agei-1)*100+1):(agei*100)] = brDF$post[1:100]
  plotValues$infoGain[((agei-1)*100+1):(agei*100)] = brDF$infoGain[1:100]
  plotValues$FPi[((agei-1)*100+1):(agei*100)] = brDF$FPi[1:100]
  plotValues$ageGroup[((agei-1)*100+1):(agei*100)] = agei
  
  
  outPlot <- ggplot(brDF, aes(x=base, y = post)) + 
    geom_line(color = '#66023C', size = 5) + 
    geom_abline(slope = 1, intercept = 0, size = 2, linetype = 'dashed') + 
    theme_classic() + 
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20),
          panel.grid.major = element_line(color = "grey",
                                          size = 2,
                                          linetype = 2)) +
    ylab('precision') +
    ggtitle(paste('posterior probability; age group: ', as.character(agei), sep = ''))
  # png(paste( 'bestPredicted_precision',predNames[vari],'_', agei, '.png', sep = ''),         # File name
  #     width=1024, height=768)
  # print(outPlot)
  # dev.off()
  
  
  outPlot <- ggplot(brDF, aes(x=base, y = infoGain)) + 
    geom_point(color = '#66023C', size = 5) + 
    scale_y_log10() +
    theme_classic() + 
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20),
          panel.grid.major = element_line(color = "grey",
                                          size = 2,
                                          linetype = 2)) +
    ylim(c(0,.25)) + 
    ylab('information gain') +
    ggtitle(paste('information gain; age group: ', as.character(agei), sep = ''))
  # png(paste( 'bestPredicted_infoGain',predNames[vari],'_', agei, '.png', sep = ''),         # File name
  #     width=1024, height=768)
  # print(outPlot)
  # dev.off()
  
  label_at <- function(n) function(x) ifelse(x %% n == 0, x, "")
  breakPoints =  c(seq(.02,.1,.02), seq(.2,1,.2), seq(2,10,2), seq(20,100,20),seq(200,1000,200))
  
  outPlot <- ggplot(brDF, aes(x=base, y = FPi)) + 
    geom_line(color = '#66023C', size = 5) + 
    scale_y_log10(breaks = breakPoints, labels = c(rep('',4), .1, 
                                                   rep('',4), 1,
                                                   rep('',4), 10,
                                                   rep('',4), 100,
                                                   rep('',4), 1000)) +
    theme_classic() + 
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20),
          panel.grid.major = element_line(color = "grey",
                                          size = 2,
                                          linetype = 2)) +
    ylab(targName) +
    xlim(c(0,.5))+ 
    ggtitle(paste('false positive index; age group: ', as.character(agei), sep = ''))
  # png(paste( 'bestPredicted_FPi',predNames[vari],'_', agei, '.png', sep = ''),         # File name
  #     width=1024, height=768)
  # print(outPlot)
  # dev.off()
  
}
plotValues$baseRate = rep(brDF$base[1:100], 3)

SDTTable %>% 
  kbl(align = 'c', digits = 2) %>% 
  kable_classic(full_width = F, 
                font_size = 20) %>%
  column_spec(1, border_right = T) %>%
  footnote(general = "Signal detection theory measures of biomarker performance. 
           AUC = area under the receiver operating characteristic curve
           IG_.0148 = information gain when base rate is .0148
           FPi_.0148 = number of false positives for each true positive when base rate is .0148
           IG_.187 and FPi_.187 = corresponding values for a base rate of .187",
           general_title = "Table 3: ",
           footnote_as_chunk = T, title_format = c("italic", "underline")
  )  

##### plotting signal detection theory measure results ####
plotValues$ageGroup = as.factor(plotValues$ageGroup)
plotValues$ageGroup = factor(plotValues$ageGroup)
outPlot <- ggplot(plotValues, aes(y=hitRate, x=FArate, color = ageGroup)) + 
  geom_line(show.legend = T, size = 5) + 
  geom_abline(slope = 1, intercept = 0, size = 2, linetype = 'dashed') + 
  ylim(c(0,1)) + 
  xlim(c(0,1)) +
  theme_classic() + 
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ylab('hit rate') +
  xlab('false alarm rate') + 
  ggtitle(paste('ROC', sep = ''))
png(paste(figPath, 'ROC_all',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()


outPlot <- ggplot(plotValues, aes(x=baseRate, y = precision, color = ageGroup)) + 
  geom_line(size = 5) + 
  geom_abline(slope = 1, intercept = 0, size = 2, linetype = 'dashed') + 
  geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  theme_classic() + 
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ylab('precision') +
  ggtitle(paste('posterior probability', sep = ''))
png(paste(figPath, 'precision_all',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

outPlot <- ggplot(plotValues, aes(x=baseRate, y = infoGain, color = ageGroup)) + 
  geom_point(size = 5) + 
  geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  theme_classic() + 
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ylim(c(0,.3)) + 
  ylab('information gain') +
  ggtitle(paste('information gain', sep = ''))
png(paste(figPath, 'infoGain_all',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()
breakPoints =  c(.1, 1, 10, 100, 1000)

outPlot <- ggplot(plotValues, aes(x=baseRate, y = FPi, color = ageGroup)) + 
  geom_line(size = 5) + 
  geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
  scale_y_log10(breaks = breakPoints, limits = c(.1,1000)) +
  theme_classic() + 
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ylab('count FP / count TP') +
  xlim(c(0,.5))+
  ggtitle(paste('false positive index', sep = ''))
png(paste(figPath, 'FPi_all',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

############## simulating at specific base rates ################

#simulating various base rates
brDF = data.frame('base' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3) ,
                  'infoGain' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3),
                  'FPi' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3),
                  'PPV' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3),
                  'post' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3),
                  'AUC' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3),
                  'acc' = rep(c(rep(.0148,100), rep(.187, 100), rep(.5,100)),3),
                  'ageGroup' = c(rep(1,300), rep(2,300), rep(3,300)))

for(agei in 1:3){
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  candidates = which(temp$Diag > .8 & tempRes$Diag > thresh)
  if(agei==1){
    target = candidates[6]
  }else{
    target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
  }
  targName = temp$dependentVariable[target]
  
  # limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
  #             limVals[2] + (limVals[2]-limVals[1])*.5)
  
  curDat = comboDat %>% filter(ageGroup==agei)
  groupID = rep(1,length(curDat$group))
  groupID[curDat$group != 'CON'] = 2
  #ALL AD AND ASD SUBS CODED AS 2
  biomark = curDat[,targName]
 
  #is the AD group higher or lower than the CON group? 
  #default is autistic < control
  if(agei==1){
    biomark = biomark * -1
  }
  limVals = quantile(biomark, c(.025,.975))
  criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
  criteria = criteria[1:100]
  
  for(bi in 1:300){
    print((agei-1)*300+bi)
    baseRate = brDF$base[bi]
    ADGroup = sample(biomark[groupID==2], round(baseRate*100000), replace = T )
    CONGroup = sample(biomark[groupID==1], round((1-baseRate)*100000), replace = T)
    simBiomark = c(ADGroup, CONGroup)
    simID = rep(1,length(simBiomark))
    simID[1:length(ADGroup)] = 2
    simROC = getROC(simBiomark, simID, criteria)
    brDF$infoGain[(agei-1)*300+bi] = simROC$post[1] - simROC$base[1]
    brDF$FPi[(agei-1)*300+bi] = simROC$FPi[1]
    brDF$PPV[(agei-1)*300+bi] = simROC$PPV[1]
    brDF$post[(agei-1)*300+bi] = simROC$post[1]
    brDF$AUC[(agei-1)*300+bi] = simROC$AUC[1] 
    brDF$acc[(agei-1)*300+bi] = max(simROC$accRaw) 
  }
 
  
  


}

brDF$base = as.factor(brDF$base)
brDF$ageGroup = as.factor(brDF$ageGroup)
outPlot <- brDF %>% #age> ageGroups[agei], age<= ageGroups[agei+1],
  ggplot(aes(x = base, y = infoGain, fill = ageGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.75)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4,
              position = position_jitterdodge(dodge.width = .75, jitter.width = .4)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ylim(c(0,.25)) + 
  ggtitle(paste('infoGain bootstrapped', sep = ''))+
  guides(color = guide_legend(override.aes = list(size=10))) 
png(paste( 'infoGainBoot',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

infoGainLM = lm(infoGain~base*ageGroup, brDF)
Anova(infoGainLM)

breakPoints =  c(.1, 1, 10, 100, 1000)
outPlot <- brDF %>% #age> ageGroups[agei], age<= ageGroups[agei+1],
  ggplot(aes(x = base, y = FPi, fill = ageGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.75)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4,
              position = position_jitterdodge(dodge.width = .75, jitter.width = .4)) +
  theme_classic() +
  scale_y_log10(breaks = breakPoints, limits = c(.1,1000)) +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ggtitle(paste('Fpi bootstrapped', sep = ''))+
  guides(color = guide_legend(override.aes = list(size=10))) 
png(paste( 'FPiBoot',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()
fpiLM = lm(FPi~base*ageGroup, brDF)
Anova(fpiLM)

outPlot <- brDF %>% #age> ageGroups[agei], age<= ageGroups[agei+1],
  ggplot(aes(x = base, y = AUC, fill = ageGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.75)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4,
              position = position_jitterdodge(dodge.width = .75, jitter.width = .4)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ggtitle(paste('AUC bootstrapped', sep = ''))+
  guides(color = guide_legend(override.aes = list(size=10))) 
png(paste( 'AUCBoot',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()
aucLM = lm(AUC~base*ageGroup, brDF)
Anova(aucLM)



outPlot <- brDF %>% #age> ageGroups[agei], age<= ageGroups[agei+1],
  ggplot(aes(x = base, y = acc, fill = ageGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, position = position_dodge(.75)) + #width = .6, coef = 0
  geom_jitter(show.legend = F, shape = 21, size = 4,
              position = position_jitterdodge(dodge.width = .75, jitter.width = .4)) +
  theme_classic() +
  ylim(c(.5,1.0)) + 
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey",
                                        size = 2,
                                        linetype = 2)) +
  ggtitle(paste('acc bootstrapped', sep = ''))+
  guides(color = guide_legend(override.aes = list(size=10))) 
png(paste( 'ACCBoot',predNames[vari],'_', agei, '.png', sep = ''),         # File name
    width=1024, height=768)
print(outPlot)
dev.off()
accLM = lm(acc~base*ageGroup, brDF)
Anova(accLM)


##################################################################################
#############################  assessing cost/benefit of biomarkers ##############
##################################################################################

baseVals = c(.0148, .187)
costRatio = logspace(-3,3,100)
for(agei in 1:3){
  

  #   limVals = quantile(comboDat[,targName], c(.025,.975))
  #   # limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
  #   #             limVals[2] + (limVals[2]-limVals[1])*.5)
  #   
  #   curDat = comboDat %>% filter(ageGroup==agei)
  
  #   biomark = curDat[,targName]
  #   criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
  #   criteria = criteria[1:100]
  #   #is the AD group higher or lower than the CON group? 
  #   #default is autistic < control

  
  
  
    temp = resultsDif[resultsDif$ageGroup == agei,]
    tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
    candidates = which(temp$Diag > .8 & tempRes$Diag > thresh)
    if(agei==1){
      target = candidates[6]
    }else{
      target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
    }
    targName = temp$dependentVariable[target]
 
  
  curDat = comboDat %>% filter(ageGroup==agei)
  groupID = rep(1,length(curDat$group))
  groupID[curDat$group != 'CON'] = 2
  #ALL AD AND ASD SUBS CODED AS 2
  biomark = curDat[,targName]
  
  #is the AD group higher or lower than the CON group? 
  #default is autistic > control
  if(agei==1){
    biomark = biomark * -1
  }
  limVals = quantile(biomark, c(.025,.975))
  criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
  criteria = criteria[1:100]
  
  ROC = getROC(biomark, groupID, criteria)
  
  
  ### prep a results data frame to try various base rates
  baseRates = c(seq(.001, 1, .01))
  results = expand.grid(base=baseRates, costRatio=costRatio)
  results$TPR = NaN
  results$FPR = NaN
  results$acc = NaN
  results$PPV = NaN
  results$infoGain = NaN
  
  
  for(bb in 1:length(baseRates)){
    print(bb)
    baseRate = baseRates[bb]
    resi = which(results$base == baseRate)
    ADGroup = sample(biomark[groupID==2], round(baseRate*1000000), replace = T )
    CONGroup = sample(biomark[groupID==1], round((1-baseRate)*1000000), replace = T)
    simBiomark = c(ADGroup, CONGroup)
    simID = rep(1,length(simBiomark))
    simID[1:length(ADGroup)] = 2
    simROC = getROC(simBiomark, simID, criteria)
    simROC$TPR[simROC$TPR==0] = 1/1000000
    simROC$FPR[simROC$FPR==0] = 1/1000000
    
    for(ii in 1:length(resi)){
      
      ratio = simROC$TPR/simROC$FPR
      costVal = ((1-baseRate)/baseRate) * results$costRatio[resi[ii]]
      
      critVal = which(abs(ratio-costVal)==min(abs(ratio-costVal)))
      if(length(critVal)>1){
        critVal = critVal[1]
      }
      if(length(critVal)==1){
        if(is.finite(critVal)){
          results$TPR[resi[ii]] = simROC$TPR[critVal]
          results$FPR[resi[ii]] = simROC$FPR[critVal]
          results$acc[resi[ii]] = simROC$acc[critVal]
          results$PPV[resi[ii]] = simROC$hits[critVal] / (simROC$hits[critVal] + simROC$FAs[critVal])
          results$infoGain[resi[ii]] = ROC$TP[critVal]*baseRate / (ROC$TP[critVal]*baseRate + ROC$FP[critVal]*(1-baseRate))
        }
      }
    }
  }
  
  breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
  outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = TPR-FPR)) + 
    geom_tile() + 
    scale_y_log10(breaks = breakPoints) +
    scale_fill_gradient(low = '#FFFFFF', high = '#3A317D', na.value = 'white') + 
    theme_classic() +
    geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
    geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
    geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) +
    ylab('FP<FN      costRatio     FP>FN') + 
    xlab('base rate') +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) + 
    ggtitle(paste('TPR-FPR as a function of costs and benefits ', as.character(agei), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10))) 
  
  locs = seq(-3.5,3.5,1)  
  for(vv in 1:8){
    outPlot = outPlot + 
      geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  locs = c(0,.25,.5,.75,1)
  for(vv in 1:length(locs)){
    outPlot = outPlot + 
      geom_vline(xintercept = locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  png(paste( 'TPR_FPR_costmap_age', agei, '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
  
  outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = FPR)) + 
    geom_tile() + 
    scale_y_log10(breaks = breakPoints) +
    scale_fill_gradient(low = '#4D7902', high = '#CE2029', na.value = 'black',
                        limits = c(0, 1)) + 
    theme_classic() +
    ylab('FP<FN      costRatio     FP>FN') + 
    xlab('base rate') +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) + 
    ggtitle(paste('FPR as a function of costs and benefits ', as.character(agei), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10))) 
  
  locs = seq(-3.5,3.5,1)  
  for(vv in 1:8){
    outPlot = outPlot + 
      geom_vline(xintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2) + 
      geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  
  png(paste( 'FPR_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
  # results$accRound = results$acc
  # results$accRound[results$acc>=.60] = .59
  # results$accRound[results$acc<=.50] = .50
  
  breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
  outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = acc)) + 
    geom_tile() + 
    scale_y_log10(breaks = breakPoints) +
    scale_fill_gradient(low = '#FFFFFF', high = '#3A317D', na.value = 'black')+
    #limits = c(.5, .75)) + 
    theme_classic() +
    ylab('FP<FN      costRatio     FP>FN') + 
    xlab('base rate') +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) + 
    geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
    geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
    geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) +
    ggtitle(paste('accuracy as a function of costs and benefits ', as.character(agei), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10))) 
  
  locs = seq(-3.5,3.5,1)  
  for(vv in 1:8){
    outPlot = outPlot + 
      geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  locs = c(0,.25,.5,.75,1)
  for(vv in 1:length(locs)){
    outPlot = outPlot + 
      geom_vline(xintercept = locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  print(outPlot)
  png(paste( 'acc_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
  outPlot <- ggplot(results, aes(x = acc)) + 
    geom_histogram()
  png(paste( 'acc_hist_age', agei, 'base', bb, '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
  
  ggplot(results, aes(x = FPR, y = TPR)) + geom_jitter()
  
  
}






##################################################################################
#############################  assessing cost/benefit of biomarkers of SEX #######
##################################################################################
baseVals = c(.0148, .187)
costRatio = logspace(-3,3,100)
costRatio = logspace(-3,3,100)
for(agei in 1:3){
  
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  candidates = which(temp$sex > .8 & tempRes$sex > thresh)
  target = candidates[which.max(tempRes[candidates,predNames[vari]])]
  targName = temp$dependentVariable[target]
  
  # limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
  #             limVals[2] + (limVals[2]-limVals[1])*.5)
  
  curDat = comboDat %>% filter(ageGroup==agei)
  groupID = rep(1,length(curDat$sex))
  groupID[curDat$sex != 'M'] = 2
  #ALL Female SUBS CODED AS 2
  biomark = curDat[,targName]
  
  #is the AD group higher or lower than the CON group? 
  #default is autistic > control
  if(agei>1){
    biomark = biomark * -1
  }
  limVals = quantile(biomark, c(.025,.975))
  criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
  criteria = criteria[1:100]
  
  ROC = getROC(biomark, groupID, criteria)
  
  
  ### prep a results data frame to try various base rates
  baseRates = c(seq(.001, 1, .01))
  results = expand.grid(base=baseRates, costRatio=costRatio)
  results$TPR = NaN
  results$FPR = NaN
  results$acc = NaN
  
  
  for(bb in 1:length(baseRates)){
    print(bb)
    baseRate = baseRates[bb]
    resi = which(results$base == baseRate)
    ADGroup = sample(biomark[groupID==2], round(baseRate*1000000), replace = T )
    CONGroup = sample(biomark[groupID==1], round((1-baseRate)*1000000), replace = T)
    simBiomark = c(ADGroup, CONGroup)
    simID = rep(1,length(simBiomark))
    simID[1:length(ADGroup)] = 2
    simROC = getROC(simBiomark, simID, criteria)
    simROC$TPR[simROC$TPR==0] = 1/1000000
    simROC$FPR[simROC$FPR==0] = 1/1000000
    
    for(ii in 1:length(resi)){
      
      ratio = simROC$TPR/simROC$FPR
      costVal = ((1-baseRate)/baseRate) * results$costRatio[resi[ii]]
      
      critVal = which(abs(ratio-costVal)==min(abs(ratio-costVal)))
      if(length(critVal)>1){
        critVal = critVal[1]
      }
      if(length(critVal)==1){
        if(is.finite(critVal)){
          results$TPR[resi[ii]] = simROC$TPR[critVal]
          results$FPR[resi[ii]] = simROC$FPR[critVal]
          results$acc[resi[ii]] = simROC$acc[critVal]
          results$PPV[resi[ii]] = simROC$hits[critVal] / (simROC$hits[critVal] + simROC$FAs[critVal])
          results$infoGain[resi[ii]] = ROC$TP[critVal]*baseRate / (ROC$TP[critVal]*baseRate + ROC$FP[critVal]*(1-baseRate))
        }
      }
    }
  }
  
  breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
  outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = TPR-FPR)) + 
    geom_tile() + 
    scale_y_log10(breaks = breakPoints) +
    scale_fill_gradient(low = '#FFFFFF', high = '#3A317D', na.value = 'white') + 
    theme_classic() +
    geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
    geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
    geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) +
    ylab('FP<FN      costRatio     FP>FN') + 
    xlab('base rate') +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) + 
    ggtitle(paste('TPR-FPR as a function of costs and benefits ', as.character(agei), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10))) 
  
  locs = seq(-3.5,3.5,1)  
  for(vv in 1:8){
    outPlot = outPlot + 
      geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  locs = c(0,.25,.5,.75,1)
  for(vv in 1:length(locs)){
    outPlot = outPlot + 
      geom_vline(xintercept = locs[vv], color = 'grey', size = 2, linetype =2)
    
  }
  png(paste( 'TPR_FPR_SEX_age', agei, '.png', sep = ''),         # File name
      width=1024, height=768)
  print(outPlot)
  dev.off()
  
  
}



# 
# 
# for(agei in 1:3){
#   temp = resultsDif[resultsDif$ageGroup == agei,]
#   tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
#   candidates = which(temp$sex > .8 & tempRes$sex > thresh)
#   target = candidates[which.max(tempRes[candidates,predNames[vari]])]
#   targName = temp$dependentVariable[target]
#   
#   # limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
#   #             limVals[2] + (limVals[2]-limVals[1])*.5)
#   
#   curDat = comboDat %>% filter(ageGroup==agei)
#   groupID = rep(1,length(curDat$sex))
#   groupID[curDat$sex != 'M'] = 2
#   #ALL Female SUBS CODED AS 2
#   biomark = curDat[,targName]
#   
#   #is the AD group higher or lower than the CON group? 
#   #default is autistic > control
#   if(agei>1){
#     biomark = biomark * -1
#   }
#   limVals = quantile(biomark, c(.025,.975))
#   criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
#   criteria = criteria[1:100]
#   
#   ROC = getROC(biomark, groupID, criteria)
#   
#   
#   ### prep a results data frame to try various base rates
#   baseRates = c(seq(.001, 1, .01))
#   results = expand.grid(base=baseRates, costRatio=costRatio)
#   results$TPR = NaN
#   results$FPR = NaN
#   results$acc = NaN
#   
#   
#   for(bb in 1:length(baseRates)){
#     print(bb)
#     baseRate = baseRates[bb]
#     resi = which(results$base == baseRate)
#     ADGroup = sample(biomark[groupID==2], round(baseRate*10000), replace = T )
#     CONGroup = sample(biomark[groupID==1], round((1-baseRate)*10000), replace = T)
#     simBiomark = c(ADGroup, CONGroup)
#     simID = rep(1,length(simBiomark))
#     simID[1:length(ADGroup)] = 2
#     simROC = getROC(simBiomark, simID, criteria)
#     simROC$TPR[simROC$TPR==0] = 1/10000
#     simROC$FPR[simROC$FPR==0] = 1/10000
# 
#     for(ii in 1:length(resi)){
#       
#       ratio = simROC$TPR/simROC$FPR
#       costVal = ((1-baseRate)/baseRate) * results$costRatio[resi[ii]]
#       
#       critVal = which(abs(ratio-costVal)==min(abs(ratio-costVal)))
#       if(length(critVal)>1){
#         critVal = critVal[1]
#       }
#       if(length(critVal)==1){
#         if(is.finite(critVal)){
#           results$TPR[resi[ii]] = simROC$TPR[critVal]
#           results$FPR[resi[ii]] = simROC$FPR[critVal]
#           results$acc[resi[ii]] = simROC$acc[critVal]
#           
#         }
#       }
#     }
#   }
#     breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
#     outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = TPR/FPR)) + 
#       geom_tile() + 
#       scale_y_log10(breaks = breakPoints) +
#       scale_fill_gradient(low = '#CE2029', high = '#4D7902', na.value = 'black')+
#                          # limits = c(0, 1)) + 
#       theme_classic() +
#       ylab('FP<FN      costRatio     FP>FN') + 
#       xlab('base rate') +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) + 
#       ggtitle(paste('TPR as a function of costs and benefits ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10))) 
#     
#     locs = seq(-3.5,3.5,1)  
#     for(vv in 1:8){
#       outPlot = outPlot + 
#         geom_vline(xintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2) + 
#         geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     png(paste( 'TPR_SEX_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     
#     outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = FPR)) + 
#       geom_tile() + 
#       scale_y_log10(breaks = breakPoints) +
#       scale_fill_gradient(low = '#4D7902', high = '#CE2029', na.value = 'black',
#                           limits = c(0, 1)) + 
#       theme_classic() +
#       ylab('FP<FN      costRatio     FP>FN') + 
#       xlab('base rate') +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) + 
#       ggtitle(paste('FPR as a function of costs and benefits ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10))) 
#     
#     locs = seq(-3.5,3.5,1)  
#     for(vv in 1:8){
#       outPlot = outPlot + 
#         geom_vline(xintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2) + 
#         geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     
#     png(paste( 'FPR_SEX_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     # results$accRound = results$acc
#     # results$accRound[results$acc>=.60] = .59
#     # results$accRound[results$acc<=.50] = .50
#     
#     breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
#     outPlot <- ggplot(results, aes(x=base, y = costRatio, fill = acc)) + 
#       geom_tile() + 
#       scale_y_log10(breaks = breakPoints) +
#       scale_fill_gradient(low = '#FFFFFF', high = '#3A317D', na.value = 'black')+
#                           #limits = c(.5, .75)) + 
#       theme_classic() +
#       ylab('FP<FN      costRatio     FP>FN') + 
#       xlab('base rate') +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) + 
#       geom_vline(xintercept = .0148, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
#       geom_vline(xintercept = .187, size = 2, linetype = 'longdash', color = 'red', alpha = .5) + 
#       geom_vline(xintercept = .5, size = 2, linetype = 'longdash', color = 'red', alpha = .5) +
#       ggtitle(paste('accuracy as a function of costs and benefits ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10))) 
#     
#     locs = seq(-3.5,3.5,1)  
#     for(vv in 1:8){
#       outPlot = outPlot + 
#         geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     locs = c(0,.25,.5,.75,1)
#     for(vv in 1:length(locs)){
#       outPlot = outPlot + 
#         geom_vline(xintercept = locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     print(outPlot)
#     png(paste( 'acc_SEX_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     outPlot <- ggplot(results, aes(x = acc)) + 
#       geom_histogram()
#     png(paste( 'acc_SEX_hist_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     
#     ggplot(results, aes(x = FPR, y = TPR)) + geom_jitter()
#     
#     
#   }
#   
  
  



#####################




setwd("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\figures\\summaryFigs")
outPlot<-ggplot(filter(results, ageGroup>0), aes(x=Diag, color = factor(ageGroup), fill = factor(ageGroup) )) + geom_histogram(binwidth = .001) +
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50)) 
png('allDiagEffectSize.png',         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

outPlot<-ggplot(filter(results, ageGroup>0), aes(x=age, color = factor(ageGroup), fill = factor(ageGroup) )) + geom_histogram(binwidth = .001) +
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50)) 
png('allAgeEffectSize.png',         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

outPlot<-ggplot(filter(results, ageGroup>0), aes(x=sex, color = factor(ageGroup), fill = factor(ageGroup) )) + geom_histogram(binwidth = .001) +
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50)) 
png('allSexEffectSize.png',         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

outPlot<-ggplot(filter(results, ageGroup>0), aes(x=IQ, color = factor(ageGroup), fill = factor(ageGroup) )) + geom_histogram(binwidth = .001) +
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50)) 
png('allIQEffectSize.png',         # File name
    width=1024, height=768)
print(outPlot)
dev.off()

ggplot(results, aes(x=age, color = factor(ageGroup), fill = factor(ageGroup) ) ) + geom_histogram(binwidth = .001)+
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50))

ggplot(results, aes(x=sex, color = factor(ageGroup), fill = factor(ageGroup) ) ) + geom_histogram(binwidth = .001)+
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50))

ggplot(results, aes(x=IQ, color = factor(ageGroup), fill = factor(ageGroup) ) ) + geom_histogram(binwidth = .001)+
  coord_cartesian(xlim=c(0,.6), ylim = c(0,50))

write.csv(results, 'results.csv')








# 
# fpfn = logspace(-3.5, 3.5, 100) #c(.001, .01, .1, 1, 10, 100, 1000)
# tptn = logspace(-3.5, 3.5, 100) #c(.001, .01, .1, 1, 10, 100, 1000)
# baseVals = c(.0148, .187)
# for(agei in 1:3){
#   temp = resultsDif[resultsDif$ageGroup == agei,]
#   tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
#   candidates = which(temp$Diag > .8 & tempRes$Diag > thresh)
#   if(agei==1){
#     target = candidates[6]
#   }else{
#     target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
#   }
#   targName = temp$dependentVariable[target]
#   limVals = quantile(comboDat[,targName], c(.025,.975))
#   # limVals = c(limVals[1] - (limVals[2]-limVals[1])*.5,
#   #             limVals[2] + (limVals[2]-limVals[1])*.5)
#   
#   curDat = comboDat %>% filter(ageGroup==agei)
#   groupID = rep(1,length(curDat$group))
#   groupID[curDat$group != 'CON'] = 2
#   #ALL AD AND ASD SUBS CODED AS 2
#   biomark = curDat[,targName]
#   criteria = rev(seq(limVals[1], limVals[2], (limVals[2] - limVals[1])/100))
#   criteria = criteria[1:100]
#   #is the AD group higher or lower than the CON group? 
#   #default is autistic < control
#   if(agei==1){
#     biomark = biomark * -1
#   }
#   ROC = getROC(biomark, groupID, criteria)
#   
#   for(bb in 1:2){
#     baseRate = baseVals[bb]
#     ADGroup = sample(biomark[groupID==2], round(baseRate*1000000), replace = T )
#     CONGroup = sample(biomark[groupID==1], round((1-baseRate)*1000000), replace = T)
#     simBiomark = c(ADGroup, CONGroup)
#     simID = rep(1,length(simBiomark))
#     simID[1:length(ADGroup)] = 2
#     simROC = getROC(simBiomark, simID, criteria)
#     simROC$TPR[simROC$TPR==0] = 1/1000000
#     simROC$FPR[simROC$FPR==0] = 1/1000000
#     results = expand.grid(tptn=tptn, fpfn=fpfn)
#     results$TPR = NaN
#     results$FPR = NaN
#     results$acc = NaN
#     for(ii in 1:length(results$TPR)){
#       
#       ratio = simROC$TPR/simROC$FPR
#       costVal = ((1-baseRate)/baseRate) * ((results$tptn[ii]-1)/(1-1*results$fpfn[ii]))
#       
#       critVal = which(abs(ratio-costVal)==min(abs(ratio-costVal)))
#       if(length(critVal)>1){
#         critVal = critVal[1]
#       }
#       if(length(critVal)==1){
#         if(is.finite(critVal)){
#           results$TPR[ii] = simROC$TPR[critVal]
#           results$FPR[ii] = simROC$FPR[critVal]
#           results$acc[ii] = simROC$acc[critVal]
#           
#         }
#       }
#     }
#     breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
#     outPlot <- ggplot(results, aes(x=tptn, y = fpfn, fill = TPR)) + 
#       geom_tile() + 
#       scale_y_log10(breaks = breakPoints) +
#       scale_x_log10(breaks = breakPoints) + 
#       scale_fill_gradient(low = '#CE2029', high = '#4D7902', na.value = 'black',
#                           limits = c(0, 1)) + 
#       theme_classic() +
#       ylab('FN<FP      FN/FP cost     FN>FP') + 
#       xlab('TN<TP      TN/TP value    TN>TP') +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) + 
#       ggtitle(paste('TPR as a function of costs and benefits ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10))) 
#     
#     locs = seq(-3.5,3.5,1)  
#     for(vv in 1:8){
#       outPlot = outPlot + 
#         geom_vline(xintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2) + 
#         geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     png(paste( 'TPR_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     
#     outPlot <- ggplot(results, aes(x=tptn, y = fpfn, fill = FPR)) + 
#       geom_tile() + 
#       scale_y_log10(breaks = breakPoints) +
#       scale_x_log10(breaks = breakPoints) + 
#       scale_fill_gradient(low = '#4D7902', high = '#CE2029', na.value = 'black',
#                           limits = c(0, 1)) + 
#       theme_classic() +
#       ylab('FN<FP      FN/FP cost     FN>FP') + 
#       xlab('TN<TP      TN/TP value    TN>TP') +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) + 
#       ggtitle(paste('FPR as a function of costs and benefits ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10))) 
#     
#     locs = seq(-3.5,3.5,1)  
#     for(vv in 1:8){
#       outPlot = outPlot + 
#         geom_vline(xintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2) + 
#         geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     
#     png(paste( 'FPR_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     results$accRound = results$acc
#     results$accRound[results$acc>=.60] = .59
#     results$accRound[results$acc<=.50] = .50
#     
#     breakPoints =  c(.001, .01, .1, 1, 10, 100, 1000)
#     outPlot <- ggplot(results, aes(x=tptn, y = fpfn, fill = accRound)) + 
#       geom_tile() + 
#       scale_y_log10(breaks = breakPoints) +
#       scale_x_log10(breaks = breakPoints) + 
#       scale_fill_gradient(low = '#FFFFFF', high = '#3A317D', na.value = 'black',
#                           limits = c(.5, .60)) + 
#       theme_classic() +
#       ylab('FN<FP      FN/FP cost     FN>FP') + 
#       xlab('TN<TP      TN/TP value    TN>TP') +
#       theme(axis.line = element_line(color = 'black', size = 3),
#             axis.ticks = element_line(colour = "black", size = 2),
#             axis.ticks.length=unit(-.25, "cm"),
#             text = element_text(size = 20)) + 
#       ggtitle(paste('accuracy as a function of costs and benefits ', as.character(agei), sep = ''))+
#       guides(color = guide_legend(override.aes = list(size=10))) 
#     
#     locs = seq(-3.5,3.5,1)  
#     for(vv in 1:8){
#       outPlot = outPlot + 
#         geom_vline(xintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2) + 
#         geom_hline(yintercept = 10^locs[vv], color = 'grey', size = 2, linetype =2)
#       
#     }
#     png(paste( 'acc_costmap_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     outPlot <- ggplot(results, aes(x = accRound)) + 
#       geom_histogram()
#     png(paste( 'acc_hist_age', agei, 'base', bb, '.png', sep = ''),         # File name
#         width=1024, height=768)
#     print(outPlot)
#     dev.off()
#     
#     
#     ggplot(results, aes(x = FPR, y = TPR)) + geom_jitter()
#     
#     
#   }
#   
#   
#   
# }