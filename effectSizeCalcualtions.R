
library(lsr)

library(tidyverse)
# library(ggpubr)
# library(rstatix)
# library(broom)
library(car)
library(circular)


test <- read.csv("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\autismBiomarkersAllData2.csv")
test$nbChanOrig[test$nbChanOrig==999] = 124

dat <- filter(test, group %in% c('CON', 'AD', 'ASD'), eyes == 'open', nbChanFinal/nbChanOrig >= .5, nbChanOrig>20)

setwd("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\figures")

#### some descriptives about the overall data set ####
png("ageAll.png",         # File name
    width=1024, height=768)# Color
ggplot(dat, aes(x = age, color = group, fill = group)) + geom_histogram()
rect(1, 5, 3, 7, col="white")
dev.off()

dat <- filter(test, group %in% c('CON', 'AD', 'ASD'), eyes == 'open', age<250, nbChanFinal/nbChanOrig >= .5, nbChanOrig>20)

png("ageYoung.png",         # File name
    width=1024, height=768)# Color
ggplot(dat, aes(x = age, color = group, fill = group)) + geom_histogram()
rect(1, 5, 3, 7, col="white")
dev.off()

png("ageDatSet.png",         # File name
    width=1024, height=768)
ggplot(dat, aes(x = age, color = dataSet, fill = dataSet)) + geom_histogram()
rect(1, 5, 3, 7, col="white")
dev.off()

png("sex.png",         # File name
    width=1024, height=768)
ggplot(dat, aes(x = sex, color = group, fill = group)) + geom_bar()
rect(1, 5, 3, 7, col="white")
dev.off()

ageGroups = quantile(filter(dat, group %in% c('AD', 'ASD'))$age , c(0,.33333,.666666,1.0))
ageGroups[1]= 0
ageGroups[4]=999

#### randomization code written on 29/11/2022 and run one time
# 
# dat <- dat %>% arrange(age)
# allIdx = c()
# for(ii in 2:length(ageGroups)){
#   temp = which(dat$age>ageGroups[ii-1] & dat$age<=ageGroups[ii] )
#   tempS = dat$sex[temp]
#   tempG = dat$group[temp]
# 
#   keepers = c()
#   #do removal separately for male and female
#   temp_f = temp[tempS=='F']
#   temp_m = temp[tempS=='M']
#   #bring diagnosis info along
#   temp_fG = tempG[tempS=='F']
#   temp_mG = tempG[tempS=='M']
#   #do removal for each diagnosis
#   for(di in c('AD', 'ASD', 'CON')){
#     curM = temp_m[temp_mG == di]
#     curF = temp_f[temp_fG == di]
#     keepers = c(keepers, sample(curM, round(length(curM)*.5)))
#     keepers = c(keepers, sample(curF, round(length(curF)*.5)))
#   }
#   allIdx = c(allIdx, is.element(temp, keepers))
# 
# }
# heldOut = dat[allIdx == F, ]
# dat = dat[allIdx,]
# write.csv(dat, "C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
# write.csv(heldOut, "C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\testSet.csv")

#### end randomization block ####

dat = read.csv("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
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
                     'ageGroup' = rep(0, 4*length(varNames)))

setwd("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\figures")
for(ii in 2:length(ageGroups)){
  #down select for age
  curDat = dat[dat$age>ageGroups[ii-1] & dat$age<ageGroups[ii], ]
  #down select for group ID
  # curDat = curDat[curDat$group == 'CON' | curDat$group == 'AD' |curDat$group == 'ASD',]
  groupID = rep(1,length(curDat$group))
  #running the analysis with AD and ASD as a single group
  groupID[curDat$group=='AD'] = 3
  groupID[curDat$group=='ASD'] = 2
  
  

  for(tt in 1:length(varNames)){
    if(!grepl('PACmi', varNames[tt] )) { #skipping the non z-scored raw PAC values
    print(tt)
    temp = curDat
    #### outlier removal ####
    #after beginning analysis, it was discovered that some variables contained
    #extreme outliers such that one or two variables were far above the mean
    #and all other variables were below the mean.
    #in these cases, simply eliminating outliers based on a standard z-score
    #would have eliminated all variables since the mean ended up in a position
    #far from any data points. To get around this, data were z scored using only
    #the middle 80% of the data. Then, based on this z-scoring, outliers were removed
    if(!grepl('phase', varNames[tt])){
      keep = !is.na(curDat[[varNames[tt]]])
      temp = temp[keep,]
      dv = temp[[varNames[tt]]]
      n = length(dv)
      lims = quantile(dv, c(.1, .90))
      dv_mean = mean(dv[dv>lims[1] & dv<lims[2]])
      dv_sd = sd(dv[dv>lims[1] & dv<lims[2]])
      dv_z = (dv - dv_mean) / dv_sd
      nonOut = which(abs(dv_z)<5)
      dv_mean = mean(dv[nonOut])
      dv_sd = sd(dv[nonOut])
      dv_z = (dv - dv_mean) / dv_sd
      
      dv = dv[abs(dv_z)<5]
      age = temp$age[abs(dv_z)<5]
      sex = temp$sex[abs(dv_z)<5]
      IQ = temp$IQ[abs(dv_z)<5]
      IQ[is.na(IQ)] = 100
      gID = groupID[keep]
      ID = gID[abs(dv_z)<5]
      #standardizing the dependent variable prior to model calculation
      #it was discovered after data analysis began that many values, particularly
      #PAC values could not easily be fit because values were so small that floating
      #point precision had a problem. To get around this, all variables were z-scored
      #prior to model fitting. Mean and standard deviation are stored so that beta
      #weights can be converted back later if necessary
      dv_mean = mean(dv)
      dv_sd = sd(dv)
      dv = (dv - dv_mean) / dv_sd
    }else{
      dv_sd = 1
      dv_mean = 0
      keep = !is.na(curDat[[varNames[tt]]])
      temp = temp[keep,]
      dv = temp[[varNames[tt]]]
      age = temp$age
      sex = temp$sex
      IQ = temp$IQ
      IQ[is.na(IQ)] = 100
      ID = groupID[keep]
    }
    
    
    
    modDat = data.frame('dv' = dv, 'age' = age, 'sex' = factor(sex), 'IQ' = IQ, 'diag' = factor(ID))
    
    
    
   


    # outPlot <- ggplot(modDat, aes(x = age, y = dv*dv_sd + dv_mean, color = diag, shape = sex, size = IQ)) +
    #   geom_jitter()+
    #   scale_color_manual(values=c("#4E554E", "#3BACDD", "#E1C271"),
    #                        labels=c('CON', 'ASD', 'AD')) +
    #   theme_classic() +
    #   theme(axis.line = element_line(color = 'black', size = 3),
    #         axis.ticks = element_line(colour = "black", size = 2),
    #         axis.ticks.length=unit(-.25, "cm"),
    #         text = element_text(size = 20)) +
    #   ylab(varNames[tt]) +
    #   ggtitle(paste(varNames[tt], '; age group: ', as.character(ii-1), sep = ''))+
    #   guides(color = guide_legend(override.aes = list(size=10)),
    #          shape = guide_legend(override.aes = list(size=10)))
    # # rect(1, 5, 3, 7, col="white")
    # png(paste(varNames[tt], '_', ii-1, '.png', sep = ''),         # File name
    #     width=1024, height=768)
    # print(outPlot)
    # dev.off()

    if(!grepl('phase', varNames[tt])){
    
   
    curLM = lm(dv ~ age + sex + diag  +IQ , data = modDat)
    aovTab = Anova(curLM, type = 3)
    
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
    ri = which(names(curLM$coefficients) == 'diag2')
    results$diag_b1[tt+ai] = curLM$coefficients[[ri]]
    ri = which(names(curLM$coefficients) == 'diag3')
    results$diag_b2[tt+ai] = curLM$coefficients[[ri]]
    
    }else { #phase needs to be treated differently because it's circular
    circdv = circular(dv, unit = 'degrees')
    dv_mean =  mean(circdv)[[1]]
    dv_sd = sd(circdv)[[1]]
    test = cbind(cos(dv*(pi/180)), sin(dv*(pi/180)))
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
    }
    
    #shared outputs for phase and other
    results$n[tt+ai] = n
    results$out[tt+ai] = n-length(dv)
    results$ageGroup[tt+ai] = ii-1
    results$dvMean = dv_mean
    results$dvSD = dv_sd
  }
  
 
  }
}

resultsTRAIN = results
#### now get the same set of results but from the test set! 
dat = read.csv("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\testSet.csv")
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
                      'ageGroup' = rep(0, 4*length(varNames)))


setwd("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\figures")
for(ii in 2:length(ageGroups)){
  #down select for age
  curDat = dat[dat$age>ageGroups[ii-1] & dat$age<ageGroups[ii], ]
  #down select for group ID
  # curDat = curDat[curDat$group == 'CON' | curDat$group == 'AD' |curDat$group == 'ASD',]
  groupID = rep(1,length(curDat$group))
  #running the analysis with AD and ASD as a single group
  groupID[curDat$group=='AD'] = 3
  groupID[curDat$group=='ASD'] = 2
  
  
  
  for(tt in 1:length(varNames)){
    if(!grepl('PACmi', varNames[tt] )) { #skipping the non z-scored raw PAC values
      print(tt)
      temp = curDat
      #### outlier removal ####
      #after beginning analysis, it was discovered that some variables contained
      #extreme outliers such that one or two variables were far above the mean
      #and all other variables were below the mean.
      #in these cases, simply eliminating outliers based on a standard z-score
      #would have eliminated all variables since the mean ended up in a position
      #far from any data points. To get around this, data were z scored using only
      #the middle 80% of the data. Then, based on this z-scoring, outliers were removed
      if(!grepl('phase', varNames[tt])){
        keep = !is.na(curDat[[varNames[tt]]])
        temp = temp[keep,]
        dv = temp[[varNames[tt]]]
        n = length(dv)
        lims = quantile(dv, c(.1, .90))
        dv_mean = mean(dv[dv>lims[1] & dv<lims[2]])
        dv_sd = sd(dv[dv>lims[1] & dv<lims[2]])
        dv_z = (dv - dv_mean) / dv_sd
        nonOut = which(abs(dv_z)<5)
        dv_mean = mean(dv[nonOut])
        dv_sd = sd(dv[nonOut])
        dv_z = (dv - dv_mean) / dv_sd
        
        dv = dv[abs(dv_z)<5]
        age = temp$age[abs(dv_z)<5]
        sex = temp$sex[abs(dv_z)<5]
        IQ = temp$IQ[abs(dv_z)<5]
        IQ[is.na(IQ)] = 100
        gID = groupID[keep]
        ID = gID[abs(dv_z)<5]
        #standardizing the dependent variable prior to model calculation
        #it was discovered after data analysis began that many values, particularly
        #PAC values could not easily be fit because values were so small that floating
        #point precision had a problem. To get around this, all variables were z-scored
        #prior to model fitting. Mean and standard deviation are stored so that beta
        #weights can be converted back later if necessary
        dv_mean = mean(dv)
        dv_sd = sd(dv)
        dv = (dv - dv_mean) / dv_sd
      }else{
        dv_sd = 1
        dv_mean = 0
        keep = !is.na(curDat[[varNames[tt]]])
        temp = temp[keep,]
        dv = temp[[varNames[tt]]]
        age = temp$age
        sex = temp$sex
        IQ = temp$IQ
        IQ[is.na(IQ)] = 100
        ID = groupID[keep]
      }
      
      
      
      modDat = data.frame('dv' = dv, 'age' = age, 'sex' = factor(sex), 'IQ' = IQ, 'diag' = factor(ID))
      
      
      
      
      
     
      
      if(!grepl('phase', varNames[tt])){
        
        
        curLM = lm(dv ~ age + sex + diag  +IQ , data = modDat)
        aovTab = Anova(curLM, type = 3)
        
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
        ri = which(names(curLM$coefficients) == 'diag2')
        results$diag_b1[tt+ai] = curLM$coefficients[[ri]]
        ri = which(names(curLM$coefficients) == 'diag3')
        results$diag_b2[tt+ai] = curLM$coefficients[[ri]]
        
      }else { #phase needs to be treated differently because it's circular
        circdv = circular(dv, unit = 'degrees')
        dv_mean =  mean(circdv)[[1]]
        dv_sd = sd(circdv)[[1]]
        test = cbind(cos(dv*(pi/180)), sin(dv*(pi/180)))
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
      }
      
      #shared outputs for phase and other
      results$n[tt+ai] = n
      results$out[tt+ai] = n-length(dv)
      results$ageGroup[tt+ai] = ii-1
      results$dvMean = dv_mean
      results$dvSD = dv_sd
    }
    
    
  }
}

resultsTEST = results
#### result difference scores #### 
resultsDif = data.frame('dependentVariable' = resultsTEST$dependentVariable, 
                       'age' = abs(resultsTEST$age - resultsTRAIN$age) / (resultsTEST$age + resultsTRAIN$age),
                       'sex' = abs(resultsTEST$sex - resultsTRAIN$sex) / (resultsTEST$sex + resultsTRAIN$sex),
                       'IQ' = abs(resultsTEST$IQ - resultsTRAIN$IQ) / (resultsTEST$IQ + resultsTRAIN$IQ),
                       'Diag' = abs(resultsTEST$Diag - resultsTRAIN$Diag) / (resultsTEST$Diag + resultsTRAIN$Diag),
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
setwd("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\figures\\summaryFigs")

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


#### plotting the effect size values ####
thresh = .00
predNames = c('age', 'sex', 'IQ', 'Diag')


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
            axis.line = element_line(colour = 'black', linewidth = 3),
            axis.ticks = element_line(colour = 'black', linewidth = 3)) + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, .65), breaks = seq(0,.6,.1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0,100,20)) + 
      geom_vline(xintercept = .05, linetype = 'dashed', size = 5)
      
    png(paste('effect_', predNames[vari], '_', as.character(agei), '.png',sep=''),         # File name
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
outPlot <-ggplot(allEffects, aes(x=agreement, y=effectSize, color = predictor)) + 
                 geom_jitter(size = 5, alpha = .5)+
                 theme_classic(base_size = 30)+
                 theme(axis.ticks.length=unit(-0.5, "cm"),
                 axis.line = element_line(colour = 'black', size = 3),
                 axis.ticks = element_line(colour = 'black', size = 3),
                 legend.position = c(0.67, 0.8)) + 
                 guides(color = guide_legend(override.aes = list(size=10),
                                             nrow=4)) + 
                 ylab('effect size (\U1D702\U00B2\U209A)') + 
                 scale_x_continuous(expand = c(0, 0), limits = c(0, 1.01), breaks = seq(0,1,.25)) +
                 scale_y_continuous(expand = c(0, 0), limits = c(0, .65), breaks = seq(0,.65,.1))
png('effect_agree_relationship.png',         # File name
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
      scale_y_continuous(expand = c(0, 0))
    
    png(paste('stability_', predNames[vari], '_', '.png',sep=''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()       
  # }  
}



#### choosing the best plot for each 
dat2 = read.csv("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
dat2 <- dat2 %>% select(-X)
comboDat = rbind(dat2, dat)
predNames2 = predNames
predNames2[4] = 'group'
setwd("C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\figures\\summaryFigs")

## PLOTS FOR AGE! 
vari = 1
  for(agei in 1:3){
    temp = resultsDif[resultsDif$ageGroup == agei,]
    tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
    candidates = which(temp[predNames[vari]] < .2 & tempRes[predNames[vari]] > thresh)
    target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
  if(length(target)>0){
    limVals = quantile(comboDat[,varNames[target]], c(.1,.9))
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.2,
                limVals[2] + (limVals[2]-limVals[1])*.2)
  outPlot <- ggplot(comboDat, aes_string(x = 'age', y = varNames[target], color = 'group', shape = 'sex', size = 'IQ')) +
    geom_jitter(alpha = .75)+
    scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                         labels=c(  'AD', 'ASD', 'CON')) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) +
    ylim(limVals) +
    ylab(varNames[target]) +
    ggtitle(paste(varNames[target], '; age group: ', as.character(agei), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10)),
           shape = guide_legend(override.aes = list(size=10))) + 
    geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
    geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
  # rect(1, 5, 3, 7, col="white")
  png(paste( 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
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
  candidates = which(temp[predNames[vari]] < .2 & tempRes[predNames[vari]] > thresh)
  target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
  if(length(target)>0){
    limVals = quantile(comboDat[,varNames[target]], c(.1,.9))
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.2,
                limVals[2] + (limVals[2]-limVals[1])*.2)
    outPlot <- ggplot(comboDat, aes_string(x = 'age', y = varNames[target], color = 'sex', shape = 'group', size = 'IQ')) +
      geom_jitter(alpha = .75)+
      # scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
      #                    labels=c(  'AD', 'ASD', 'CON')) +
      theme_classic() +
      theme(axis.line = element_line(color = 'black', size = 3),
            axis.ticks = element_line(colour = "black", size = 2),
            axis.ticks.length=unit(-.25, "cm"),
            text = element_text(size = 20)) +
      ylim(limVals) +
      ylab(varNames[target]) +
      ggtitle(paste(varNames[target], '; age group: ', as.character(agei), sep = ''))+
      guides(color = guide_legend(override.aes = list(size=10)),
             shape = guide_legend(override.aes = list(size=10))) + 
      geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
      geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
    # rect(1, 5, 3, 7, col="white")
    png(paste( 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }
}

## PLOTS FOR IQ! 
vari = 3
for(agei in 1:3){
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  candidates = which(temp[predNames[vari]] < .2 & tempRes[predNames[vari]] > thresh)
  target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
  if(length(target)>0){
    limVals = quantile(comboDat[,varNames[target]], c(.1,.9))
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.2,
                limVals[2] + (limVals[2]-limVals[1])*.2)
    outPlot <- ggplot(comboDat, aes_string(x = 'IQ', y = varNames[target], color = 'age', shape = 'group')) +
      geom_jitter(alpha = .75, size = 7)+
      # scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
      #                    labels=c(  'AD', 'ASD', 'CON')) +
      theme_classic() +
      theme(axis.line = element_line(color = 'black', size = 3),
            axis.ticks = element_line(colour = "black", size = 2),
            axis.ticks.length=unit(-.25, "cm"),
            text = element_text(size = 20)) +
      ylim(limVals) +
      ylab(varNames[target]) +
      ggtitle(paste(varNames[target], '; age group: ', as.character(agei), sep = ''))+
      guides(color = guide_legend(override.aes = list(size=10)),
             shape = guide_legend(override.aes = list(size=10)))  
      # geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
      # geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
    # rect(1, 5, 3, 7, col="white")
    png(paste( 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }
}

## PLOTS FOR Diag! 
vari = 4
for(agei in 1:3){
  temp = resultsDif[resultsDif$ageGroup == agei,]
  tempRes = resultsTRAIN[resultsDif$ageGroup == agei,]
  candidates = which(temp[predNames[vari]] < .2 & tempRes[predNames[vari]] > thresh)
  target = candidates[which(tempRes[candidates,predNames[vari]]==max(tempRes[candidates,predNames[vari]]))]
  if(length(target)>0){
    limVals = quantile(comboDat[,varNames[target]], c(.1,.9))
    limVals = c(limVals[1] - (limVals[2]-limVals[1])*.25,
                limVals[2] + (limVals[2]-limVals[1])*.25)
    outPlot <- ggplot(comboDat, aes_string(x = 'age', y = varNames[target], color = 'group', shape = 'sex', size = 'IQ')) +
      geom_jitter(alpha = .75)+
      scale_color_manual(values=c( "#E1C271", "#3BACDD","#4E554E"),
                         labels=c(  'AD', 'ASD', 'CON')) +
      theme_classic() +
      theme(axis.line = element_line(color = 'black', size = 3),
            axis.ticks = element_line(colour = "black", size = 2),
            axis.ticks.length=unit(-.25, "cm"),
            text = element_text(size = 20)) +
      ylim(limVals) +
      scale_size(range = c(3,10))+
      ylab(varNames[target]) +
      ggtitle(paste(varNames[target], '; age group: ', as.character(agei), sep = ''))+
      guides(color = guide_legend(override.aes = list(size=10)),
             shape = guide_legend(override.aes = list(size=10))) + 
      geom_vline(xintercept = ageGroups[2], linetype = 'dashed', linewidth = 2, alpha = .75)+ 
      geom_vline(xintercept = ageGroups[3], linetype = 'dashed', linewidth = 2, alpha = .75)
    # rect(1, 5, 3, 7, col="white")
    png(paste( 'bestPredicted_',predNames[vari],'_', agei, '.png', sep = ''),         # File name
        width=1024, height=768)
    print(outPlot)
    dev.off()
  }
}










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

