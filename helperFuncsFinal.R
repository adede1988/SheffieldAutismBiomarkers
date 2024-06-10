
#helper functions 


getROC_prob <- function(biomark, groupID, criteria) {
  ROC = data.frame('crit' = criteria, 
                   'hits' = criteria,
                   'misses'=criteria,
                   'CRs' = criteria, 
                   'FAs' = criteria,
                   'TP' = criteria, 
                   'TN' = criteria,
                   'FP' = criteria, 
                   'FN' = criteria, 
                   'accRaw' = criteria)
  
  AUC = 0
  for(ci in 1:length(criteria)){
    guess = rep(1,length(biomark))
    guess[biomark<criteria[ci]] = 2
    ROC$hits[ci] = sum(guess==2 & groupID==2)
    ROC$misses[ci] = sum(guess==1 & groupID==2)
    ROC$CRs[ci] = sum(guess==1 & groupID==1)
    ROC$FAs[ci] = sum(guess==2 & groupID==1)
    ROC$TP[ci] = sum(guess==2 & groupID==2) / sum(groupID==2)
    ROC$TN[ci] = sum(guess==1 & groupID==1) / sum(groupID==1)
    ROC$FP[ci] = sum(guess==2 & groupID==1) / sum(groupID==1)
    ROC$FN[ci] = sum(guess==1 & groupID==2) / sum(groupID==2)
    ROC$accRaw[ci] = (sum(guess==2 & groupID==2) + sum(guess==1 & groupID==1)) / length(groupID)
    if(ci>1){
      baseDist = ROC$FP[ci] - ROC$FP[ci-1]
      #top triangle area:
      tria = ((ROC$TP[ci] - ROC$TP[ci-1]) * baseDist) / 2
      #base rectangle area:
      reca = ROC$TP[ci-1] * baseDist
      AUC = AUC + tria + reca
    } else {
      baseDist = ROC$FP[ci]
      #top triangle area:
      tria = ((ROC$TP[ci]) * baseDist) / 2
      AUC = AUC + tria
    }
    
  }
  baseDist = 1 - ROC$FP[ci]
  #top triangle area: 
  tria = ((1 - ROC$TP[ci]) * baseDist) / 2
  #base rectangle area: 
  reca = ROC$TP[ci] * baseDist
  AUC = AUC + tria + reca
  ROC$AUC = AUC
  
  
  #sensitivity, true positive rate
  ROC$TPR = ROC$hits/(ROC$hits + ROC$misses)
  #fallout, false positive rate
  ROC$FPR = ROC$FAs/(ROC$FAs + ROC$CRs)
  #specificity, true negative rate
  ROC$TNR = ROC$CRs/(ROC$CRs + ROC$FAs)
  #accuracy
  ROC$acc = (ROC$TP + ROC$TN) / (ROC$TP + ROC$TN + ROC$FP + ROC$FN)
  #what's the optimal criterion? 
  ci = which(ROC$acc == max(ROC$acc))
  if(length(ci)>1){
    ci = ci[1]
  }
  
  baseRate = sum(groupID==2) / length(groupID)
  #bayesian posterior: 
  post = ROC$TP[ci]*baseRate / (ROC$TP[ci]*baseRate + ROC$FP[ci]*(1-baseRate))
  ROC$base = baseRate
  ROC$post = post
  
  #false positive index
  guess = rep(1,length(biomark))
  guess[biomark>criteria[ci]] = 2
  ROC$FPi = sum(guess==2 & groupID==1) / sum(guess==2 & groupID==2)
  #precision, positive predictive value
  ROC$PPV = ROC$hits[ci] / (ROC$hits[ci] + ROC$FAs[ci])
  #negative predictive value
  ROC$NPV = ROC$CRs[ci] / (ROC$CRs[ci] + ROC$misses[ci])
  
  
  
  return(ROC)
}

getResults <- function(dat, varNames, ageGroups){
  #results for this bootstrap sample only 
  results = data.frame('dependentVariable' = rep('A', 4*length(varNames)),
                       'age' = rep(0, 4*length(varNames)),
                       'sex' = rep(0, 4*length(varNames)),
                       'IQ' = rep(0, 4*length(varNames)),
                       'Diag' = rep(0, 4*length(varNames)),
                       'diag_gvif'= rep(0, 4*length(varNames)),
                       'DiagBinom' = rep(0, 4*length(varNames)),
                       'Diag_p' = rep(0, 4*length(varNames)),
                       'age_b' = rep(0, 4*length(varNames)),
                       'sex_b' = rep(0, 4*length(varNames)),
                       'IQ_b' = rep(0, 4*length(varNames)),
                       'diag_b1' = rep(0, 4*length(varNames)),
                       'diag_b2' = rep(0, 4*length(varNames)),
                       'n' = rep(0, 4*length(varNames)),
                       'out' = rep(0, 4*length(varNames)),
                       'outlierFirst' = rep(0, 4*length(varNames)),
                       'ageGroup' = rep(0, 4*length(varNames)),
                       'logBase' = rep(0, 4*length(varNames)),
                       'logAug' = rep(0, 4*length(varNames)),
                       'nonLin' = rep(0, 4*length(varNames)),
                       'dvMean' = rep(0, 4*length(varNames)),
                       'dvSD' = rep(0, 4*length(varNames)),
                       'RsqBase' = rep(NA, 4*length(varNames)),
                       'RsqALL' = rep(NA, 4*length(varNames)))
  
  
  
  for(ii in 2:length(ageGroups)){
    #down select for age
    curDat = dat[dat$age>ageGroups[ii-1] & dat$age<=ageGroups[ii], ]
    #down select for group ID
    # curDat = curDat[curDat$group == 'CON' | curDat$group == 'AD' |curDat$group == 'ASD',]
    groupID = rep('CON',length(curDat$group))
    
    groupID[curDat$group=='AD'] = 'AD'
    groupID[curDat$group=='ASD'] = 'ASD'
    # 
    outliers = matrix(0,length(varNames), 3)
    
    for(tt in 1:length(varNames)){
      tryCatch({
      if(!grepl('PACmi', varNames[tt] )) { #skipping the non z-scored raw PAC values
        # print(paste('ageGroup: ', ii-1, 'var: ', tt))
        temp = curDat
        #### outlier removal ####
        #after beginning analysis, it was discovered that some variables contained
        #extreme outliers such that one or two observations were far above the mean
        #and all other variables were below the mean.
        #in these cases, simply eliminating outliers based on a standard z-score
        #would have eliminated all variables since the mean ended up in a position
        #far from any data points. To get around this, data were z scored using only
        #the middle 80% of the data. Then, based on this z-scoring, outliers were removed
        tmpOut <- getModDat(dat, varNames, tt, temp, ageGroups, ii, groupID)
        tmpDat <- tmpOut[[1]]
        dat <- tmpOut[[2]]
        modDat <- tmpDat
        
        modDat$IQ[is.nan(modDat$IQ)] = mean(modDat$IQ, na.rm = T)
        modDat$IQ[is.na(modDat$IQ)] = mean(modDat$IQ, na.rm = T)
        #optional plot for this variable
        # modelPlot(modDat, varNames[tt], ii)
        
        
        if(!grepl('phase', varNames[tt])){
          
          modDat$diag <- factor(modDat$diag, levels = c("CON", "ASD", "AD"))
          modDat$sex <- factor(modDat$sex, levels = c( "F","M"))  
          
          
          
          #fit basic non-interactive multiple linear regression
          curLM = lm(dv ~ age + sex + diag  +IQ , data = modDat)
          #fit model with diag*sex interaction
          curLMSex = lm(dv~age+sex+diag+diag*sex, data = modDat)
          #fit model with diag*age interaction
          curLMAge = lm(dv~age+sex+IQ+diag+diag*age, data = modDat)
          #fit model with diag*age*sex interaction
          curLMInt = lm(dv~age*sex*diag+IQ, data = modDat)
          
          #fit a base model with age*sex*diag interaction then predict residuals with diag
          curLMBase = lm(dv~age*sex*I(age^2)+IQ, data = modDat)
          curLMALL = lm(dv~age*sex*I(age^2)+IQ+diag, data = modDat)
          
          all_base_test = anova(curLMBase, curLMALL)
          SStot = sum((modDat$dv - mean(modDat$dv))^2)
          SSresBase = all_base_test$RSS[1]
          SSresALL = all_base_test$RSS[2]
          RsqBase = 1 - SSresBase/SStot
          RsqALL = 1 - SSresALL/SStot
          
          
          curLMresiduals = lm(curLMBase$residuals ~ modDat$diag)
          
         
          
          
          
          
          modDat$testDiag = modDat$diag
          modDat$testDiag[modDat$testDiag == "ASD"] = 'AD'
          modDat$testDiag = modDat$testDiag == "AD"
          
          #auxilary modeling to see if diagnosis can be modeled as a function 
          #of predictors and whether adding in the EEG variable helps!
          auxLM = glm(testDiag ~ age + sex + IQ, data = modDat, family = "binomial")
          auxLM_aug = glm(testDiag ~ age + sex + dv + IQ, data = modDat, family = "binomial")
          test = lrtest(auxLM, auxLM_aug)
          groupIDROC = modDat$testDiag
          groupIDROC[groupIDROC == T] = 2
          groupIDROC[groupIDROC ==0] = 1
          bioMark = auxLM$fitted.values
          bioMark2 = auxLM_aug$fitted.values
          ROC = getROC(bioMark, groupIDROC, seq(.99,.01,-.01), F)
          ROC2 = getROC(bioMark2, groupIDROC, seq(.99, .01, -.01), F)
          baseAcc = ROC$AUC[1]
          augAcc = ROC2$AUC[1]
          
        
          
          
          #is a non-linear fit better for age? 
          curLM_nonLin = lm(dv~ age + I(age^2), data = modDat)
          curLM_ageAlone = lm(dv~age, data = modDat)
          
          nonLinTest = lrtest(curLM_ageAlone, curLM_nonLin)
          
          if( nonLinTest$`Pr(>Chisq)`[2] <.05) { #use the non linear if it's better
            curLM = lm(dv ~ age + sex + diag  +IQ + I(age^2), data = modDat)
            curLMSex = lm(dv ~ age + sex + diag  +IQ + I(age^2)+diag*sex, data = modDat)
            curLMAge = lm(dv ~ age + sex + diag  +IQ + I(age^2)+diag*age, data = modDat)
            curLMInt = lm(dv ~ age*sex*diag  +IQ + I(age^2), data = modDat)
            
          }
          
          aovTab = Anova(curLM, type = 3)
          
          ai = (ii-2)*length(varNames)
          
          results$RsqBase[tt+ai] = RsqBase
          results$RsqALL[tt+ai] = RsqALL
          
          results$logBase[tt+ai] = baseAcc
          results$logAug[tt+ai] = augAcc
          results$DiagBinom[tt+ai] = test$`Pr(>Chisq)`[2]
          
        
          
          results$nonLin[tt+ai] = nonLinTest$`Pr(>Chisq)`[2] #non lin better? 
          
          sink(file = "nul")
          
          # Call gvif and store its output
          vifTab <- gvif(curLM)  # Ensure to replace 'curLM' with your actual model variable
          
          # Restore console output
          sink()
          
          ri = which(grepl('diag', rownames(vifTab)))
          results$diag_gvif[tt+ai] = vifTab[ri,3]
          
          
          results[tt+ai,1] = varNames[tt]
          ei = which(row.names(aovTab) == 'Residuals')
          ri = which(row.names(aovTab) == 'diag')
          results$Diag[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          results$Diag_p[tt+ai] = aovTab$`Pr(>F)`[ri] 
          
          ri = which(row.names(aovTab) == 'sex')
          results$sex[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          
          ri = which(row.names(aovTab) == 'IQ')
          results$IQ[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          ri = which(row.names(aovTab) == 'age')
          results$age[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          
          
          ri = which(names(curLM$coefficients) == 'age')
          results$age_b[tt+ai] = curLM$coefficients[[ri]]
          ri = which(names(curLM$coefficients) == 'sexM')
          results$sex_b[tt+ai] = curLM$coefficients[[ri]]
          ri = which(names(curLM$coefficients) == 'IQ')
          results$IQ_b[tt+ai] = curLM$coefficients[[ri]]
          ri = which(names(curLM$coefficients) == 'diagASD')
          if(length(ri) == 0){
            results$diag_b1[tt+ai] = NA
          } else { 
            results$diag_b1[tt+ai] = curLM$coefficients[[ri]]
          }
          ri = which(names(curLM$coefficients) == 'diagAD')
          if(length(ri) == 0 ){
            results$diag_b2[tt+ai] = NA
          } else {
            results$diag_b2[tt+ai] = curLM$coefficients[[ri]]
          }
          
          
          aovTab = Anova(curLMAge, type = 3, singular.ok = T)
          ei = which(row.names(aovTab) == 'Residuals')
          ri = which(row.names(aovTab) == 'age:diag')
          results$Diag_Age[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          aovTab = Anova(curLMSex, type = 3, singular.ok = T)
          ei = which(row.names(aovTab) == 'Residuals')
          ri = which(row.names(aovTab) == 'sex:diag')
          results$Diag_Sex[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          aovTab = Anova(curLMInt, type = 3, singular.ok = T)
          ei = which(row.names(aovTab) == 'Residuals')
          ri = which(row.names(aovTab) == 'age:sex:diag')
          results$Diag_Age_Sex[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          
      
          
          #look at the residual fit
          aovTab = Anova(curLMresiduals, type = 3, singular.ok = T)
          ei = which(row.names(aovTab) == 'Residuals')
          ri = which(grepl('diag', row.names(aovTab)))
          results$Diag_resid[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          results$Diag_residR2[tt+ai] = summary(curLMresiduals)$r.squared
          
          results$n[tt+ai] = length(temp[,1])
          results$out[tt+ai] = length(temp[,1])-length(modDat$dv)
          results$ageGroup[tt+ai] = ii-1
          results$dvMean[tt+ai] = modDat$dv_mean[1]
          results$dvSD[tt+ai] = modDat$dv_sd[1]
          
        }else { #phase needs to be treated differently because it's circular
          # circdv = circular(modDat$dv, unit = 'degrees')
          # dv_mean =  mean(circdv)[[1]]
          # dv_sd = sd(circdv)[[1]]
          # test = cbind(cos(modDat$dv*(pi/180)), sin(modDat$dv*(pi/180)))
          # modDat$cos_sin = test
          # 
          # 
          # 
          # 
          # 
          # curLM = lm(cos_sin ~ age + sex + diag  + IQ, data = modDat)
          # #fit model with diag*sex interaction
          # curLMSex = lm(cos_sin~age+sex+diag+diag*sex, data = modDat)
          # #fit model with diag*age interaction
          # curLMAge = lm(cos_sin~age+sex+IQ+diag+diag*age, data = modDat)
          # #fit model with diag*age*sex interaction
          # curLMInt = lm(cos_sin~age*sex*diag+IQ, data = modDat)
          # 
          # 
          # #fit a base model with age*sex*diag interaction then predict residuals with diag
          # curLMBase = lm(cos_sin~age*sex*IQ, data = modDat)
          # 
          # 
          # 
          # 
          # curLMresiduals = lm(curLMBase$residuals ~ modDat$diag)
          # 
          # #is a non-linear fit better for age? 
          # curLM_nonLin = lm(cos_sin~ age + I(age^2), data = modDat)
          # curLM_ageAlone = lm(cos_sin~age, data = modDat)
          # 
          # nonLinTest = anova(curLM_ageAlone, curLM_nonLin)
          # 
          # if( nonLinTest$`Pr(>F)`[2] <.05) { #use the non linear if it's better
          #   curLM = lm(cos_sin ~ age + sex + diag  +IQ + I(age^2), data = modDat)
          #   curLMSex = lm(cos_sin~age+sex+diag+diag*sex + I(age^2), data = modDat)
          #   curLMAge = lm(cos_sin~age+sex+IQ+diag+diag*age + I(age^2), data = modDat)
          #   curLMInt = lm(cos_sin~age*sex*diag+IQ + I(age^2), data = modDat)
          # }
          # 
          # 
          # #logistic regression based modeling
          # modDat$testDiag = modDat$diag
          # modDat$testDiag[modDat$testDiag == "ASD"] = 'AD'
          # modDat$testDiag = modDat$testDiag == "AD"
          # 
          # auxLM = glm(testDiag ~ age + sex + IQ, data = modDat, family = "binomial")
          # auxLM_aug = glm(testDiag ~ age + sex + IQ + cos_sin, data = modDat, family = "binomial")
          # test = lrtest(auxLM, auxLM_aug)
          # groupIDROC = modDat$testDiag
          # groupIDROC[groupIDROC == T] = 2
          # groupIDROC[groupIDROC ==0] = 1
          # bioMark = auxLM$fitted.values
          # bioMark2 = auxLM_aug$fitted.values
          # ROC = getROC(bioMark, groupIDROC, seq(.99,.01,-.01), F)
          # ROC2 = getROC(bioMark2, groupIDROC, seq(.99, .01, -.01), F)
          # baseAcc = ROC$AUC[1]
          # augAcc = ROC2$AUC[1]
          # 
          # auxLM = glm(testDiag ~ age + sex + IQ, data = modDat, family = "binomial")
          # auxLM_aug = glm(testDiag ~ age + sex + IQ + cos_sin + cos_sin*age, data = modDat, family = "binomial")
          # testAGE = lrtest(auxLM, auxLM_aug)
          # bioMark = auxLM$fitted.values
          # bioMark2 = auxLM_aug$fitted.values
          # ROC = getROC(bioMark, groupIDROC, seq(.99,.01,-.01), F)
          # ROC2 = getROC(bioMark2, groupIDROC, seq(.99, .01, -.01), F)
          # baseAccAGE = ROC$AUC[1]
          # augAccAGE = ROC2$AUC[1]
          # 
          # #create function to extract anova table
          # outtests <- car:::print.Anova.mlm
          # body(outtests)[[16]] <- quote(invisible(tests))
          # body(outtests)[[15]] <- NULL
          # #get anova table
          # tab <- lapply(c("Pillai"), 
          #               function(i)  outtests(Anova(curLM, test.statistic=i)))
          # tab <- do.call(rbind, tab)
          # 
          # 
          # 
          # 
          # 
          # 
          # ai = (ii-2)*length(varNames)
          # 
          # results$logBase[tt+ai] = baseAcc
          # results$logAug[tt+ai] = augAcc
          # results$DiagBinom[tt+ai] = test$`Pr(>Chisq)`[2] #did EEG var improve logistic prediction? 
          # 
          # 
          # results$logBaseAGE[tt+ai] = baseAccAGE
          # results$logAugAGE[tt+ai] = augAccAGE
          # results$DiagBinomAGE[tt+ai] = testAGE$`Pr(>Chisq)`[2] #did EEG var improve logistic prediction? 
          # 
          # results$nonLin[tt+ai] = nonLinTest$`Pr(>F)`[2] #non lin better? 
          # ri = which(grepl('diag', row.names(tab)))
          # results[tt+ai,1] = varNames[tt]
          # 
          # #working out effect size for diag
          # if(length(ri)==0) { 
          #   results$Diag[tt+ai] = NA
          # } else {
          #   results$Diag[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #     (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          #   results$Diag_p[tt+ai] = tab$`Pr(>F)`[ri]
          #   
          # }
          # ri = which(grepl('sex', row.names(tab)))
          # results$sex[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # 
          # ri = which(grepl('IQ', row.names(tab)))
          # results$IQ[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # ri = which(grepl('age', row.names(tab)))
          # results$age[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # 
          # #extract interaction terms: 
          # tab <- lapply(c("Pillai"), 
          #               function(i)  outtests(Anova(curLMAge, test.statistic=i)))
          # tab <- do.call(rbind, tab)
          # ri = which(row.names(tab) == 'age:diag')
          # results$Diag_Age[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # 
          # tab <- lapply(c("Pillai"), 
          #               function(i)  outtests(Anova(curLMSex, test.statistic=i)))
          # tab <- do.call(rbind, tab)
          # ri = which(row.names(tab) == 'sex:diag')
          # results$Diag_Sex[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # 
          # tab <- lapply(c("Pillai"), 
          #               function(i)  outtests(Anova(curLMInt, test.statistic=i)))
          # tab <- do.call(rbind, tab)
          # ri = which(row.names(tab) == 'age:sex:diag')
          # results$Diag_Age_Sex[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # 
          # 
          # 
          # #look at the anova table for fitting the residuals: 
          # tab <- lapply(c("Pillai"), 
          #               function(i)  outtests(Anova(curLMresiduals, test.statistic=i)))
          # tab <- do.call(rbind, tab)
          # ri = which(grepl('diag', row.names(tab)))
          # results$Diag_resid[tt+ai] = (tab$`num Df`[ri] * tab$`approx F`[ri]) / 
          #   (tab$`num Df`[ri] * tab$`approx F`[ri] + tab$`den Df`[ri])
          # 
          # results$n[tt+ai] = length(temp[,1])
          # results$out[tt+ai] = length(temp[,1])-length(modDat$dv)
          # results$ageGroup[tt+ai] = ii-1
          # results$dvMean[tt+ai] = dv_mean
          # results$dvSD[tt+ai] = dv_sd
        }
        
        
      }
      
        }, error = function(e) {
          # If an error occurs, this code will run
          print(paste("Error encountered with input:", tt, "; Skipping to next iteration."))
          ai = (ii-2)*length(varNames)
          
          results[tt+ai,1] = varNames[tt]
          results[tt+ai,2:38] = NA
        })
    }
  }
  results <- results %>% filter(dependentVariable != 'A')
  
  results$type = 'A'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('pow', x, ignore.case = T))] = 'power'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('slope', x, ignore.case = T))] = 'slope'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('relalpha', x, ignore.case = T) | 
                       grepl('logalpha', x, ignore.case = T))] = 'alpha'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('sampEnt', x, ignore.case = T))] = 'sampEnt'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('ispc', x, ignore.case = T))] = 'ispc'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('PAC', x, ignore.case = T))] = 'PAC'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('phase', x, ignore.case = T))] = 'phase'
  
  
  results$binomImprove = (results$logAug - results$logBase) / (1-results$logBase)

  return(list(results, dat))
}



makeDemoTable <- function(dat){
  dataSumTable = data.frame('group' = rep(c('AD', 'ASD', 'CON')), 
                            'data set' = as.vector(apply(as.matrix(unique(dat$dataSet)), 1, function(x) rep(x, 3))),
                            'IQ' = rep('', 15),
                            'IQ metric' = c(rep('DAS GCA', 3), rep('MSEL', 3), rep('WTAR', 3), rep('DAS GCA', 3), rep('WTAR', 3)),
                            'age' = rep('', 15),
                            'n female' = rep('', 15),
                            'n total' = rep('', 15),
                            'orig channels' = rep('', 15),
                            'final channels' = rep('', 15),
                            'orig epochs' = rep('', 15),
                            'final epochs' = rep('', 15))
  
  
  
  for(ii in 1:15) { 
    cur = filter(dat, group == dataSumTable$group[ii], dataSet == dataSumTable$data.set[ii])
    if(length(cur$age) > 0){
      dataSumTable$IQ[ii] = paste(round(mean(cur$IQ, na.rm = T)), ' (', round(sd(cur$IQ, na.rm = T)), ')', sep='')
      dataSumTable$age[ii] = paste(round(mean(cur$age)), ' (', round(sd(cur$age)), ')', sep='')
      dataSumTable$n.female[ii] = sum(cur$sex == 'F')
      dataSumTable$n.total[ii] = length(cur$sex)
      dataSumTable$orig.channels[ii] = paste(round(mean(cur$nbChanOrig)), ' (', round(sd(cur$nbChanOrig)), ')', sep='')
      dataSumTable$final.channels[ii] = paste(round(mean(cur$nbChanFinal)), ' (', round(sd(cur$nbChanFinal)), ')', sep='')
      dataSumTable$orig.epochs[ii] = paste(round(mean(cur$nbTrialOrig)), ' (', round(sd(cur$nbTrialOrig)), ')', sep='')
      dataSumTable$final.epochs[ii] = paste(round(mean(cur$nbTrialFinal)), ' (', round(sd(cur$nbTrialFinal)), ')', sep='')
      
      #replace na IQ values
      dat$IQ[dat$group == dataSumTable$group[ii] & dat$dataSet == dataSumTable$data.set[ii] & is.na(dat$IQ)] = mean(cur$IQ, na.rm = T)
    }
  }
  
  write.csv(dataSumTable, "G:\\My Drive\\Milne\\clusterDifferences\\table1.csv")
  
  dataSumTable %>% 
    kbl(align = 'c') %>% 
    kable_classic(full_width = F, 
                  font_size = 20) %>%
    row_spec(1, align = 'c')%>%
    footnote(general = "DAS GCA = Differential Ability Scales General Conceptual Ability
                      MSEL = Mullen Scales of Early Learning
                      WTAR = Weschler Test of Adult Reading
                      biomarkCon = The Autism Biomarkers Consortium for Clinical Trials
                      biomarkDev = Biomarkers of Developmental Trajectories and Treatment in ASD
                      bpSZ = Bipolar & Schizophrenia Consortium for Parsing Intermediate Phenotypes 
                      femaleASD = Multimodal Developmental Neurogenetics of Females with ASD
                      socBrain = The Social Brain in Schizophrenia and Autism Spectrum Disorders
                      Numeric values indicate mean and (standard deviation).",
             general_title = "Table 1: ",
             footnote_as_chunk = T, title_format = c("italic", "underline")
    )  
}


myNMI <- function(set1, set2){
  
  # Get unique cluster IDs and remove -1 (noise or outliers)
  clustIDs <- unique(set1[set1 != 0])
  clustIDs2 <- unique(set2[set2 != 0])
  
  # Calculate the normalized mutual information
  nominSum <- 0
  denomh <- 0
  denoml <- 0
  n <- length(set1)
  
  # Iterate over each cluster in set1
  for (c in clustIDs) {
    vec1 <- as.numeric(set1 == c)
    nh <- sum(vec1)
    
    # Iterate over each cluster in set2
    for (cc in clustIDs2) {
      vec2 <- as.numeric(set2 == cc)
      nl <- sum(vec2)
      vec3 <- as.numeric(vec1 == 1 & vec2 == 1)
      nhl <- sum(vec3)
      
      # Correction for zero overlap
      if (nhl == 0) {
        nhl <- .000000001
      }
      
      nominSum <- nominSum + nhl * log10((n * nhl) / (nh * nl))
    }
    
    denomh <- denomh + nh * log10(nh / n)
  }
  
  # Calculate the denominator for clusters in set2
  for (cc in clustIDs2) {
    vec2 <- as.numeric(set2 == cc)
    nl <- sum(vec2)
    denoml <- denoml + nl * log10(nl / n)
  }
  
  # Calculate NMI
  nmi <- nominSum / sqrt(denomh * denoml)
  
  # Check for values outside [0,1] range
  if (nmi < 0) {
    # warning("Check inputs for cluster scheme of all -1")
    nmi <- 0
  } else if (nmi > 1) {
    # warning("Check inputs for cluster scheme of all -1")
    nmi <- 1
  }
  
  # Print the NMI result
  return(nmi)
}

getModDat <- function(dat, varNames, tt, temp, ageGroups, ii, groupID){
  if(!grepl('phase', varNames[tt])){
    ageIdx = (dat$age > ageGroups[ii-1] & dat$age <= ageGroups[ii])
    ageIdx = which(ageIdx)
    keep = !is.na(dat[ageIdx,varNames[tt]] )
    keepi = which(keep)
    temp = temp[keep,]
    dv = temp[[varNames[tt]]]
    if(is.character(dv)){
      dv <- as.numeric(dv)
    }
    n = length(dv)
    lims = quantile(dv, c(.1, .90))
    dv_mean = mean(dv[dv>lims[1] & dv<lims[2]])
    dv_sd = sd(dv[dv>lims[1] & dv<lims[2]])
    dv_z = (dv - dv_mean) / dv_sd
    nonOut = which(abs(dv_z)<5)
    dv_mean = mean(dv[nonOut])
    dv_sd = sd(dv[nonOut])
    dv_z = (dv - dv_mean) / dv_sd
    
    
    #what group has outliers? 
    # outliers[tt,1] <<- sum(groupID[abs(dv_z)>5]==1)
    # outliers[tt,2] <<- sum(groupID[abs(dv_z)>5]==2)
    # outliers[tt,3] <<- sum(groupID[abs(dv_z)>5]==3)
    
    dv = dv[abs(dv_z)<5]
    age = temp$age[abs(dv_z)<5]
    sex = temp$sex[abs(dv_z)<5]
    IQ = temp$IQ[abs(dv_z)<5]
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
    # if(length(keep) != length(dv_z)){
    # print(paste(length(keep), length(dv_z), tt))
    # }
    keepz = abs(dv_z)<5
    
  }else{
    ageIdx = (dat$age > ageGroups[ii-1] & dat$age <= ageGroups[ii])
    ageIdx = which(ageIdx)
    dv_sd = 1
    dv_mean = 0
    keep = !is.na(dat[ageIdx,varNames[tt]])
    temp = temp[keep,]
    dv = temp[[varNames[tt]]]
    age = temp$age
    sex = temp$sex
    IQ = temp$IQ
    ID = groupID[keep]
    keepz = rep(T, length(dv))
  }
  
  dat[ageIdx[!keep], varNames[tt]] = NA #replace original NA values
  dat[ageIdx[!keepz], varNames[tt]] = NA #replace outlier values
  
  
  
  modDat = data.frame('dv' = dv, 'age' = age, 'sex' = factor(sex), 'IQ' = IQ, 'diag' = factor(ID))
  modDat$dv_sd = dv_sd
  modDat$dv_mean = dv_mean
  return(list(modDat, dat))
  
}

modelPlot <- function(modDat, varName, ii){
  ggplot(modDat, aes(x = age, y = dv*dv_sd + dv_mean, color = diag, shape = sex, size = IQ)) +
    geom_jitter()+
    scale_color_manual(values=c("#4E554E", "#3BACDD", "#E1C271"),
                       labels=c('CON', 'ASD', 'AD')) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) +
    ylab(varName) +
    ggtitle(paste(varName, '; age group: ', as.character(ii-1), sep = ''))+
    guides(color = guide_legend(override.aes = list(size=10)),
           shape = guide_legend(override.aes = list(size=10)))
}


randomSplit <- function(dat, ageGroups, prop = .5){
  dat <- dat %>% arrange(age)
  allIdx = c()
  for(ii in 2:length(ageGroups)){
    temp = which(dat$age>ageGroups[ii-1] & dat$age<=ageGroups[ii] )
    tempS = dat$sex[temp]
    tempG = dat$group[temp]

    keepers = c()
    #do removal separately for male and female
    temp_f = temp[tempS=='F']
    temp_m = temp[tempS=='M']
    #bring diagnosis info along
    temp_fG = tempG[tempS=='F']
    temp_mG = tempG[tempS=='M']
    #do removal for each diagnosis
    for(di in c('AD', 'ASD', 'CON')){
      curM = temp_m[temp_mG == di]
      curF = temp_f[temp_fG == di]
      keepers = c(keepers, sample(curM, round(length(curM)*prop)))
      keepers = c(keepers, sample(curF, round(length(curF)*prop)))
    }
    allIdx = c(allIdx, is.element(temp, keepers))

  }
  heldOut = dat[allIdx == F, ]
  dat = dat[allIdx,]
  return(list(heldOut, dat))
  
}






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
  


getBalanced <- function(modDat){
  modDat$ID = seq(1,length(modDat$group), 1)
  #go through the model Data and create a new dataFrame that has matched samples
  outDat <- modDat[0,]
  
  ADgroup = modDat %>% filter(group == 'AD')
  ASDgroup = modDat %>% filter(group == 'ASD')
  CONgroup = modDat %>% filter(group == 'CON')
  oi = 1
  for(jj in 1:length(ADgroup$group)){
   
    curAD = ADgroup[jj,]
    sampTmp = CONgroup %>% filter(sex == curAD$sex, 
                                  age>curAD$age-5,
                                  age<curAD$age+5, 
                                  IQ>curAD$IQ-10,
                                  IQ<curAD$IQ+10)
    #check that we have at least one match option: 
    if(length(sampTmp$group)>0){
      outDat[oi,] = sample_n(sampTmp, 1)
      CONgroup <- CONgroup %>% filter(ID != outDat$ID[oi])
      oi = oi+1
      outDat[oi,] = curAD
      oi = oi+1
    }
    
    
  }
  
  for(jj in 1:length(ASDgroup$group)){
    curAD = ASDgroup[jj,]
    sampTmp = CONgroup %>% filter(sex == curAD$sex, 
                                  age>curAD$age-5,
                                  age<curAD$age+5, 
                                  IQ>curAD$IQ-10,
                                  IQ<curAD$IQ+10)
    #check that we have at least one match option: 
    if(length(sampTmp$group)>0){
      outDat[oi,] = sample_n(sampTmp, 1)
      CONgroup <- CONgroup %>% filter(ID != outDat$ID[oi])
      oi = oi+1
      outDat[oi,] = curAD
      oi = oi+1
    }
    
    
  }
  outDat<- outDat[ , !(names(outDat) %in% c('ID'))]
  return(outDat)
  
  
}






getROC <- function(biomark, groupID, criteria, flipLogic = F) {
  ROC = data.frame('crit' = criteria,
                   'hits' = criteria,
                   'misses'=criteria,
                   'CRs' = criteria,
                   'FAs' = criteria,
                   'TP' = criteria,
                   'TN' = criteria,
                   'FP' = criteria,
                   'FN' = criteria,
                   'accRaw' = criteria)
  
  if(flipLogic){
    criteria = rev(criteria)
    ROC$crit = criteria
  }

  AUC = 0
  for(ci in 1:length(criteria)){
    guess = rep(1,length(biomark))
    if(flipLogic){
      guess[biomark<criteria[ci]] = 2
    } else {
      guess[biomark>criteria[ci]] = 2
    }
    ROC$hits[ci] = sum(guess==2 & groupID==2)
    ROC$misses[ci] = sum(guess==1 & groupID==2)
    ROC$CRs[ci] = sum(guess==1 & groupID==1)
    ROC$FAs[ci] = sum(guess==2 & groupID==1)
    ROC$TP[ci] = sum(guess==2 & groupID==2) / sum(groupID==2)
    ROC$TN[ci] = sum(guess==1 & groupID==1) / sum(groupID==1)
    ROC$FP[ci] = sum(guess==2 & groupID==1) / sum(groupID==1)
    ROC$FN[ci] = sum(guess==1 & groupID==2) / sum(groupID==2)
    ROC$accRaw[ci] = (sum(guess==2 & groupID==2) + sum(guess==1 & groupID==1)) / length(groupID)
    if(ci>1){
      baseDist = ROC$FP[ci] - ROC$FP[ci-1]
      #top triangle area:
      tria = ((ROC$TP[ci] - ROC$TP[ci-1]) * baseDist) / 2
      #base rectangle area:
      reca = ROC$TP[ci-1] * baseDist
      AUC = AUC + tria + reca
    } else {
      baseDist = ROC$FP[ci]
      #top triangle area:
      tria = ((ROC$TP[ci]) * baseDist) / 2
      AUC = AUC + tria
    }

  }
  baseDist = 1 - ROC$FP[ci]
  #top triangle area:
  tria = ((1 - ROC$TP[ci]) * baseDist) / 2
  #base rectangle area:
  reca = ROC$TP[ci] * baseDist
  AUC = AUC + tria + reca
  ROC$AUC = AUC


  #sensitivity, true positive rate
  ROC$TPR = ROC$hits/(ROC$hits + ROC$misses)
  #fallout, false positive rate
  ROC$FPR = ROC$FAs/(ROC$FAs + ROC$CRs)
  #specificity, true negative rate
  ROC$TNR = ROC$CRs/(ROC$CRs + ROC$FAs)
  #accuracy
  ROC$acc = (ROC$TP + ROC$TN) / (ROC$TP + ROC$TN + ROC$FP + ROC$FN)
  #what's the optimal criterion?
  ci = which(ROC$acc == max(ROC$acc))
  ci = ci[1]

  baseRate = sum(groupID==2) / length(groupID)
  #bayesian posterior:
  post = ROC$TP[ci]*baseRate / (ROC$TP[ci]*baseRate + ROC$FP[ci]*(1-baseRate))
  ROC$base = baseRate
  ROC$post = post

  #false positive index
  guess = rep(1,length(biomark))
  guess[biomark>criteria[ci]] = 2
  ROC$FPi = sum(guess==2 & groupID==1) / sum(guess==2 & groupID==2)
  #precision, positive predictive value
  ROC$PPV = ROC$hits[ci] / (ROC$hits[ci] + ROC$FAs[ci])
  #negative predictive value
  ROC$NPV = ROC$CRs[ci] / (ROC$CRs[ci] + ROC$misses[ci])



  return(ROC)
}