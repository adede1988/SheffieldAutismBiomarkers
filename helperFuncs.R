
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


getModDat <- function(dat, varNames, tt, temp, ageGroups, ii){
  if(!grepl('phase', varNames[tt])){
    ageIdx = (dat$age > ageGroups[ii-1] & dat$age <= ageGroups[ii])
    ageIdx = which(ageIdx)
    keep = !is.na(dat[ageIdx,varNames[tt]] )
    keepi = which(keep)
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
    
    
    #what group has outliers? 
    outliers[tt,1] <<- sum(groupID[abs(dv_z)>5]==1)
    outliers[tt,2] <<- sum(groupID[abs(dv_z)>5]==2)
    outliers[tt,3] <<- sum(groupID[abs(dv_z)>5]==3)
    
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
  









getROC <- function(biomark, groupID, criteria) {
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
    guess[biomark>criteria[ci]] = 2
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