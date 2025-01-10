
source("G:/My Drive/GitHub/MasterStatsUsingR/courseTheme.R")

figPath = 'G:\\My Drive\\Milne\\clusterDifferences\\ageEffects\\'

for( var in varNames){
  print(var)
  test = lm(logPow_beta_lF ~ poly(age,3), data = dat)
  
  if (sum(summary(test)$coefficients[3:4, "Pr(>|t|)"]<.05)>0){
  
    
    dv = dat[[var]]
    if(is.character(dv)){
      dv <- as.numeric(dv)
    }
    age = dat$age[!is.na(dv)]
    dv = dv[!is.na(dv)]
    lims = quantile(dv, c(.1, .90))
    dv_mean = mean(dv[dv>lims[1] & dv<lims[2]])
    dv_sd = sd(dv[dv>lims[1] & dv<lims[2]])
    dv_z = (dv - dv_mean) / dv_sd
    nonOut = which(abs(dv_z)<5)
    dv_mean = mean(dv[nonOut])
    dv_sd = sd(dv[nonOut])
    dv_z = (dv - dv_mean) / dv_sd
    
    
    dv = dv[abs(dv_z)<5]
    age = age[abs(dv_z)<5]
   
    df = data.frame(x = age, y = dv)
    
    
  
  outplot <- df %>% ggplot(aes(x = x, y = y)) + 
    geom_point(size = 3, alpha = .5, color = '#4ED6B2') + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), 
                linewidth = 3,
                color = "#FB9564", se = FALSE)+
    xlab('age') +
    ylab(var) + 
    myTheme 
  
    png(paste(figPath, 'AgeScatter_', var, '.png',sep=''),         # File name
      width=1024, height=768)
    print(outplot)
    dev.off()
  }
  
  
}