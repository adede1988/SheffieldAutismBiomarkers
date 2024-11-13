

library(tidyverse)
sampleSizes = seq(5,1000, 20)
reps = 100

results = data.frame(d = rep(0, length(sampleSizes)*reps), 
                     n = rep(0, length(sampleSizes)*reps))
ri = 1

for(sampSize in sampleSizes){ 
  for(ii in 1:reps){
    samp1 = rnorm(sampSize)
    samp2 = rnorm(sampSize)
    
    m1 = mean(samp1)
    m2 = mean(samp2)
    
    ss1 = sum((samp1 - m1)^2)
    ss2 = sum((samp2 - m2)^2)
    
    spool = sqrt((ss1+ss2) / (sampSize*2 - 2))
    
    d = (m1 - m2)
    results$d[ri] = d
    results$n[ri] = sampSize
    ri = ri+1
    
  }
  
  }

ggplot(results, aes(x = n, y = d)) + geom_point() + 
  geom_hline(yintercept = .2)

results %>% group_by(n) %>% summarize(prop.2 = sum(d>.2)/length(d)) %>%
  ggplot(aes(x = n, y = prop.2)) + geom_line() + 
  ylab('prop effect sizes greater than .2')

