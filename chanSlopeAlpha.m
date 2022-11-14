function [slopeValsLog, slopeValsRel, alphaPeakLog, alphaPeakRel, out] = chanSlopeAlpha(out, frex)


warning('off','all')
if isfield(out, 'slopeValsLog')
    slopeValsLog = out.slopeValsLog;
    slopeValsRel = out.slopeValsRel;
    alphaPeakLog = out.alphaPeakLog;
    alphaPeakRel = out.alphaPeakRel;
else
    allPow = out.power; 
    logPow = log10(allPow) ; 
    relPow = allPow ./ sum(allPow,2);
    slopeValsLog = zeros(size(logPow,1),2);
    slopeValsRel = zeros(size(logPow,1),2);
    alphaPeakLog = zeros(size(logPow,1),1);
    alphaPeakRel = zeros(size(logPow,1),1);
    opt = statset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
    
    for chan = 1:size(logPow,1)
        try
        modelfun = @(b,x) b(1)*x + b(2);
        beta0 = [-1, 1];
        modPow = logPow(chan, frex<7 | (frex>14 & frex<24) );
        modFrex = frex(frex<7 | (frex>14 & frex<24));
        mdl = fitnlm(modFrex,modPow, modelfun, beta0);
        slopeValsLog(chan,:) = mdl.Coefficients.Estimate; 
        catch
            slopeValsLog(chan,:) = [999, 999]; 
        end
        try
        modPow = relPow(chan, frex<7 | (frex>14 & frex<24) );
        modFrex = frex(frex<7 | (frex>14 & frex<24));
        mdl2 = fitnlm(modFrex,modPow, modelfun, beta0, 'options', opt);
        slopeValsRel(chan,:) = mdl2.Coefficients.Estimate;
        catch
            slopeValsRel(chan,:) = [999, 999]; 
        end
        
        try
        %find alpha peak at all electrodes with log transformed power
        afrex = frex(frex>6 & frex<14);
        gausMod = @(b,x) b(1)*exp(-( (x-b(2)).^2 ./ (2*b(3)^2) ) );
        beta0 = [1,10,2]; 
 
        detrendDat = logPow(chan,frex>6 & frex<14) -...
            modelfun(mdl.Coefficients.Estimate, frex(frex>6 & frex<14)); 
        gausmdl = fitnlm([1:length(detrendDat)], detrendDat, gausMod, beta0, 'options', opt); 
        if round(gausmdl.Coefficients.Estimate(2)) > 0 && round(gausmdl.Coefficients.Estimate(2)) < length(afrex)
            alphaPeakLog(chan) = afrex(round(gausmdl.Coefficients.Estimate(2)));
        else
         
            alphaPeakLog(chan) = 999; 
        end
        catch
            alphaPeakLog(chan) = 999; 
        end
        %find alpha peak at all electrodes with relative power

        try
        detrendDat = relPow(chan,frex>6 & frex<14) -...
            modelfun(mdl2.Coefficients.Estimate, frex(frex>6 & frex<14)); 
        gausmdl = fitnlm([1:length(detrendDat)], detrendDat, gausMod, beta0, 'options', opt); 
        if round(gausmdl.Coefficients.Estimate(2)) > 0 && round(gausmdl.Coefficients.Estimate(2)) < length(afrex)
            alphaPeakRel(chan) = afrex(round(gausmdl.Coefficients.Estimate(2)));
        else %what to do if there's no available alpha peak! 
          
            
            alphaPeakRel(chan) = 999; 
        end
        catch
            alphaPeakRel(chan) = 999; 
        end
    end  
    out.slopeValsLog =  slopeValsLog;
    out.slopeValsRel = slopeValsRel;
    out.alphaPeakLog = alphaPeakLog;
    out.alphaPeakRel = alphaPeakRel;
    
    
end
end