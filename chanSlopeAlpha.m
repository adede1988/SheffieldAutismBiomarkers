function [slopeValsLog, slopeValsRel, alphaPeakLog, alphaPeakRel] = chanSlopeAlpha(out, frex)
%input: 
%       chanDat         a struct with the fields data, srate, and power 
%                       data is a channels X time X epochs array
%                       srate is a scalar (e.g. 500 or 1000)
%                       power is a channels X frequencies X epochs array of
%                       power values 
%       frex            vector of frequency values at which power has
%                       already been calculated

%output: 
%       slopeValsLog    1 X  2 matrix of slope and intercept
%                       fitted values resulting from fitting a line to the 
%                       1/f curve of the log transformed power spectrum
%       slopeValsRel    similar to slopeValsLog except resulting from a fit
%                       to the relative power spectrum
%       alphaPeakLog    scalar value of alpha peak frequency in Hz as found
%                       when estimated from the log transformed power
%                       spectrum
%       alphaPeakRel    Similar to alphaPeakLog except resulting from an
%                       estimate using the relative power spectrum

%NOTE: failures of estimation of the alpha peak frequency result in a
%returned value of 999 for either estimate. 

%calculations for 1/f slope fit follow those presented by: 
% Age-Related Changes in 1/f Neural Electrophysiological Noise
% Bradley Voytek, Mark A. Kramer, John Case, Kyle Q. Lepage, Zechari R. Tempesta, Robert T. Knight, Adam Gazzaley
% Journal of Neuroscience 23 September 2015, 35 (38) 13257-13265

%calculations for alpha peak frequency follow those presented by: 
% Dickinson A, DiStefano C, Senturk D, Jeste SS. Peak alpha frequency is a 
% neural marker of cognitive function across the autism spectrum. 
% Eur J Neurosci. 2018 Mar;47(6):643-651. doi: 10.1111/ejn.13645. 
% Epub 2017 Aug 1. PMID: 28700096; PMCID: PMC5766439.

%written by Adam Dede (adam.osman.dede@gmail.com)
%fall 2022


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

    
    
end
end