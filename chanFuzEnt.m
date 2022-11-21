function [fuzEnt] = chanFuzEnt(out)
%input: 
%       out     a struct with the fields data and srate (not used)
%               assumes that data is a channels X time X epochs array
%               and that each epoch is 2000ms long regardless of sampling
%               rate. This is critical, because all data are resampled to
%               1000 Hz sampling rate, but the way the code is written,
%               this only works with the assumption of 2s long epochs. 

%output: 
%       fuzEnt Channels X epochs X scale; scale is hard coded to range
%       from 1 to 20. 

%function will first check to see if sampleEnt has already been calcualted
%for this data file. If yes, then don't waste time recalcualting it.
%Otherwise, do the calculation. 

%Note, the function sampEnt.m relies on functions MFE_mu.m and
%FuzEn.m. these functions were obtained from 
% Ref:
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", Medical & Biological Engineering &
% Computing, 2016.
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk

%this wrapper function was written by Adam Dede (adam.osman.dede@gmail.com)
%fall, 2022


if isfield(out, 'fuzEnt')
   fuzEnt = out.fuzEnt; 
   
else
  
   fuzEnt = nan([size(out.data,1), size(out.data,3), 20 ]);

   
   for snip = 1:size(out.data,3) %loop epochs
        curSnip = out.data(:,:,snip); 
        %resample to 1000 Hz sampling rate
        times = round(linspace(1,2000, size(out.data,2))); 
        curSnip = timeseries(curSnip', times ); 
        curSnip = resample(curSnip, [1:2000]);
        curSnip = curSnip.Data'; 
      
          
               
       fuzEnt(1,snip, :) = MFE_mu(curSnip, 2, .15, 2, 1, 20);
               
  

   end

  
 

end
end