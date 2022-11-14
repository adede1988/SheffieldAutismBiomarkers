function [sampEnt] = chanEntropyVals(out)

if isfield(out, 'sampEnt')
   sampEnt = out.sampEnt; 
   fuzEnt = out.fuzEnt; 
   compSampEnt = out.compSampEnt; 
   compFuzEnt = out.compFuzEnt; 
else
  
   sampEnt = nan([size(out.data,1), size(out.data,3), 20 ]);
   fuzEnt = sampEnt; 
%    onsets = [1:size(out.dataSTD,2):prod(size(out.data,[2,3]))];
   
   for snip = 1:size(out.data,3) %loop epochs
        curSnip = out.data(:,:,snip); 
        %resample to 1000 Hz sampling rate
        times = round(linspace(1,2000, size(out.data,2))); 
        curSnip = timeseries(curSnip', times ); 
        curSnip = resample(curSnip, [1:2000]);
        curSnip = curSnip.Data'; 
      
          
               
       sampEnt(1,snip, :) = msentropy(curSnip,2, .3, 20);
               
  

   end

  
 

end
end