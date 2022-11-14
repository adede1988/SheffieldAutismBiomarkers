function [nbChanOrig, nbChanFinal, nbTrialOrig, nbTrialFinal, EEG] = removeNoiseChansVolt(EEG, path)

if ~isfield(EEG, 'noiseRemoved')
    datMean = squeeze(std(EEG.data,[],[1,2])); 
    EEG.data(:,:,datMean==0) = []; 
    
    nbChanOrig = size(EEG.data,1); 
    nbTrialOrig = size(EEG.data,3); 

    maxDeflection = zeros(size(EEG.data,1), size(EEG.data,3)); 
    window = 80/ (1000 / EEG.srate);
    for chan = 1:size(EEG.data,1)
        for trial = 1:size(EEG.data,3)
            for ii = 1:size(EEG.data,2) - window
                cur = max(EEG.data(chan,ii:ii+window, trial)) - min(EEG.data(chan,ii:ii+window, trial));
                if cur > maxDeflection(chan, trial)
                    maxDeflection(chan,trial) = cur; 
                end
            end
        end
    end
    
%     overallSTD = std(EEG.data,[], [1,2,3]);
%     overallMean = mean(EEG.data, [1,2,3]);
    
    %z-score it
%     maxDeflection = (maxDeflection - overallMean) ./ overallSTD; 


    %remove channels where over 50% of trials involve a max deflection of
    %greater than 100 microvolts
    noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);


    if length(noiseChans) > size(EEG.data,1)/4 %if over a quarter of channels are about to be removed, then try doing rough trial removal first
        %remove trials where over 75% of channels have 100 microvolt
        %deflections
        noiseTrials = find(sum(maxDeflection>100,1)> (size(EEG.data,1)*3/4) );
        if ~isempty(noiseTrials)
            EEG.data(:,:,noiseTrials) = []; 
            maxDeflection(:,noiseTrials) = []; 
        end
    
        nbChanFinal = size(EEG.data,1); 
        nbTrialFinal = size(EEG.data,3);
        % then go back and do channel removal
        noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);
        if ~isempty(noiseChans) 
             %checking to see if there are channels with super high power
             %anywhere
             EEG.chanlocs(noiseChans) = []; 
             EEG.nbchan = EEG.nbchan - length(noiseChans); 
             EEG.data(noiseChans,:,:) = []; 
             maxDeflection(noiseChans,:) = []; 
        end
         %remove trials where over 25% of channels have 100 microvolt
        %deflections
        noiseTrials = find(sum(maxDeflection>100,1)> (size(EEG.data,1)/4) );
        if ~isempty(noiseTrials)
            EEG.data(:,:,noiseTrials) = []; 
            maxDeflection(:,noiseTrials) = []; 
        end
    
        nbChanFinal = size(EEG.data,1); 
        nbTrialFinal = size(EEG.data,3);

    else 

    
        if ~isempty(noiseChans) 
             %checking to see if there are channels with super high power
             %anywhere
             EEG.chanlocs(noiseChans) = []; 
             EEG.nbchan = EEG.nbchan - length(noiseChans); 
             EEG.data(noiseChans,:,:) = []; 
             maxDeflection(noiseChans,:) = []; 
        end
         %remove trials where over 25% of channels have 100 microvolt
        %deflections
        noiseTrials = find(sum(maxDeflection>100,1)> (size(EEG.data,1)/4) );
        if ~isempty(noiseTrials)
            EEG.data(:,:,noiseTrials) = []; 
            maxDeflection(:,noiseTrials) = []; 
        end
    
        nbChanFinal = size(EEG.data,1); 
        nbTrialFinal = size(EEG.data,3);
    end

%% rereference the data to average
   
    EEG.data = EEG.data - mean(EEG.data, 1); 

   


    [folder,name,~] = fileparts(path);
    EEG.noiseRemoved = true; 
    EEG.nbChanOrig = nbChanOrig; 
    EEG.nbTrialOrig = nbTrialOrig;
    save(join([folder '/' name '_processed.mat'],''), "EEG")
else
    nbChanFinal = size(EEG.data,1); 
    nbTrialFinal = size(EEG.data,3);
    try
        nbChanOrig = EEG.nbChanOrig; 
        nbTrialOrig = EEG.nbTrialOrig; 
    catch
        nbChanOrig = 999; 
        nbTrialOrig = 999; 
    end
end

end