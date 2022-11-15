function [EEG] = readEEGdat(varargin)

%This function handles initial import of several formats of EEG data
%It can take in optional channel location information if that is provided
%from a separate file. It will sense the data file type depending on the
%file extension and use the appropriate EEGLAB readin function. It is
%assumed that EEGLAB has already been started 

%main actions performed: 
%The file type is detected
%The appropriate readin function is used
%If necessary the chanLocs from the input are aligned to the labels in the
%   loaded data.
%If necessary, epoch the data
%lowpass, highpass, and notch filter the data. 

%Adam Dede, Summer 2022

%% get out the inputs

if nargin == 3
    path = varargin{1}; 
    datSet = varargin{2}; 
    startVal = varargin{3}; 
    needLocs = false; 
else
    path = varargin{1}; 
    datSet = varargin{2}; 
    startVal = varargin{3}; 
    chanLocs = varargin{4}; 
    needLocs = true; 
end

%% detect file type and load in the EEG data structure
[folder,name,ext] = fileparts(path);

%local: 
% prefix = 'Z:/User/pc1aod/';
%HPC: 
% prefix = '/shared/dede_group/User/pc1aod/';
prefix = split(path, 'pc1aod'); 
prefix = [prefix{1}, 'pc1aod/']; 



switch datSet
    case 'KatChild'
            try % looking for previously processed data
                
                EEG = load(join([folder '/' name '_processed.mat'],'')); 
               
                disp('found previous work!')
           
                while isfield(EEG, 'EEG')
                    EEG = EEG.EEG; 
                end
                
            catch 
    
                EEG = pop_loadcnt(path, 'dataformat', 'auto', 'memmapfile', ''); 
                EEG=pop_chanedit(EEG, 'lookup', [prefix 'CODE/GEDbounds_clusterImprove/eeglab2022.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']);
                EEG=pop_select(EEG, 'nochannel', {'EKG', 'EMG'});
                
                restStart = EEG.event(find([EEG.event.type]==3004)).latency; 
                restEnd = EEG.event(find([EEG.event.type]==3006)).latency;
                EEG.data = EEG.data(:,find(EEG.times>restStart,1):find(EEG.times>restEnd,1)-1 ); 

            epochs = floor(size(EEG.data,2)/(2*EEG.srate));
            if epochs>100
                epochs = 100; 
            end
            EEG.data = EEG.data(:,1:epochs*(2*EEG.srate)); 
            EEG.data = reshape(EEG.data(:, 1:(2*EEG.srate)*epochs), [size(EEG.data,1),(2*EEG.srate),epochs ]);
            EEG.probFlag = 0; 
            EEG.times = [1:1000/EEG.srate:2000]; 

          

             %design notch filter
            wo = 50/(EEG.srate/2); 
            bw = wo/35; 
            [b,a] = iirnotch(wo,bw); 
          
            for seg = 1:epochs
                nextSeg = EEG.data(:,:,seg);
                nextSeg = nextSeg - mean(nextSeg,2); 
                padDat = flip(nextSeg,2);
                padDat = [padDat, nextSeg, padDat]; 
                test =  highpass(double(padDat)', .5, EEG.srate); 
                test = lowpass(test, 200, EEG.srate); 
                test = filtfilt(b,a,test); 
                EEG.data(:,:, seg) = test(EEG.srate*2+1:EEG.srate*4,:)';       
            end

            end

    case 'bpSZ'
        try % looking for previously processed data
         
            EEG = load(join([folder '/' name '_processed.mat'],'')); 
          
            while isfield(EEG, 'EEG')
                EEG = EEG.EEG; 
            end
            disp('found previous work!')
       
        catch % load raw data if not previously processed
            EEG = pop_loadcnt(path, 'dataformat', 'auto', 'memmapfile', ''); 
            EEG=pop_chanedit(EEG, 'lookup', [prefix 'CODE/GEDbounds_clusterImprove/eeglab2022.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']);
            EEG=pop_select(EEG, 'nochannel', {'EKG', 'EMG'});
            %eliminate channels without theta/radius info
            elim = zeros(size(EEG.data,1),1); 
            for chan=1:size(EEG.data,1)
                if isempty(EEG.chanlocs(chan).theta)
                    elim(chan) =  1; 
                end
            end
            EEG.chanlocs(elim==1) = [];
            EEG.data(elim==1, :) = []; 
            EEG.nbchan = size(EEG.data,1); 
            needLocs = false; 
            

            epochs = floor(size(EEG.data,2)/2000);
            if epochs>100
                epochs = 100; 
            end
            EEG.data = EEG.data(:,1:epochs*2000); 
            EEG.data = reshape(EEG.data(:, 1:2000*epochs), [size(EEG.data,1),2000,epochs ]);
            EEG.probFlag = 0; 
            EEG.times = [1:1:2000]; 
            


           %design notch filter
                wo = 60/(EEG.srate/2); 
                bw = wo/35; 
                [b,a] = iirnotch(wo,bw); 
              
                for seg = 1:epochs
                    nextSeg = EEG.data(:,:,seg);
                    nextSeg = nextSeg - mean(nextSeg,2); 
                    padDat = flip(nextSeg,2);
                    padDat = [padDat, nextSeg, padDat]; 
                    test =  highpass(double(padDat)', .5, EEG.srate); 
                    test = lowpass(test, 200, EEG.srate); 
                    test = filtfilt(b,a,test); 
                    EEG.data(:,:, seg) = test(EEG.srate*2+1:EEG.srate*4,:)';       
                end


        end

    case 'femaleASD'
        
        try % looking for previously processed data
           
            EEG = load(join([folder '/' name '_processed.mat'],'')); 
           
            while isfield(EEG, 'EEG')
                    EEG = EEG.EEG; 
            end
        catch % load raw data if not previously processed
            EEG = load(join([folder '/' name '.mat'],'')); 
            while isfield(EEG, 'EEG')
                    EEG = EEG.EEG; 
            end
            epochs = size(EEG.data,3); 
            triali = size(EEG.data,2); 


            %design notch filter
            wo = 60/(EEG.srate/2); 
            bw = wo/35; 
            [b,a] = iirnotch(wo,bw); 
          
            for seg = 1:epochs
                nextSeg = EEG.data(:,:,seg);
                nextSeg = nextSeg - mean(nextSeg,2); 
                padDat = flip(nextSeg,2);
                padDat = [padDat, nextSeg, padDat]; 
                test =  highpass(double(padDat)', .5, EEG.srate); 
                test = lowpass(test, 200, EEG.srate); 
                test = filtfilt(b,a,test); 
                EEG.data(:,:, seg) = test(triali+1:triali*2,:)';        
            end


            if ~isfield(EEG, 'probFlag')
                EEG.probFlag = 0; 
            end
        end

    case 'biomarkDev'
        try % looking for previously processed data
         
            
            EEG = load(join([folder '/' name '_processed.mat'],'')); 
            
            while isfield(EEG, 'EEG')
                EEG = EEG.EEG; 
            end
            disp('found previous work!')
        catch % load raw data if not previously processed
            delete(join([folder '/' name '_processed.mat'],'')); 
            EEG = ft_read_header(join([folder '/' name '.bin'],''));
            EEG = EEG.orig; 
            EEG.probFlag = 0; 

            Sessi = find(cellfun(@(x) strcmp(x, 'SESS'), {EEG.event.code}));
            SessType = {EEG.event(Sessi).label};
            restBool = cellfun(@(x) strcmp(x, 'Resting EEG'), SessType);
            restOn = EEG.event(Sessi(find(restBool, 1))).latency; 
            restOff = EEG.event(Sessi(find(restBool==0, 1))).latency; 


            try 



                dat = EEG.data; 
                dat = dat(:,[EEG.times>restOn & EEG.times<restOff]);

                epochs = floor(size(dat,2)/1000);
                if epochs>100
                    epochs = 100; 
                end
    
                EEG.data = dat(:,1:epochs*1000); 
                EEG.data = reshape(EEG.data, [size(EEG.data,1),1000,epochs ]);
                EEG.probFlag = 0; 
                EEG.times = [1:2:2000]; 
                EEG.data(125:128,:, :) = [];
                EEG.nbchan = 125; 
                EEG.chanlocs(125:128) = []; 
              
                %design notch filter
                wo = 60/(EEG.srate/2); 
                bw = wo/35; 
                [b,a] = iirnotch(wo,bw); 
              
                for seg = 1:epochs
                    nextSeg = EEG.data(:,:,seg);
                    nextSeg = nextSeg - mean(nextSeg,2); 
                    padDat = flip(nextSeg,2);
                    padDat = [padDat, nextSeg, padDat]; 
                    test =  highpass(double(padDat)', .5, EEG.srate); 
                    test = lowpass(test, 200, EEG.srate); 
                    test = filtfilt(b,a,test); 
                    EEG.data(:,:, seg) = test(EEG.srate*2+1:EEG.srate*4,:)';        
                end
             
            catch
                EEG.probFlag = 1; 
            end

         
            needLocs = false; 
            
          
            

        end


    case 'socBrain'
        try % looking for previously processed data
           
            
            EEG = load(join([folder '/' name '_processed.mat'],'')); 
          
            while isfield(EEG, 'EEG')
                EEG = EEG.EEG; 
            end
            disp('found previous work!')
        catch % load raw data if not previously processed
            EEG = pop_loadcnt(path, 'dataformat', 'auto', 'memmapfile', ''); 
            EEG=pop_chanedit(EEG, 'lookup', [prefix 'CODE/GEDbounds_clusterImprove/eeglab2022.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']);
            EEG=pop_select(EEG, 'nochannel', {'EKG', 'EMG'});
            %eliminate channels without theta/radius info
            elim = zeros(size(EEG.data,1),1); 
            for chan=1:size(EEG.data,1)
                if isempty(EEG.chanlocs(chan).theta)
                    elim(chan) =  1; 
                end
            end
            EEG.chanlocs(elim==1) = [];
            EEG.data(elim==1, :) = []; 
            EEG.nbchan = size(EEG.data,1); 
            needLocs = false; 
            EEG.probFlag = 0;
            try
                epochs = floor(size(EEG.data,2)/2000);
                EEG.data = reshape(EEG.data(:, 1:2000*epochs), [size(EEG.data,1),2000,epochs ]);
                 
                EEG.times = [1:2000]; 
               %design notch filter
                wo = 60/(EEG.srate/2); 
                bw = wo/35; 
                [b,a] = iirnotch(wo,bw); 
              
                for seg = 1:epochs
                    nextSeg = EEG.data(:,:,seg);
                    nextSeg = nextSeg - mean(nextSeg,2); 
                    padDat = flip(nextSeg,2);
                    padDat = [padDat, nextSeg, padDat]; 
                    test =  highpass(double(padDat)', .5, EEG.srate); 
                    test = lowpass(test, 200, EEG.srate); 
                    test = filtfilt(b,a,test); 
                    EEG.data(:,:, seg) = test(EEG.srate*2+1:EEG.srate*4,:)';        
                end
            catch
                EEG.probFlag = 1; 
            end

        end

    case 'biomarkCon'
        try % looking for previously processed data
           
            
            EEG = load(join([folder '/' name '_processed.mat'],'')); 
           
            while isfield(EEG, 'EEG')
                EEG = EEG.EEG; 
            end
            disp('found previous work!')
        catch % load raw data if not previously processed
           EEG = load(path); 
            while isfield(EEG, 'EEG')
                EEG = EEG.EEG; 
            end
       
            pnts = size(EEG.EEG_Resting,2); 
            dat = zeros(size(EEG.EEG_Resting,1), prod(size(EEG.EEG_Resting,[2,3])) );
            for block = 1:size(EEG.EEG_Resting,3)
                dat(:,(block-1)*size(EEG.EEG_Resting,2)+1: ... 
                    block*size(EEG.EEG_Resting,2) ) = ...
                    squeeze(EEG.EEG_Resting(:,:,block)); 
            end
            dat(:,pnts+1:end) = dat(:,pnts+1:end) - (dat(:,pnts+1) - dat(:,pnts)); 
            dat(:,pnts*2+1:end) = dat(:,pnts*2+1:end) - (dat(:,pnts*2+1) - dat(:,pnts*2)); 
            sRate = EEG.samplingRate; 

            EEG = load([prefix 'CODE/GEDbounds_clusterImprove/emptyEEG.mat']).EEG; 
            EEG.data = dat; 
            EEG.probFlag = 0; 
            try
                chanLocs(str2num(chanLocs.label)>124,:) = [];
                
                for chan=1:124
                    EEG.chanlocs(chan).X = 0; 
                end
                trodeNums = str2num(chanLocs.label); 
                for chan=1:124
                    loc = trodeNums(chan); 
                    lab = chanLocs.labels(chan);
                    EEG.chanlocs(loc).labels =  lab{1}; 
                    EEG.chanlocs(loc).theta = chanLocs.theta(chan); 
                    EEG.chanlocs(loc).radius = chanLocs.radius(chan); 
                    EEG.chanlocs(loc).X = chanLocs.X(chan); 
                    EEG.chanlocs(loc).Y = chanLocs.Y(chan); 
                    EEG.chanlocs(loc).Z = chanLocs.Z(chan); 
                    EEG.chanlocs(loc).trodeNum = loc; 
    
                end
                
                EEG.nbchan = 124; 
                if size(EEG.data,1)>EEG.nbchan
                    EEG.data(EEG.nbchan+1:end,:) = []; 
                end
                EEG.srate = sRate; 
                EEG.times = [1:2000]; 

                epochs = floor(size(EEG.data,2)/2000);
                EEG.data = reshape(EEG.data(:, 1:2000*epochs), [size(EEG.data,1),2000,epochs ]);
                 
               
                 %design notch filter
                wo = 60/(EEG.srate/2); 
                bw = wo/35; 
                [b,a] = iirnotch(wo,bw); 
              
                for seg = 1:epochs
                    nextSeg = EEG.data(:,:,seg);
                    nextSeg = nextSeg - mean(nextSeg,2); 
                    padDat = flip(nextSeg,2);
                    padDat = [padDat, nextSeg, padDat]; 
                    test =  highpass(double(padDat)', .5, EEG.srate); 
                    test = lowpass(test, 200, EEG.srate); 
                    test = filtfilt(b,a,test); 
                    EEG.data(:,:, seg) = test(EEG.srate*2+1:EEG.srate*4,:)';        
                end


                needLocs = false; 


            catch
                EEG.probFlag = 1; 
            end
            
            

        end


    otherwise
        warning('unrecognized type')
end



end