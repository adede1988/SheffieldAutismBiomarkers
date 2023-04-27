function [] = setReadInAggregate(varargin)  %data, saveLoc, chanLocs

switch nargin
    case 4
        data = varargin{1}; 
        saveLoc = varargin{2}; 
        startVal = varargin{3}; 
        prefix = varargin{4};
        chanLocs = 5; 
        needChanLocs = 0; 
    case 5
        data = varargin{1}; 
        saveLoc = varargin{2}; 
        startVal = varargin{3};
        prefix = varargin{4};
        chanLocs = varargin{5}; 
        needChanLocs = 1; 
    otherwise
        "Error: incorrect number of inputs!"
end

%frequency params
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);



    disp(['subject: '  data.fn])
    tic
   
    %% 

  
    disp('reading in data')
    if needChanLocs == 1
        EEG = readEEGdat(join([data.dir data.fn],''), data.dataSet, startVal, chanLocs); 
    else
        EEG = readEEGdat(join([data.dir data.fn],''), data.dataSet, startVal); 
    end
   
    if EEG.probFlag == 0
    

        %% cleaning off noise epochs and channels
        [data.nbChanOrig, data.nbChanFinal, data.nbTrialOrig, data.nbTrialFinal,EEG]...
            = removeNoiseChansVolt(EEG, join([data.dir data.fn],''));

        %% interpolation to standard 32-channel montage 
        standardTrodes; 
        EEG = convertCoordinates(EEG, standardEEGlocs, data.dataSet); 
    
        %% do connectivity analysis since this is easier to do before splitting the participant
        [data.ispc, EEG] = getISPC(EEG, frex, numfrex, stds, join([data.dir data.fn],''));

        %% split the data into individual files for each channel to facilitate single channel HPC analysis
        disp('splitting the file')
        fileSplitter(EEG, join([data.dir data.fn],'')); 

       else
            disp(['subject: ' data.fn ' corrupt data!! time: ' num2str(toc)])
            data.donotsave = true; 
    end

if ~isfield(data, 'donotsave')
    disp(['trying to save' data.fn 'time: ' num2str(toc)])
    disp(saveLoc)
    if isfile(saveLoc)
        disp('summary Dat already saved')
    else
        save(saveLoc,'data')
    end
end

end

