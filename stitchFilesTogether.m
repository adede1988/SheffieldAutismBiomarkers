%data combine


prefix = 'Z:/User/pc1aod/'; %subNum = 1; chan = 1; 
%HPC: 
% prefix = '/shared/dede_group/User/pc1aod/';


addpath([prefix 'CODE/GEDbounds_clusterImprove'])
addpath([prefix 'CODE/subNetworkDynamics'])
addpath([prefix 'CODE/export_fig_repo'])
savedir = [prefix 'FIGURES/'];
summaryDatSave = [prefix 'SUMDAT/'];

filenames = dir(fullfile(summaryDatSave,'*.mat'));

%frequency params
frex = logspace(log10(2),log10(80),100);
transforms = [1,2,3,4]; 
numfrex = length(frex); 
stds = linspace(2,5,numfrex);


parfor ii = 1:length(filenames)
    tic
    filenames(ii).errorFlag = 0; 
    data = load([filenames(ii).folder '/' filenames(ii).name]).data;

    switch data.dataSet
        case 'biomarkDev'
            fn = split(data.fn, '.bin'); 
        case 'biomarkCon'
            fn = split(data.fn, '.mat'); 
        case 'bpSZ'
            fn = split(data.fn, '.cnt'); 
        case 'femaleASD'
            fn = split(data.fn, '.mat'); 
        case 'KatChild'
            fn = split(data.fn, '.cnt'); 
        case 'socBrain'
            fn = split(data.fn, '.cnt'); 
        otherwise
            fn = split(data.fn, '.mat'); 
    end
    
    folderLoc = data.dir;
    
    folderLoc = split(folderLoc, 'pc1aod');
    folderLoc = [prefix folderLoc{2}];


    trialCount = data.nbTrialFinal;
    missingData = false; 
    %preallocate power variable
    if ~isfield(data, 'power')
        
        data.power = zeros(32,numfrex,trialCount); 
        missingData = true; 
    end
    
    %preallocate PAC variables
    if ~isfield(data, 'MINull')
        data.PAC = zeros(10,21,18,32); %low freq X high freq X phase bins X channels
        data.PACMI = cell(32,1); 
        data.PACsig = cell(32,1); 
        data.phasePref = cell(32,1); 
        data.MINull = cell(32,1); 
        missingData = true; 
      
    end

    %preallocate entropy variable
    if ~isfield(data, 'sampEnt')
        data.sampEnt = zeros(32,trialCount, 20); 
        missingData = true; 
    end

    if ~isfield(data, 'fuzEnt')
        data.fuzEnt = zeros(32, trialCount, 20); 
        missingData = true; 
    end

    %preallocate 1/f slope and peak alpha fits

    if ~isfield(data, 'slopeValsLog')
        data.slopeValsLog = zeros(32,2); 
        data.slopeValsRel = zeros(32,2); 
        data.alphaPeakLog = zeros(32,1); 
        data.alphaPeakRel = zeros(32,1); 
        missingData = true; 
    end





    if missingData %no point in opening all the channel files if there's no missing data
    for chan = 1:32
        try %try for the channel file in general
            out = load([folderLoc fn{1} '_processed_temp' num2str(chan) '.mat']).out; 
           
            
    %% try for power values
            try
                data.power(chan,:,:) = out.power; 
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan]; %+0 give errors for power
            end
    %% try for PAC values
            try 
                if data.PAC(1,1,1,chan) == 0 %do we need to fill in the PAC values?
                    data.PAC(:,:,:,chan) = out.PAC(:,:,:,chan);
                    data.PACMI{chan} = out.PACMI{chan};
                    data.PACsig{chan} = out.PACsig{chan};
                    data.phasePref{chan} = out.phasePref{chan};
                    data.MINull{chan} = out.MINull{chan};
                end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+100]; %+100 give errors for PAC
            end
            
    %% try for sample entropy values
            try 
                if data.sampEnt(chan,1,1) == 0 
                    data.sampEnt(chan,:,:) = out.sampEnt; 
                end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+200]; %+200 give errors for sampEnt
            end

    %% try for 1/f slopes and peak alpha values
            try
               if data.slopeValsLog(chan,1)==0
                   data.slopeValsLog(chan,:) = out.slopeValsLog; 
                   data.slopeValsRel(chan,:) = out.slopeValsRel; 
                   data.alphaPeakLog(chan) = out.alphaPeakLog; 
                   data.alphaPeakRel(chan) = out.alphaPeakRel; 
               end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+300]; %+300 give errors for slope/alpha
            end

    %% try for fuzzy entropy 
            try 
                if data.fuzEnt(chan,1,1) == 0 
                    data.fuzEnt(chan,:,:) = out.fuzEnt; 
                end
            catch
                filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+600]; %+600 give errors for fuzEnt
            end

        catch
        filenames(ii).errorFlag = [filenames(ii).errorFlag, chan+400]; %+400 gives global error
        
        end


    end
    end

  

    %% do a check for ISPC
    if ~isfield(data, 'ispc')
        filenames(ii).errorFlag = [filenames(ii).errorFlag, 500]; %500 means ispc missing
    end

    


% save the summary data file back out
parSave([filenames(ii).folder '/' filenames(ii).name], data); 
% save([filenames(ii).folder '/' filenames(ii).name], 'data');


disp(['subject: ' num2str(ii) ' errors: ' num2str(length(filenames(ii).errorFlag)-1) ' time: ' num2str(round(toc))])

end






