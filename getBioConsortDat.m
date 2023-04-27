
%local: 
% prefix = 'Z:/User/pc1aod\'; jobID = 3;  
%HPC: 
prefix = '/shared/dede_group/User/pc1aod/';

addpath([prefix 'CODE/GEDbounds_clusterImprove'])
addpath([prefix 'CODE/subNetworkDynamics'])
addpath([prefix 'CODE/export_fig_repo'])
addpath([prefix 'CODE/GEDbounds_clusterImprove/eeglab2022.0'])
savedir = [prefix "FIGURES/"];
summaryDatSave = [prefix 'SUMDAT/'];

%channel location information for any sets that need it
chanLocs = load([prefix 'CODE/GEDbounds_clusterImprove/emptyEEG.mat']);
chanLocs = chanLocs.EEG.chanlocs;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


%% NIMH Autism biomarkers consortium for clinical trials


data = struct; 

%metadata about participants contained in a .csv file
subInfo = readtable([prefix 'CODE/GEDbounds_clusterImprove/biomarkConDat.csv']);
%folder in which the raw EEG data can be found. These data can be
%downloaded from the NIMH data archive 
datFolder = [prefix 'Package_1202361/eeg_sub_files01/'];

%loop over participants and put metadata into a struct
for ii = 1:height(subInfo)
    %col1: subject key
    data(ii).key = subInfo.subjectkey{ii};
    %col2: dir (directory of the data)
    curFolder = subInfo.file{ii};
    [~,curFolder] = fileparts(curFolder);
    data(ii).dir = datFolder; 
    %col3: fn (file name of data)
    data(ii).fn = join([curFolder '_ORIG.mat'],'');
    %col4: sex
    data(ii).sex = subInfo.sex{ii};
    %col5: age
    data(ii).age = subInfo.age(ii);
    %col6: task
    data(ii).task = subInfo.task{ii};
    %col7: diag (group)
    data(ii).group = subInfo.diag{ii};
    %col8: iq1
    data(ii).iq1 = subInfo.iq1(ii);
    %col9: iq2
    data(ii).iq2 = subInfo.iq2(ii);
    %col10: iq1_measure
    data(ii).iq1_measure = subInfo.iq1_measure{ii};
    %col11: iq2_measure
    data(ii).iq2_measure = subInfo.iq2_measure{ii};
    %col12: ADOS
    data(ii).ADOS = subInfo.ADOS(ii);
    %col13: ADOS_version
    data(ii).ADOS_version = subInfo.ADOS_version(ii);
    %col14: eyes
    data(ii).eyes = subInfo.eyes{ii};
    %col15: dataSet
    data(ii).dataSet = subInfo.dataSet{ii};

end

%output message placed into the HPC output text file associated with this
%job
disp(join(['going for subs: ' num2str(jobID)  ' of ' num2str(length(data))],''))

%downselect to the current target subject. This set of metadata will serve
%as the initialization of the summary stats for this participant
data = data(jobID);


%EGI channel locations were not included with raw data, so a standard set
%of channel locations are used for all participants in this dataset
chanLocs = load([prefix 'CODE/GEDbounds_clusterImprove/EGIchanLocs3.mat']).chanLocs; 

%add in the channel numbers as labels! 
for chan =1:128
    temp = chanLocs.labels(chan);
    temp = temp{1}; 
    if isempty(temp)
        temp = cell(1,0); 
        temp{1} = chanLocs.label(chan,:);
        chanLocs.labels(chan) = temp; 
    end
end

%check that the subject has complete metadata and then send to the analysis
%pipeline
if ~isempty(data.sex)
setReadInAggregate(data, join([summaryDatSave 'biomarkCon_' num2str(jobID) '_' num2str(jobID) '.mat']), start,prefix, chanLocs);
end











