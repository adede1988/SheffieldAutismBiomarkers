


%% batting quick clean up


%local: 
prefix = 'Z:/User/pc1aod/'; 
%HPC: 
% prefix = '/shared/dede_group/User/pc1aod/';


addpath([prefix 'CODE/GEDbounds_clusterImprove'])
addpath([prefix 'CODE/subNetworkDynamics'])
addpath([prefix 'CODE/export_fig_repo'])
savedir = [prefix 'FIGURES/'];
summaryDatSave = [prefix 'SUMDAT/'];

filenames = dir(fullfile(summaryDatSave,'*.mat'));
% PACmissing = load([ prefix 'CODE/GEDbounds_clusterImprove/pacMissing221114.mat']).PACmissing;
% targets = find(auditResults==1); 
PACmissing = zeros(88160,1); 
errors = zeros(88160,1); 
parfor inputNum = 1:88160
    disp(inputNum)
%     if PACmissing(inputNum) == 1
        try
            PACmissing(inputNum) = PACaudit(filenames, inputNum, prefix);
        catch
            PACmissing(inputNum) = 3; %some reason there was an error in even setting up to load 
        end
   
%     end
end

