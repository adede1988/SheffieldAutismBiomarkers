

%% recreating topo plots: 


datLoc = 'F:\SUMDAT';
codeLoc = 'G:\My Drive\GitHub\SheffieldAutismBiomarkers';
addpath(codeLoc)
sumDat = dir(datLoc);
sumDat(1:2) = []; 
frex = logspace(log10(2),log10(80),100);
ageGroupCuts = [0,97,127,250];
EEGlocs = readtable([codeLoc '/' 'standardEEGlocs.csv']);
theta = EEGlocs.theta+180; 
rd = cos(pi*(EEGlocs.phi)/180);
%store the values for generating a topo plot
%channels X subjects
allAlphaTopo = zeros([32,length(sumDat)]);
%store age group info: 
allAgeGroups = zeros([1, length(sumDat)]); 
%store diag group info: 
%1 = CON, 2 = ASD, 3 = AD
allDiagGroups = zeros([1, length(sumDat)]); 

allErrors = zeros(length(sumDat), 1); 
parfor ii = 1:length(sumDat)
    
    try
    curDat = load([sumDat(ii).folder '/' sumDat(ii).name]).data;
    age = curDat.age; 
    group = curDat.group;
   
       
        pow = curDat.power; 
        relPow = pow ./ sum(pow,2); 
        logPow = log10(pow); 

        psd = squeeze(mean(pow,3));
        psdRel = squeeze(mean(relPow,3));
        psdLog = squeeze(mean(logPow,3));
    
       

    %store channel specific psds 
    for chan = 1:32
        %raw
        T = table(psd(chan,:)', frex');

        % name columns
        T.Properties.VariableNames = {'power', 'frex'};

        % Save the table as a CSV file
        writetable(T, ['F:\testPSDs2/' curDat.key ...
            '_' curDat.eyes '_rawPow_' num2str(chan) '.csv']);

        %rel
        T = table(psdRel(chan,:)', frex');

        % name columns
        T.Properties.VariableNames = {'power', 'frex'};

        % Save the table as a CSV file
        writetable(T, ['F:\testPSDs2/' curDat.key ...
            '_' curDat.eyes '_relPow_' num2str(chan) '.csv']);

        %log
        T = table(psdLog(chan,:)', frex');

        % name columns
        T.Properties.VariableNames = {'power', 'frex'};

        % Save the table as a CSV file
        writetable(T, ['F:\testPSDs2/' curDat.key ...
            '_' curDat.eyes '_logPow_' num2str(chan) '.csv']);

    end


    catch
        disp(ii)
        allErrors(ii) = 1; 
    end
end
tempAlpha = allAlphaTopo; 
tempAge = allAgeGroups; 
tempDiag = allDiagGroups;
allAlphaTopo(:,allAlphaTopo(1,:)==0) = []; 
allAgeGroups(allAgeGroups ==0) = []; 
allDiagGroups(allDiagGroups==0) = []; 

ageGroupVals = unique(allAgeGroups); 
diagGroupVals = unique(allDiagGroups); 
diagLabs = {'CON', 'ASD', 'AD'}; 

for ai = 1:3
    for gi = 1:3
        curAlpha = allAlphaTopo(:,allAgeGroups==ageGroupVals(ai) &...
                                    allDiagGroups == diagGroupVals(gi)); 
        figure('position', [0,0,500,510], 'visible', false)
        quickTopo2(squeeze(mean(curAlpha,2)), EEGlocs.theta+180, ...
            cos(pi*(EEGlocs.phi)/180), 1);
        colorbar
        if ai == 1
            caxis([.005, .007])
            caxis([.003, .008])
        elseif ai == 2
            caxis([.007, .011])
            caxis([.007, .012])
        else 
            caxis([.009, .013])
            caxis([.009, .014])
        end
        % title([diagLabs(gi) ' age: ' num2str(ai)])
        set(gcf, 'Color', 'w');
        export_fig(fullfile('G:', 'My Drive',...
            'Milne', 'pubFigPanels', ...
     ['Fig1_'  diagLabs{gi} '_' num2str(ai) '.jpg']),...
        '-r300')


    end
end