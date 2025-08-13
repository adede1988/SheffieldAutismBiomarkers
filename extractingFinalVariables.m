% now I need to extract key variables from individual participants and
% store it into a final spreadsheet.





prefix = 'F:\'; %subNum = 1; chan = 1; 
%HPC: 
% prefix = '/shared/dede_group/User/pc1aod/';


addpath([prefix 'CODE/GEDbounds_clusterImprove'])
addpath([prefix 'CODE/subNetworkDynamics'])
addpath([prefix 'CODE/export_fig_repo'])
savedir = [prefix 'FIGURES/'];
summaryDatSave = [prefix 'SUMDAT/'];

filenames = dir(fullfile(summaryDatSave,'*.mat'));
allFooofRes = dir(fullfile(prefix, 'testPSDs2\Complete', '*.csv'));
splitNames = cellfun(@(x) strsplit(x, '_'), {allFooofRes.name},...
    'UniformOutput',false)';
for ii =1:length(allFooofRes)
    nameBits = splitNames{ii};
    allFooofRes(ii).ID = nameBits{1}; 
    allFooofRes(ii).eyes = nameBits{2}; 
    allFooofRes(ii).type = nameBits{3}; 
    chan = strsplit(nameBits{4}, '.csv'); 
    allFooofRes(ii).chan = str2double(chan{1});
end
%frequency params
frex = logspace(log10(2),log10(80),100);
transforms = [1,2,3,4]; 
numfrex = length(frex); 
stds = linspace(2,5,numfrex);

exampleDat = load('F:\SUMDAT\biomarkCon_1_1.mat').data; 
labels = exampleDat.labels; 
Th = exampleDat.Th; 
Rd = exampleDat.Rd; 


%regions to average 
rightFrontal = ["Fp2", "AF4", "F4", "F8"]; 
leftFrontal = ["Fp1", "AF3", "F3", "F7"]; 
rightCentroParietal = ["FC2", "FC6", "C4", "CP2", "CP6"];
leftCentroParietal = ["FC1", "FC5", "C3", "CP1", "CP5"];
rightOccipitoParietal = ["P4", "P8", "PO4", "O2"];
leftOccipitoParietal = ["P3", "P7", "PO3", "O1"];
frontal = ["Fp1", "Fp2", "AF3", "AF4", "F4", "Fz", "F3"];
occipital = ["PO4", "PO3", "O2", "Oz", "O1"];
central = ["Fz", "Cz", "Pz", "Oz"]; 
leftLateral = ["F7", "FC5", "T7", "CP5", "P7"];
rightLateral = ["F8", "FC6", "T8", "CP6", "P8"];
rightHemisphere = ["Fp2", "AF4", "F4", "F8", "FC6", "FC2", "T8", "C4", "CP6", "CP2", "P8", "P4", "PO4", "O2"]; 
leftHemisphere = ["Fp1", "AF3", "F3", "F7", "FC5", "FC1", "T7", "C3", "CP5", "CP1", "P7", "P3", "PO3", "O1"];

%regionNames
regNames = {'rF', 'lF', 'rCP', 'lCP', 'rOP', 'lOP', 'F', 'O', 'C', 'lL', 'rL', 'rH', 'lH'}; 
regVars = {rightFrontal, leftFrontal, rightCentroParietal, leftCentroParietal, rightOccipitoParietal, leftOccipitoParietal,...
    frontal, occipital, central, leftLateral, rightLateral, rightHemisphere, leftHemisphere}; 


%asymmetries
interhemispheric = {["Fp1", "F3", "F7", "C3", "T7", "P3", "P7", "O1"], ["Fp2", "F4", "F8", "C4", "T8", "P4", "P8", "O2"]};
leftRostroCaudal = {["O1", "P3", "P7", "C3", "T7", "CP1"], ["P3", "C3", "T7", "F3", "F7", "FC1"]};
rightRostroCaudal = {["O2", "P4", "P8", "C4", "T8", "CP2"], ["P4", "C4", "T8", "F4", "F8", "FC2"]};
leftmedioLateral = {["P7", "CP5", "T7", "FC5", "F7"], ["P3", "CP1", "C3", "FC1", "F3"]};
rightmedioLateral= {["P8", "CP6", "T8", "FC6", "F8"], ["P4", "CP2", "C4", "FC2", "F4"]};

%asym names
asymNames = {'iH', 'lRC', 'rRC', 'lML', 'rML'}; 
asymVars = {interhemispheric, leftRostroCaudal, rightRostroCaudal, leftmedioLateral, rightmedioLateral}; 


%frequency bounds greater than lower through to less than or equal to upper
freqBounds = [0, 19, 38, 53, 73, 87, 100];
freqNames = {'delta', 'theta', 'alpha', 'beta', 'gam1', 'gam2'};

%graining boundaries for sample entropy
grainBounds = [0, 7, 13, 20];
grainNames = {'fine', 'med', 'course'}; 

%distances for ISPC
 dists = zeros(32); 
 [y,x] = pol2cart(Th,Rd); 
    for chan1 = 1:32
    for chan2 = 1:32
        dists(chan1, chan2) = sqrt( (y(chan1) - y(chan2))^2 + (x(chan1) - x(chan2))^2 );
    end
    end
distVals = triu(dists,1); 
distVals = distVals(:); 
distVals(distVals==0) = []; 

%frequencies for PAC
LF_pac = [0, 2, 4, 7, 10];
HF_pac = [0, 4, 9, 21];
LF_names = {'del', 'the', 'al', 'bet'}; 
HF_names = {'bet', 'gL', 'gH'}; 

errorLog = []; 

for ii = 1:length(filenames)
    ii
  try  
  %% load the file
  data = load([filenames(ii).folder '/' filenames(ii).name]).data;

  %% get fooof based 1/f stats
    idx = cellfun(@(x, y) strcmp(x, data.key) & strcmp(y, data.eyes), ...
                    {allFooofRes.ID}, {allFooofRes.eyes}); 
    curFooof = allFooofRes(idx,:); 
    if length(curFooof) ~= 32
        errorLog =[errorLog, ii]; 
    else
        for chan = 1:32
            T = readtable([curFooof(chan).folder '/' curFooof(chan).name]);
            data.slopeValsLog(curFooof(chan).chan, 1) = T.exponent(1); 
            data.slopeValsLog(curFooof(chan).chan, 2) = T.offset(1); 
        end
    end
  %% basic demographics
  %data set
  filenames(ii).dataSet = data.dataSet; 
  %sex
  filenames(ii).sex = data.sex; 
  %age
  filenames(ii).age = data.age; 
  %eyes
  filenames(ii).eyes = data.eyes; 
  %group
  filenames(ii).group = data.group; 
  %sub ID
  filenames(ii).key = data.key; 
  %original file name
  filenames(ii).origfile = data.fn; 
  %original trial count
  filenames(ii).nbTrialOrig = data.nbTrialOrig; 
  %final trial count
  filenames(ii).nbTrialFinal = data.nbTrialFinal; 
  %original chan count
  filenames(ii).nbChanOrig = data.nbChanOrig; 
  %final chan count
  filenames(ii).nbChanFinal = data.nbChanFinal; 
    

  switch data.dataSet
      case 'bpSZ'
          filenames(ii).IQ = data.iq1; 
          filenames(ii).IQ_measure = data.iq1_measure; 
      case 'femaleASD'
          filenames(ii).IQ = data.iq1; 
          filenames(ii).IQ_measure = data.iq1_measure;
      case 'biomarkCon'
          if isnan(data.iq1)
            filenames(ii).IQ = data.iq2;
            filenames(ii).IQ_measure = data.iq2_measure; 
          else
            filenames(ii).IQ = data.iq1;
            filenames(ii).IQ_measure = data.iq1_measure; 
          end
      case 'socBrain'
          if isnan(data.iq1)
            filenames(ii).IQ = data.iq2;
            filenames(ii).IQ_measure = data.iq2_measure; 
          else
            filenames(ii).IQ = data.iq1;
            filenames(ii).IQ_measure = data.iq1_measure; 
          end
      case 'KatChild'
          filenames(ii).IQ = data.iq1; 
          filenames(ii).IQ_measure = data.iq1_measure;
      case 'biomarkDev'
          if isnan(data.iq2)
            filenames(ii).IQ = data.iq1;
            filenames(ii).IQ_measure = data.iq1_measure; 
          else
            filenames(ii).IQ = data.iq2;
            filenames(ii).IQ_measure = data.iq2_measure; 
          end
  end

  %% regional power 
  pow = data.power; 
  relPow = pow ./ sum(pow,2); 
  logPow = log10(pow); 
  
  for reg = 1:length(regNames)
  for fi = 2:7
        regVar = regVars{reg}; 
        idx = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), regVar); 
        meanPow = mean(pow(idx, freqBounds(fi-1)+1:freqBounds(fi), : ), 'all' );
        filenames(ii).(['pow_' freqNames{fi-1} '_' regNames{reg}]) = meanPow; 
        meanPow =  mean(relPow(idx, freqBounds(fi-1)+1:freqBounds(fi), : ), 'all' );
        filenames(ii).(['relPow_' freqNames{fi-1} '_' regNames{reg}]) = meanPow; 
        meanPow =  mean(logPow(idx, freqBounds(fi-1)+1:freqBounds(fi), : ), 'all' );
        filenames(ii).(['logPow_' freqNames{fi-1} '_' regNames{reg}]) = meanPow; 
  end
  end

  %% power asymmetries 
  for asi = 1:length(asymNames)
  for fi = 2:7  
        asVar = asymVars{asi}; 
        idx1 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{1}); 
        idx2 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{2}); 
        meanPow = mean((pow(idx1, freqBounds(fi-1)+1:freqBounds(fi), : ) - ...
                       pow(idx2, freqBounds(fi-1)+1:freqBounds(fi), : )) ./ ...
                       (pow(idx1, freqBounds(fi-1)+1:freqBounds(fi), : ) + ...
                       pow(idx2, freqBounds(fi-1)+1:freqBounds(fi), : )), 'all' );
        filenames(ii).(['powAsym_' freqNames{fi-1} '_' asymNames{asi}]) = meanPow;

        meanPow = mean((relPow(idx1, freqBounds(fi-1)+1:freqBounds(fi), : ) - ...
                       relPow(idx2, freqBounds(fi-1)+1:freqBounds(fi), : )) ./ ...
                       (relPow(idx1, freqBounds(fi-1)+1:freqBounds(fi), : ) + ...
                       relPow(idx2, freqBounds(fi-1)+1:freqBounds(fi), : )), 'all' );
        filenames(ii).(['relPowAsym_' freqNames{fi-1} '_' asymNames{asi}]) = meanPow; 

        meanPow = mean((logPow(idx1, freqBounds(fi-1)+1:freqBounds(fi), : ) - ...
                       logPow(idx2, freqBounds(fi-1)+1:freqBounds(fi), : ))./ ...
                       (logPow(idx1, freqBounds(fi-1)+1:freqBounds(fi), : ) + ...
                       logPow(idx2, freqBounds(fi-1)+1:freqBounds(fi), : )), 'all' );
        filenames(ii).(['logPowAsym_' freqNames{fi-1} '_' asymNames{asi}]) = meanPow; 
  end
  end

  %% regional 1/f slope fit
  logSlope = data.slopeValsLog(:,1);
  relSlope = data.slopeValsRel(:,1); 

    for reg = 1:length(regNames)
    
        regVar = regVars{reg}; 
        idx = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), regVar); 
        meanSlope = mean(logSlope(idx));
        filenames(ii).(['logSlope_' regNames{reg}]) = meanSlope; 
        meanSlope = mean(relSlope(idx));
        filenames(ii).(['relSlope_' regNames{reg}]) = meanSlope; 

    end

  %% 1/f slope asymmetries
 for asi = 1:length(asymNames)
 
        asVar = asymVars{asi}; 
        idx1 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{1}); 
        idx2 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{2}); 
        meanSlope = mean((logSlope(idx1) - logSlope(idx2)) ./ (logSlope(idx1) + logSlope(idx2)));
        filenames(ii).(['logSlopeAsym_' asymNames{asi}]) = meanSlope; 

        meanSlope = mean((relSlope(idx1) - relSlope(idx2)) ./ (logSlope(idx1) + logSlope(idx2)));
        filenames(ii).(['relSlopeAsym_' asymNames{asi}]) = meanSlope; 
 end

  %%  regional peak alpha frequency
  logAlpha = data.alphaPeakLog; 
  relAlpha = data.alphaPeakRel; 

  filenames(ii).(['missedLogAlpha']) = sum(logAlpha==999);
  filenames(ii).(['missedRelAlpha']) = sum(relAlpha==999);


  %nan out the 999 values
  logAlpha(logAlpha == 999) = nan; 
  relAlpha(relAlpha == 999) = nan; 
    for reg = 1:length(regNames)
    
        regVar = regVars{reg}; 
        idx = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), regVar); 
        meanAlpha = mean(logAlpha(idx), 'omitnan');
        filenames(ii).(['logAlpha_' regNames{reg}]) = meanAlpha; 
        meanAlpha = mean(relAlpha(idx), 'omitnan');
        filenames(ii).(['relAlpha_' regNames{reg}]) = meanAlpha; 

    end

  %% alpha asymmetries
 for asi = 1:length(asymNames)
 
        asVar = asymVars{asi}; 
        idx1 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{1}); 
        idx2 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{2}); 
        meanAlpha = mean((logAlpha(idx1) - logAlpha(idx2)) ./ (logAlpha(idx1) + logAlpha(idx2)), 'omitnan');
        filenames(ii).(['logAlphaAsym_' asymNames{asi}]) = meanAlpha; 
        meanAlpha = mean((relAlpha(idx1) - relAlpha(idx2)) ./ (relAlpha(idx1) + relAlpha(idx2)), 'omitnan');
        filenames(ii).(['relAlphaAsym_' asymNames{asi}]) = meanAlpha; 
 end


 %% regional sample entropy  
  sampEnt = data.sampEnt; 
% %graining boundaries for sample entropy
% grainBounds = [0, 7, 13, 20];
% grainNames = {'fine', 'med', 'course'}; 

  for reg = 1:length(regNames)
      %get all as a separate thing
      regVar = regVars{reg}; 
      idx = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), regVar);
      meanSE = mean(sampEnt(idx, :,  :), 'all' );
      filenames(ii).(['sampEnt_' 'all' '_' regNames{reg}]) = meanSE; 
  for fi = 2:4
         
        meanSE = mean(sampEnt(idx, :,  grainBounds(fi-1)+1:grainBounds(fi)), 'all' );
        filenames(ii).(['sampEnt_' grainNames{fi-1} '_' regNames{reg}]) = meanSE; 
       
  end
  end

  %% sample entropy asymmetries 
  for asi = 1:length(asymNames)
        asVar = asymVars{asi}; 
        idx1 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{1}); 
        idx2 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{2}); 
        meanSE = mean(sampEnt(idx1, :,  :) - sampEnt(idx2, :,  :), 'all' );
        filenames(ii).(['sampEntAsym_' 'all' '_' asymNames{asi}]) = meanSE; 

  for fi = 2:4 

        meanSE = mean((sampEnt(idx1, :,  grainBounds(fi-1)+1:grainBounds(fi)) - ...
                      sampEnt(idx2, :,  grainBounds(fi-1)+1:grainBounds(fi))) ./ ...
                      (sampEnt(idx1, :,  grainBounds(fi-1)+1:grainBounds(fi)) + ...
                      sampEnt(idx2, :,  grainBounds(fi-1)+1:grainBounds(fi))), 'all' );
        filenames(ii).(['sampEntAsym_' grainNames{fi-1} '_' asymNames{asi}]) = meanSE; 
  end
  end

  %% ISPC connectivity
    ispc = data.ispc; 



  for fi = 2:7  
      %special for the long and short distances
    temp = ispc(freqBounds(fi-1)+1:freqBounds(fi), :,:);
    temp = squeeze(mean(temp, 1));
    temp = triu(temp,1); 
    temp = temp(:); 
    temp(temp==0) = []; 
    shortISPC = mean( temp(distVals<median(distVals)) );
    longISPC = mean( temp(distVals>median(distVals)) );
    filenames(ii).(['ispc_' freqNames{fi-1} '_' 'short']) = shortISPC;
    filenames(ii).(['ispc_' freqNames{fi-1} '_' 'long']) = longISPC;

    temp = ispc(freqBounds(fi-1)+1:freqBounds(fi), :,:);
    temp = squeeze(mean(temp, 1));
  for asi = 1:length(asymNames)
    temp = ispc(freqBounds(fi-1)+1:freqBounds(fi), :,:);
    temp = squeeze(mean(temp, 1));
    asVar = asymVars{asi}; 
    idx1 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{1}); 
    idx2 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{2});
    meanispc = mean(arrayfun(@(x,y)  temp(x, y), idx1, idx2)); 
    filenames(ii).(['ispc_' freqNames{fi-1} '_' asymNames{asi}]) = meanispc;
     


  end




  end


  %% regional PAC
  PAC_cell = data.PACMI; 
  null = data.MINull; 
  phase = data.phasePref; 

  %special indices for low and high frequency for PAC calculations
%   LF_pac = [0, 2, 4, 7, 10];
%   HF_pac = [0, 4, 9, 21];
%   LF_names = {'del', 'the', 'al', 'bet'}; 
%   HF_names = {'bet', 'gL', 'gH'}; 


    %measure 1: PAC modulation index: channel X LF X HF
  PAC = zeros(32, 10, 21); 
  for chan = 1:32
      PAC(chan,:,:) = PAC_cell{chan}; 
  end

    %measure 2: PAC z-scored to their nulls
  PACz = zeros(32,10,21); 
  for chan = 1:32
      chanNull = null{chan}; 
      chanObs = squeeze(PAC(chan,:,:)); 
      nullMean = squeeze(mean(chanNull, 3));
      nullStd = squeeze(std(chanNull, [], 3)); 
      PACz(chan,:,:) = (chanObs - nullMean) ./ nullStd; 
  end

    %measure 3: PAC phase preference
  PACphase = zeros(32,10,21); 
  for chan = 1:32
    PACphase(chan,:,:) = phase{chan};

  end



    for reg = 1:length(regNames)
        regVar = regVars{reg}; 
        idx = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), regVar);
        for LFi = 2:5
        for HFi = 2:4
         if HFi ==2 && LFi==5 %beta to beta, adjust bounds
             HF_pac(1) = 2;
             meanPAC = mean(PAC(idx, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)), 'all' );
             meanPACz = mean(PACz(idx, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)), 'all' );
             tempPhase = PACphase(idx, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi));
             meanPhase= meanangle(tempPhase(:), [0,360]); 
             HF_pac(1) = 0; 
         else
             meanPAC = mean(PAC(idx, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)), 'all' );
             meanPACz = mean(PACz(idx, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)), 'all' );
             tempPhase = PACphase(idx, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi));
             meanPhase= meanangle(tempPhase(:), [0,360]);
         end
        filenames(ii).(['PACmi_' LF_names{LFi-1} '_' HF_names{HFi-1} '_' regNames{reg}]) = meanPAC; 
        filenames(ii).(['PACz_' LF_names{LFi-1} '_' HF_names{HFi-1} '_' regNames{reg}]) = meanPACz; 
        filenames(ii).(['phase_' LF_names{LFi-1} '_' HF_names{HFi-1} '_' regNames{reg}]) = meanPhase; 
        
        end
        end
    end




  %%  asymmetries in PAC

   for asi = 1:length(asymNames)
        asVar = asymVars{asi}; 
        idx1 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{1}); 
        idx2 = arrayfun(@(y) find(cellfun(@(x) strcmp(y, x), labels)), asVar{2}); 
        for LFi = 2:5
        for HFi = 2:4
         if HFi ==2 && LFi==5 %beta to beta, adjust bounds
             HF_pac(1) = 2;
             meanPAC = mean(PAC(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) - ...
                            PAC(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)), 'all' );
             meanPACz = mean(PACz(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) - ...
                             PACz(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)), 'all' );
             tempPhase = PACphase(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) - ...
                         PACphase(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi));
             tempPhase = 180 - abs(abs(tempPhase) - 180);
             meanPhase= mean(tempPhase(:)); 
             HF_pac(1) = 0; 
         else
             meanPAC = mean((PAC(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) - ...
                            PAC(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi))) ./ ...
                            (PAC(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) + ...
                            PAC(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi))), 'all' );
             meanPACz = mean((PACz(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) - ...
                             PACz(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi))) ./ ...
                             (PACz(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) + ...
                             PACz(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi))), 'all' );
             tempPhase = PACphase(idx1, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi)) - ...
                         PACphase(idx2, LF_pac(LFi-1)+1:LF_pac(LFi), HF_pac(HFi-1)+1:HF_pac(HFi));
             tempPhase = 180 - abs(abs(tempPhase) - 180);
             meanPhase= mean(tempPhase(:));
         end
        filenames(ii).(['PACmiAsym_' LF_names{LFi-1} '_' HF_names{HFi-1} '_' asymNames{asi}]) = meanPAC; 
        filenames(ii).(['PACzAsym_' LF_names{LFi-1} '_' HF_names{HFi-1} '_' asymNames{asi}]) = meanPACz; 
        filenames(ii).(['phaseAsym_' LF_names{LFi-1} '_' HF_names{HFi-1} '_' asymNames{asi}]) = meanPhase; 
        
        end
        end
    end



  catch

    errorLog = [errorLog ii]; 
  end




end



%% save out the final struct as a .csv


writetable(struct2table(filenames), 'G:\My Drive\GitHub\SheffieldAutismBiomarkers\autismBiomarkersAllData4.csv')























