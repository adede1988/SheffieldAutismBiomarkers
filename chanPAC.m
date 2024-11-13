function [rawMeans,MI,pacSig,phasePref, MINull] = chanPAC(chanDat)
    
%input: 
%       chanDat     a struct with the fields data and srate 
%                   data is a channels X time X epochs array
%                   srate is a scalar (e.g. 500 or 1000)

%output: 
%       rawMeans    low frequency X high frequency X phase normalized
%                   average amplitude. Low frequency (LF) values range from 2 to
%                   20 in steps of 2. High frequency (HF) values range from 20
%                   to 100 in steps of 4. Phase is divided into 18 bins.
%                   For every LF - HF pair, phase values sum to 1.0
%       MI          low frequency X high frequency KL divergence values or
%                   modulation index values (same thing, just terminology
%                   difference). KL divergence is calculated relative to a
%                   uniform distribution across the 18 phase bins. Input to
%                   calculate MI is the rawMeans array
%       pacSig      low frequency X high frequency boolean matrix where
%                   true indicates that the modulation index is significantly 
%                   different from chance and false indicates that the 
%                   modulation index is not different from chance.   
%       phasePref   low frequency X high frequency matrix of radian values
%                   indicating the weighted mean phase of the amplitude distribution
%                   for each LF - HF pair. 
%       MINull      200 X low frequency X high frequency array of KL
%                   divergence values calculated for data where the LF and
%                   HF data were first temporally rotated relative to one
%                   another. 

%calculation generally follows phase amplitude coupling analysis presented
%in:
%Pack, Naples, Webb, ..., Levin, ABC-CT Network, preprint.
%Phase-amplitude coupling in autism spectrum disorder: results from the
%autism biomarkers consortium for clinical trials

%with conceptual recommendations drawn from: 
%Tort, Komorowski, Eichenbaum, Kopell. 2010. Measuring phase-amplitude
%coupling between neuronal oscillations of different frequencies. J
%Neurophysiology.

%This function calculates phase amplitude coupling across a range of
%frequencies within a single channel of electrophysiological data. 

%Written by Adam Dede (adam.osman.dede@gmail.com)
%fall, 2022




    lowFrex = [2:2:20]; 
    highFrex = [20:4:100];
    lfSTDs = linspace(2,3.5, length(lowFrex)); 
    hfSTDs = linspace(3.5,6, length(highFrex)); 
    radBounds = [0:2*pi/18:2*pi]; 
    %store the binned, normalized mean amplitude by phase values as
    % trial X LF X HF X phase 
    results = single(zeros(size(chanDat.data,3), length(lowFrex), length(highFrex), 18));
    nullRes = single(zeros([200, size(results)]));

    for snip = 1:size(chanDat.data,3)
        snip
        tic
        hz = linspace(0, chanDat.srate, size(chanDat.data,2)*3 );
        snipDat = chanDat.data(:,:,snip);
        padDat = flip(snipDat,2);
        padDat = [padDat, snipDat, padDat]; 
        fftDat = fft(padDat', [], 1);
            for LF =1:length(lowFrex)
              
                s  = lfSTDs(LF)*(2*pi-1)/(4*pi); % normalized width
                x  = hz-lowFrex(LF);                 % shifted frequencies (pre insert the mean for the gaussian)
                fx = exp(-.5*(x/s).^2);    % gaussian
                fx = fx./abs(max(fx));     % gain-normalized
    
                %filter 
                fDat = ifft( fftDat.*fx').';
                fDat = angle(fDat(:,size(chanDat.data,2)+1:size(chanDat.data,2)*2))+pi; 

                for HF = 1:length(highFrex)
                    s  = hfSTDs(HF)*(2*pi-1)/(4*pi); % normalized width
                    x  = hz-highFrex(HF);                 % shifted frequencies (pre insert the mean for the gaussian)
                    fx = exp(-.5*(x/s).^2);    % gaussian
                    fx = fx./abs(max(fx));     % gain-normalized
        
                    %filter 
                    fDatHF = ifft( fftDat.*fx').';
                    fDatHF = abs(2*fDatHF(:, size(chanDat.data,2)+1:size(chanDat.data,2)*2));
                    
                    %get the mean binned amplitude of the HF signal in
                    %phase bins of the LF signal, organized as 
                    %18 (phase bins) 
                    binnedHFAmp = arrayfun(@(x,y) mean(fDatHF(1, fDat>x & fDat<y)), ...
                       radBounds(1:18), radBounds(2:19)); 
                    binnedHFAmp = binnedHFAmp ./ sum(binnedHFAmp); %normalize by total amp across bins
                    results(snip, LF, HF, :) = single(binnedHFAmp); 

                    for ni = 1:200 %null distribution construction
                        
                        test = circshift(fDatHF, round(size(chanDat.data,2)*rand(1)), 2);
                        binnedHFAmp = arrayfun(@(x,y) mean(test(1, fDat>x & fDat<y)), ...
                            radBounds(1:18), radBounds(2:19));
                        binnedHFAmp = binnedHFAmp ./ sum(binnedHFAmp); %normalize by total amp across bins
                        nullRes(ni,snip, LF, HF, :) = single(binnedHFAmp);
                    end
                   
                end
            end
        toc
    end
            
    rawMeans = squeeze(mean(results, 1)); 
    nullMeans = squeeze(mean(nullRes, 2));

    %convert to KL divergence (modulation index: MI)
    %x = LF
    %y = HF
    %z = chan
    %use arrayfuns to loop across frequencies and channels 
    MI =    reshape(cell2mat(arrayfun(@(y) ...
                                    arrayfun(@(x) sum(rawMeans(x,y,:).*log(rawMeans(x,y,:) ./ (1/18))),...
                                    1:length(lowFrex)),...
                        1:length(highFrex), 'uniformoutput', false)), [10,21]);

    %n = null shuffle add n as another looping variable to represent the
    %null shuffle
    %calculate the KL divergence for all null shuffles
    MINull =  reshape(cell2mat(arrayfun(@(n)  ...
                            reshape(cell2mat(arrayfun(@(y) ...
                                    arrayfun(@(x) sum(nullMeans(n,x,y,:).*log(nullMeans(n,x,y,:) ./ (1/18))),...
                                    1:length(lowFrex)),...
                            1:length(highFrex), 'uniformoutput', false)), [10,21]), ...
                         1:200, 'uniformoutput', false)), [10,21,200]);

    
    %ask whether thet observed value is greater than the 99th percentile. 
    pacSig = MI  > prctile(MINull,99,3); 
   

    %generate representative set of angles to take the mean of. 
    %Basically, creating a distribution of angle observations that is 
    %proportional to the mean values observed, then take the meanangle of it in the next computation.
    %Probs a dumb way to do it, but oh well, it works. 
    phasePref = arrayfun(@(y) ...
                                    arrayfun(@(x)  ...
                                    arrayfun(@(a,b) repmat(b,[a,1]), ...
                                    squeeze(round(1000*rawMeans(x,y,:)))', radBounds(1:18)/pi .* 180 +10 ,'uniformoutput', false),...
                                    1:length(lowFrex), 'uniformoutput', false),...
                        1:length(highFrex), 'uniformoutput', false);

    %take the mean phase angle preference
     phasePref = reshape(cell2mat(arrayfun(@(y) ...
                                    arrayfun(@(x)  ...
                                     meanangle(vertcat(phasePref{x}{y}{:}), [0,360]) ,...
                                    1:length(highFrex)),...
                        1:length(lowFrex), 'uniformoutput', false)), [21,10])';



end




