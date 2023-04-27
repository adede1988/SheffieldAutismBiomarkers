function [ispc, EEG] = getISPC(EEG, frex, numfrex, stds, path)

    if isfield(EEG, 'ispc')
        ispc = EEG.ispc; 
    else
    
        ispc = zeros([numfrex, size(EEG.dataSTD,1), size(EEG.dataSTD,1)]);
        for fi = 1:numfrex
            tempISPC = zeros([size(EEG.dataSTD,1), size(EEG.dataSTD,1), size(EEG.dataSTD,3)]);
        
            for snip = 1:size(EEG.dataSTD,3)
                hz = linspace(0, EEG.srate, size(EEG.dataSTD,2)*3 );
                snipDat = EEG.dataLap(:,:,snip);
                padDat = flip(snipDat,2);
                padDat = [padDat, snipDat, padDat]; 
                fftDat = fft(padDat', [], 1);
        
                % create Gaussian
                s  = stds(fi)*(2*pi-1)/(4*pi); % normalized width
                x  = hz-frex(fi);                 % shifted frequencies (pre insert the mean for the gaussian)
                fx = exp(-.5*(x/s).^2);    % gaussian
                fx = fx./abs(max(fx));     % gain-normalized
        
                %filter and turn into phase time series
                tempFiltDat = angle(2*( ifft( fftDat.*fx') )') ;
                tempFiltDat = tempFiltDat(:,size(EEG.dataSTD,2)+1: size(EEG.dataSTD,2)*2); 
                tempISPC(:,:,snip) =  cell2mat(arrayfun(@(x) ...
                        abs(sum(exp(1i* bsxfun(@minus, tempFiltDat, tempFiltDat(x,:))),2)./size(EEG.dataSTD,2) ), ...
                        [1:size(snipDat,1)], 'UniformOutput', false)) ;
        
        
        
            end
            
            ispc(fi,:,:) = squeeze(mean(tempISPC,3)); 
        
        end
        
        EEG.ispc = ispc; 
        
        
        [folder,name,~] = fileparts(path);
        save(join([folder '/' name '_processed.mat'],''), "EEG")
    
    
    
    end




end














