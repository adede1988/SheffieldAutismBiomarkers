function [allPow] = chanPower(chanDat, frex, numfrex, stds)
    

    allPow = zeros(1, numfrex, size(chanDat.data,3)); 
     
    for snip = 1:size(chanDat.data,3)
        
        hz = linspace(0, chanDat.srate, size(chanDat.data,2)*3 );
        snipDat = chanDat.data(:,:,snip);
        padDat = flip(snipDat,2);
        padDat = [padDat, snipDat, padDat]; 
        fftDat = fft(padDat', [], 1);
        
        for fi = 1:numfrex
            % create Gaussian
            s  = stds(fi)*(2*pi-1)/(4*pi); % normalized width
            x  = hz-frex(fi);                 % shifted frequencies (pre insert the mean for the gaussian)
            fx = exp(-.5*(x/s).^2);    % gaussian
            fx = fx./abs(max(fx));     % gain-normalized

            %filter 
            fDat = ifft( fftDat.*fx')';
            %take the mean power over trials
            allPow(:,fi,snip) = squeeze(mean( abs(2*fDat( : , size(chanDat.data,2)+1:size(chanDat.data,2)*2 ).^2),2 ) );
           
        end
        
    end



end



