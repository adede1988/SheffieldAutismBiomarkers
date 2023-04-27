function [errorCode] = PACaudit(fileInfo, inputNum, prefix)



    subNum = ceil(inputNum/32); 
    chan = inputNum - 32*floor(inputNum/32); 
    if chan==0
        chan = 32; 
    end

  

        %overall subject summary data
        data = load([fileInfo(subNum).folder '/' fileInfo(subNum).name]).data;
        %filename and folder location
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
        clear data 
        
        
        %frequency params
        frex = logspace(log10(2),log10(80),100);
        transforms = [1,2,3,4]; 
        numfrex = length(frex); 
        stds = linspace(2,5,numfrex);
        
        
        
        
        
        %loop on channels, only loading data for each channel as needed
        % parpool('Processes', 8)
        % for chan = 1:32
        try
            out = load([folderLoc fn{1} '_processed_temp' num2str(chan) '.mat']).out; 
            errorCode = 0;  %everything looks good
            if ~isfield(out, 'MINull')
                errorCode = 1; %missing PAC data
            end
        catch
            errorCode = 2; %file corruption for some reason
        end



end