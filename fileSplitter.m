function [] = fileSplitter(EEG, path)
    %loop over channels and split the data file for storage all remaining
    %analyses will be carried out on one channel at a time

    disp('tried to split!')
    [folder,name,~] = fileparts(path);
    for ii = 1:32
        curChan = EEG.dataSTD(ii,:,:); 
        srate = EEG.srate; 
        out = struct; 
        out.data = curChan; 
        out.srate = srate; 
        save(join([folder '/' name '_processed_temp' num2str(ii) '.mat'],''), "out")

    end


end