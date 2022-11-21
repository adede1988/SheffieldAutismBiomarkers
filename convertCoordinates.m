function [EEG] = convertCoordinates(EEG, standardLocs, dataSet)

%input: 
%       EEG             eeglab struct containing one subject's eeg data
%       standardLocs    32 X 3 table containing standard electrode labels,
%                       theta, and phi values. Theta and phi are assumed to
%                       be in degrees. 
%       dataSet         String labeling the dataSet. Some datasets have 
%                       the assumption that the nose is 0 degrees, but 
%                       others are rotated by 90 degrees. Specifying dataSet 
%                       allows for customization of rotation if needed  

%output:
%       EEG             eeglab struct augmented with standardized electrode
%                       locations, interpolated data to the new standard
%                       montage, and Laplacian transformed data

%NOTE: interpolation and Laplacian transform rely on functions
%interpolate_perrinX and laplacian_perrinX which were written by Michael X.
%Cohen

%written by Adam Dede (adam.osman.dede@gmail.com)
%fall 2022


    if ~isfield(EEG, 'dataSTD')
    %% transform coordinates for input channels to a standard circle with radius 1
    
    
    X = [EEG.chanlocs.X]; 
    Y = [EEG.chanlocs.Y]; 
    Z = [EEG.chanlocs.Z]; 
    
    theta = zeros(length(X),1); 
    phi = zeros(length(Y),1); 
    
    %calculate distance between theta and phi based coordinates on unit sphere
    %to the X,Y,Z coordinates given in the montage 
    distFun = @(b,x) sqrt( ( cos(b(2))*sin(b(1)) - x(1) )^2 + ...
                            ( sin(b(2))*sin(b(1)) - x(2) )^2 + ...
                            ( cos(b(1))           - x(3) )^2);
    %use distance function, find theta and phi that minimize distance to input
    %coordinates for each input electrode
    beta0 = [0,1]; 
    for ii = 1:length(X)
        cartCord = [X(ii), Y(ii), Z(ii)];
        tempFun = @(b) distFun(b, cartCord);
        temp = fminsearch(tempFun, beta0);
        theta(ii) = temp(1);
        phi(ii) = temp(2);  

        EEG.chanlocs(ii).Xorig = cartCord(1); 
        EEG.chanlocs(ii).Yorig = cartCord(2); 
        EEG.chanlocs(ii).Zorig = cartCord(3); 

        EEG.chanlocs(ii).X = cos(phi(ii)) .* sin(theta(ii));
        EEG.chanlocs(ii).Y = sin(phi(ii)) .* sin(theta(ii)); 
        EEG.chanlocs(ii).Z = cos(theta(ii)); 
    end

    %% convert coordinates for wanted standard electrode positions to X,Y,Z on unit circle

    %'standard' locations were found to be misaligned. I fixed this with a
    %hack: rotate the standard map by pi/2 for some datasets. This needs to be checked for each dataset

    if strcmp(dataSet, 'femaleASD')
        adjFact = 0; 
    elseif strcmp(dataSet, 'biomarkCon')
        adjFact = 0; 
    else
        adjFact = pi/2; 
    end

    %note assumption here that standardLocs contains theta and phi
    %coordinates in degrees rather than radians
   standardLocs.theta = standardLocs.theta .*pi ./ 180; 
   standardLocs.phi = standardLocs.phi .*pi ./ 180; 
    Xstd = cos(standardLocs.phi-adjFact) .* sin(standardLocs.theta); 
    Ystd = sin(standardLocs.theta) .* sin(standardLocs.phi-adjFact); 
    Zstd = cos(standardLocs.theta); 

    for ii=1:length(Xstd)
        EEG.chanlocsSTD(ii).labels = standardLocs.label(ii); 
        EEG.chanlocsSTD(ii).X = Xstd(ii); 
        EEG.chanlocsSTD(ii).Y = Ystd(ii); 
        EEG.chanlocsSTD(ii).Z = Zstd(ii); 
        EEG.chanlocsSTD(ii).theta = standardLocs.theta(ii);
        EEG.chanlocsSTD(ii).phi = standardLocs.phi(ii);
        
        %top right quadrant
        if Xstd(ii)>0 && Ystd(ii)>0
            EEG.chanlocsSTD(ii).theta = atan(Xstd(ii)/Ystd(ii)); 
        %bottom right quadrant
        elseif Xstd(ii)>0 && Ystd(ii)<0
            EEG.chanlocsSTD(ii).theta = pi/2 + atan(-Ystd(ii) / Xstd(ii)); 
        %bottom left quadrant
        elseif Xstd(ii)<0 && Ystd(ii)<0
            EEG.chanlocsSTD(ii).theta = -pi/2 - atan(-Ystd(ii) / -Xstd(ii)); 
        %top left quadrant
        elseif Xstd(ii)<0 && Ystd(ii)>0
            EEG.chanlocsSTD(ii).theta = -atan(-Xstd(ii)/Ystd(ii)); 
        %along an axis
        elseif Xstd(ii)<0 && Ystd(ii)==0
            EEG.chanlocsSTD(ii).theta = -pi/2; 
        elseif Xstd(ii)>0 && Ystd(ii)==0
            EEG.chanlocsSTD(ii).theta = pi/2;  
        elseif Xstd(ii)==0 && Ystd(ii)>0
            EEG.chanlocsSTD(ii).theta = 0; 
        elseif Xstd(ii)==0 && Ystd(ii)<0
            EEG.chanlocsSTD(ii).theta = pi;  
        end


        cz = [0,0,1]; 
        ct = [Xstd(ii), Ystd(ii), Zstd(ii)]; 
        d = sqrt( (cz(1) - ct(1))^2 + (cz(2) - ct(2))^2 + (cz(3) - ct(3))^2 );
        EEG.chanlocsSTD(ii).radius = 2*asin(d/2); %this value isn't the radius on a sphere, it's the radius for polar plotting on a flat circle

    end

    %% interpolate the data

    allDat = zeros([size(EEG.data,1) + length(Xstd), size(EEG.data,2),size(EEG.data,3)]); 
    allDat(1:EEG.nbchan, :, :) = EEG.data ; 
    allX = [[EEG.chanlocs.X] [EEG.chanlocsSTD.X]];
    allY = [[EEG.chanlocs.Y] [EEG.chanlocsSTD.Y]];
    allZ = [[EEG.chanlocs.Z] [EEG.chanlocsSTD.Z]];
    badElecs = [size(EEG.data,1)+1:size(EEG.data,1)+length(Xstd)];

    test = interpolate_perrinX(allDat,allX,allY,allZ,badElecs); 
    
    datStd = test(badElecs,:,:); 
%     datStd = reshape(datStd, [32, size(EEG.data, [2,3])] ); 

    EEG.dataSTD = datStd; 

    %% apply laplacian on converted data
    EEG.dataLap = laplacian_perrinX(datStd, Xstd, Ystd, Zstd); 

    end
end