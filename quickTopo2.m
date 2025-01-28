function [Zi, Xi, Yi] = quickTopo2(m, Th, Rd, plotit)

  
    Rd = Rd ./ (max(Rd)); %weird to need this adjustment, but there it is

%     sharpness = 0; 
    headrad = 0.70;
    plotrad = .8;
    intrad      = min(1.0,max(Rd)*1.02); 
    GRID_SCALE = 67;        % plot map on a 67X67 grid
    AXHEADFAC = 1.05; 
    BLANKINGRINGWIDTH = .035;% width of the blanking ring
    HEADRINGWIDTH    = .007;% width of the cartoon head ring
    SHADING = 'interp';
    CIRCGRID   = 201;       % number of angles to use in drawing circles
    HEADCOLOR = [0 0 0];    % default head color (black)
    HLINEWIDTH = 1.7;         % default linewidth for head, nose, ears
%     cMap = [linspace(0,1,10000); zeros(10000,1)'; linspace(1,0,10000)]';
    %% convert angle and radius into x and y values

    if max(Th)>2*pi
        Th = pi/180*Th; 
    end
    [y,x] = pol2cart(Th,Rd); 
    colors = m - min(m); 
    colors = round((colors / max(colors))*10000); 
%     colors(colors==0) = 1; 

    

    %% Find plotting channels
    
%     pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
    intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square

    intValues = m(intchans);
    intx  = x(intchans);
    inty  = y(intchans);
    squeezefac = headrad/plotrad;
    intx  = intx*squeezefac;
    inty  = inty*squeezefac;
    xmin = min(-headrad,min(intx)); xmax = max(headrad,max(intx));
    ymin = min(-headrad,min(inty)); ymax = max(headrad,max(inty));
    xi   = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
    yi   = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)
    [Xi,Yi,Zi] = griddata(intx,inty, intValues,xi,yi','v4'); % interpolate data

    mask = (sqrt(Xi.^2 + Yi.^2) <= headrad); % mask outside the plotting circle
    Zi(mask == 0) = NaN;                  % mask non-plotting voxels with NaNs
%     grid = plotrad;                       % unless 'noplot', then 3rd output arg is plotrad
%     delta = xi(2)-xi(1); % length of grid entry

    if plotit>0

%         figure
        set(gca,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC)
        surface(Xi,Yi,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor','interp')
        axis off
        axis equal

        %% Plot filled ring to mask jagged grid boundary
        hwidth = HEADRINGWIDTH;                   % width of head ring
        hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius
        
        if strcmp(SHADING,'interp')
            rwidth = BLANKINGRINGWIDTH*1.3;             % width of blanking outer ring
        else
            rwidth = BLANKINGRINGWIDTH;         % width of blanking outer ring
        end
        rin    =  headrad*(1-rwidth/2);              % inner ring radius
        if hin>rin
            rin = hin;                              % dont blank inside the head ring
        end
        
        circ = linspace(0,2*pi,CIRCGRID);
        rx = sin(circ);
        ry = cos(circ);
        ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
        ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];
        ringh = patch(ringx,ringy,0.01*ones(size(ringx)),get(gcf,'color'),'edgecolor','none','hittest','off'); hold on
        
        %% Plot cartoon head, ears, nose
        headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
        heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
        ringh = patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off'); hold on
        
        % Plot ears and nose
        base  = headrad-.0046;
        basex = 0.18*headrad;                   % nose width
        tip   = 1.15*headrad;
        tiphw = .04*headrad;                    % nose tip half width
        tipr  = .01*headrad;                    % nose tip rounding
        q = .04; % ear lengthening
        EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % headrad = 0.5
        EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
        sf    = headrad/plotrad;                                          % squeeze the model ears and nose
        % by this factor
        plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,2*ones(size([basex;tiphw;0;-tiphw;-basex])),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');                 % plot nose
        plot3(EarX*sf+.18,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')    % plot left ear
        plot3(-EarX*sf - .18,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')   % plot right ear



    end
end