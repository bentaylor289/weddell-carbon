%%% Plot DIC across a line of latitude.
% need to fix this to actually get the proper longitude etc! 
lat = 65;% line of lat, assuming south
lon1 = 180; % starting longitude (0-360, counting East)
lon2 = 179 ;
time = 1; % which month jan 2008 thru dec 2012 (60)?

yc= (80-lat)*8; % need the actual function
if (lon1< lon2)
end
    
if(lon1> lon2)
    xc1 = lon1*3;
    xc2 = 1080;
    xcdif1 = xc2-xc1;
    xc3 = 1;
    xc4 = lon2*3;
    xcdif2 = xc4-xc3;
    F = ncread('monDIC.nc', 'TRAC01', [xc1 yc 1 time ], [xcdif1 1 Inf 1]);
    F = squeeze(F);
    G = ncread('monDIC.nc', 'TRAC01', [xc3 yc 1 time ], [xcdif2 1 Inf 1]);
    G = squeeze(G);
    H = cat(1,F,G);
    H = H';
end

%dep = -10:-80:-4090

dep= ncread('monDIC.nc', 'Z', 1 , Inf);
lon = (lon1-360):1/3:lon2-2/3;
lev = 2.20:0.008:2.40;
contourf(lon,dep, H, lev)
xlabel('Longitude (E)');
ylabel('Depth');
colorbar();
title('DIC at '+ string(lat)+ ' S at month ' + string(time));

