%%% Plot DIC along a line of longitude.
% and now a bunch of other stuff! (neutral d and velocities)
figure(76)

lon = 30;% line of longitude, assuming east (0-360)
lat1 = 30; % northern latitude (S)
lat2 = 80; %southern
time = 2; % which month, jan 2008 (1) thru dec 2012 (60)?
DICfolder = '../../../data/bSOSE/iter105/monthly/monDIC.nc';
lats = ncread(DICfolder, 'YC');
[min1,yc1] = min(abs(lats+lat1));
[min2,yc2] = min(abs(lats+lat2));

longs = ncread(DICfolder, 'XC');
[min3, xc]= min(abs(longs-lon));
    
    F = ncread(DICfolder, 'TRAC01', [xc yc2 1 time ], [1 yc1-yc2+1 Inf 1]);
    F = squeeze(F);
    H = F;
    H = H';
    H = H./1.03;
    
% learn how to remove extraneous labels.
% learn how to remove 0's from the deal. 
% figure out if the rho_0 is wrong! 
figure(2)
dep= ncread(DICfolder, 'Z', 1 , Inf);
lat = lats(yc2:yc1) 
lev = cat(2, 1.94:0.02:2.2, 2.2:0.01:2.23, 2.23:0.005:2.26)
contourf(lat,dep, H, lev)
xlabel('Latitude (S)');
ylabel('Depth');
colorbar();
%cmap = colormap(cmap2)
title('DIC at '+ string(lon)+ ' E at month ' + string(time));
saveas(gcf,string(lon) + 'DIC.png')

%% meridional velocity goes here!
figure(106)
    
    %V = ncread('', 'VVEL', [xc yc2 1 time ], [1 yc1-yc2+1 Inf 1]);
    %V = squeeze(V)';
    %vlevels = -0.03:0.005:0.03
    %contourf(lat,dep,V,vlevels)
    xlabel('Latitude (S)');
    ylabel('Depth')
    colorbar();
    title('VVEL (+ means North) along '+ string(lon)+ ' E at month ' + string(time));
   
    %% neutral density
    
    figure(126)
    
    D = ncread('../../../data/bSOSE/iter105/3day/neutrald.nc', 'GAMMA', [xc yc2 1 time ], [1 (yc1-yc2)+1 Inf 1]);
    D = squeeze(D)';
    dlevels = 26:0.1:29
    contourf(lat,dep,D,dlevels)
    xlabel('Latitude (S)');
    ylabel('Depth')
    colorbar();
    title('Neutr density along '+ string(lon)+ ' E at month ' + string(time));
    saveas(gcf,string(lon) + 'GAMMA.png')
    %% vertical velocity
    figure(115)
    
    %W = ncread('monwvel.nc', 'WVEL', [xc yc2 1 time ], [1 yc1-yc2+1 Inf 1]);
    %W = squeeze(W)';
    wlevels = -0.00003:0.000005:0.00003;
    %contourf(lat,dep,W,wlevels)
    xlabel('Latitude (S)');
    ylabel('Depth');
    colorbar();
    title('WVEL along '+ string(lon)+ ' E at month ' + string(time));
