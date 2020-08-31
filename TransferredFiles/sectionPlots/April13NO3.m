%%% Plot DIC along a line of longitude.
% and now a bunch of other stuff! (neutral d and velocities)
adhoc = [145, 190, 210, 270];
for k = 1:4
lon = adhoc(k);% line of longitude, assuming east (0-360)
lat1 = 30; % northern latitude (S)
lat2 = 80; %southern
time = 2; % which month, jan 2008 (1) thru dec 2012 (60)?

density = 1.00; % to convert from volume to weight concentrations
folder = '../../../../data/bSOSE/iter105/monthly/';
DICfile = strcat(folder, 'monDIC.nc');
lats = ncread(DICfile, 'YC');
[min1,yc1] = min(abs(lats+lat1));
[min2,yc2] = min(abs(lats+lat2));

longs = ncread(DICfile, 'XC');
[min3, xc]= min(abs(longs-lon));

% here we find some contours to overlay. 
overlay = [27.1, 27.4, 27.8, 28.04, 28.27];
gammaIn = ncread('../../../../data/bSOSE/iter105/3day/neutrald.nc','GAMMA', [xc yc2 1 floor(10.1*time - 9)], [1 yc1-yc2+1 Inf 10]);
gammaAvg = squeeze(sum(gammaIn, 4))/10;  


    F = ncread(DICfolder+'bsose, 'TRAC01', [xc yc2 1 time ], [1 yc1-yc2+1 Inf 1]);
    F = squeeze(F);
size(gammaIn);
size(gammaAvg);
size(F);
gammaAvg(F < 1) = NaN; 
gammaAvg(F == NaN) = NaN;

% ugly extract contour code .. 
% for each contour get that contour 
% while loop to iterate thru the output
% store the contour x-y data in a 4D array actuall (for all contours for all values). 


%for i = 1:squeeze(size(overlay))
% could use pcolor. 

    H = F;
    H = H';
    H = H./density;
    
% learn how to remove extraneous labels.
% learn how to remove 0's from the deal. 
% figure out if the rho_0 is wrong! 
figure(2)
dep= ncread(DICfile, 'Z', 1 , Inf);
lat = lats(yc2:yc1);
%contour(lat,dep, gammaAvg', overlay)
lev = cat(2, 1.94:0.02:2.2, 2.2:0.01:2.23, 2.23:0.005:2.26)
contourf(lat,dep, H, lev, 'LineColor', 'none')
xlabel('Latitude (S)');
ylabel('Depth');
ax1 = gca;
colorbar(ax1);
ax2= gca;
%contour(ax2, lat, dep, gammaAvg', overlay);
%cmap = colormap(cmap2)
title('DIC (mol/m^3) at' + string(lon)+ ' E at month ' + string(time));
%% nitrate
saveas(gcf, string(lon) + 'DIC.png');

%GamCon(overlay);
hold on
figure(8)
dep= ncread(DICfile, 'Z', 1 , Inf);
lat = lats(yc2:yc1); 
%contour(lat,dep, gammaAvg', overlay)
lev = cat(2, 1.94:0.02:2.2, 2.2:0.01:2.23, 2.23:0.005:2.26)
contourf(lat,dep, H, lev, 'LineColor', 'none')
xlabel('Latitude (S)');
ylabel('Depth');
ax1 = gca;
colorbar(ax1);
ax2= gca;
%contour(ax2, lat, dep, gammaAvg', overlay);
%cmap = colormap(cmap2)
title('DIC (mol/m^3) at '+ string(lon)+ ' E at month ' + string(time));
%saveas(gcf,string(lon) + 'DIC.png')


%% nitrate


%% zonal velocity goes here!
figure(106)

    
    U = ncread(strcat(folder, 'monuvel.nc'), 'UVEL', [xc yc2 1 time ], [1 yc1-yc2+1 Inf 1]);
    U = squeeze(U)';
    ulevels = -0.03:0.005:0.03
    %contourf(lat,dep,U,ulevels)
    xlabel('Latitude (S)');
    ylabel('Depth')
    colorbar();
    title('UVEL (+ means East) along '+ string(lon)+ ' E at month ' + string(time));
    
    %% what we actually need

    figure(7+k)
    load('../grid.mat', 'DYG','DRF', 'hFacW');
    DYG = DYG(xc, yc2:yc1);
    Utrans = DRF.*U.*DYG.*hFacW(xc, yc2:yc1, :);
    size(Utrans);
   Utransd = sum(Utrans);
   for i = 2:1:(yc1-yc2+1)
	   Utransd(i) = Utransd(i) + Utransd(i-1);
   end
   plot(lat, 1e-6*Utransd)
   title('Cumulative Zonal Transport across '+ string(lon)+ ', Sv')
   %saveas(gcf, string(lon)+'Utrans.png');
    %% neutral density
    
    figure(126)
    
    D = ncread('../../../../data/bSOSE/iter105/3day/neutrald.nc', 'GAMMA', [xc yc2 1 time ], [1 (yc1-yc2)+1 Inf 1]);
    D = squeeze(D)';
    dlevels = 26:0.1:29;
    contourf(lat,dep,D,dlevels)
    xlabel('Latitude (S)');
    ylabel('Depth');
    colorbar();
    title('Neutr density along '+ string(lon)+ ' E at month ' + string(time));
  %  saveas(gcf,string(lon) + 'GAMMA.png')
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
   
   end
   function whoCares = GamCon(lays);               
                                                    
            dep = ncread(DICfile, 'Z', 1, Inf);     
	            lat = lats(yc2:yc1);                    
		            whoCares =0 %= get(gca);                
			            contour(lat, dep, gammaAvg, lays)       
			    end                                             
