%%% Plot DIC across a line of latitude.
% need to fix this to actually get the proper longitude etc! 
lat = 60;% line of lat, assuming south
lon1 = 300; % starting longitude (0-360, counting East)
lon2 = 360 ;
time = 2; % which month jan 2008 thru dec 2012 (60)?
 vfile = '../../../../data/bSOSE/iter105/monthly/monvvel.nc';

lats= ncread(vfile, 'YG');
[junk, yc] = min(abs(lats+lat));% need the actual function
if (lon1< lon2)
   xc1 = lon1*3;
   xc2 = lon2*3;
   vfile = '../../../../data/bSOSE/iter105/monthly/monvvel.nc';
   V = squeeze(ncread(vfile, 'VVEL', [xc1 yc 1 time], [xc2-xc1+1 1 Inf 1]));
   load('../grid.mat', 'DRF');
   dxG = squeeze(ncread(vfile, 'dxG', [xc1 yc], [xc2-xc1+1 1]));
  hFacS = squeeze(ncread(vfile, 'hFacS', [xc1 yc 1], [xc2-xc1+1 1 Inf]));
    Vtrans = sum(V.*dxG.*DRF'.*hFacS,2);
    for i = 2: (xc2-xc1+1)
	    Vtrans(i ) = Vtrans(i) + Vtrans(i-1);
    end
    figure(5)
    plot([lon1:1/3:lon2], 1e-6*Vtrans)
    title('Meridional Transport across '+string(lat) + 'S, Sv')
    saveas(gcf,string(lat)+'Vtrans2.png')
 end   

