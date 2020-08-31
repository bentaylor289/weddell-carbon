%%% 5day CO2 flux output from the Weddell Polynya region 
%%% for the Iter133 2013-2018 5day carbon data

%%% I'll average over a specified length of time, and a specified region

clear all

load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');


wkstart = 1;
wkend = 438;
budgetfolder = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day';

Nlat = 63.2
Slat = 66.5
lats = ncread('/local/data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,ycn] = min(abs(lats+Nlat));
[min2,ycs] = min(abs(lats+Slat));
lonW = 358;
lonE = 8;
xc1 = [6*lonW+1:2160];
xcW = 6*lonW+1;
xc2 = [1:6*lonE];
xcE = 6*lonE;
%xcFull = [xc1:2160 1:xc2-1];

 % number of iterations.
 %bla = size(xcFull);
wbox = squeeze(ncread(strcat(budgetfolder,'_surfCO2flx.nc'), 'BLGCFLX', [xcW ycs wkstart], [Inf ycn-ycs+1 wkend-wkstart+1]));
%exTot = zeros(bla(2), yc-1);
ebox = squeeze(ncread(strcat(budgetfolder,'_surfCO2flx.nc'), 'BLGCFLX', [1 ycs wkstart], [xcE ycn-ycs+1 wkend-wkstart+1]));
    %t_array = 1;
    %t_real = t;
    %l = t- monstart+1;
wflux = sum(sum(wbox.*RAC(xc1, ycs:ycn)));
eflux = sum(sum(ebox.*RAC(xc2, ycs:ycn)));

polyflux = wflux+eflux;


save('polygasex'+string(wkstart)+'to'+string(wkend)+'.mat', 'polyflux')
%%
% plot(budlats, advTot-corrTot','b')
% plot(budlats, tendTot-dilutTot,'r')
% plot(budlats, ones(1,24)*(sum(tendTot)-sum(dilutTot))/24, '--k');
% ylim([-8e6, 4e6])
