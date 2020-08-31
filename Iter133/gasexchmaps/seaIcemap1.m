%%% Averaging depth of a density layer over a length of time in the Weddell 
%%% for the Iter133 2013-2018 5day carbon data

%%% I'll average over a specified length of time, and a specified region

clear all

load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%
depth = zeros(1,52);
for k=1:52
 %volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
 depth = sum(DRF(1:k)) - DRF(k)/2;
end



z1 = 20; %190 m
z2 = 33; %750 m
zd = z2-z1+1;
budgetfolder = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';

Nlat = 53;
lats = ncread('/local/data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 290;
lonE = 35;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcF = [xc1:2160 1:xc2-1];


% number of iterations.

wkstart = 1;
wkend = 73;
bla = size(xcF);
load('/local/projects/bSOSE_carbon_Ben/Iter133/BlockyLayerBudget1to438.mat', 'isINS'); 

SIfirstM = zeros(bla(2), yc-1);
% all I'm doing is just taking the mean for 
for t = wkstart:wkend
    

% this will be much harder now! - each one needs three calculations!
% north boundary of the box
SItmp = squeeze(ncread(strcat(budgetfolder,'SeaIceArea.nc'), 'SIarea', [1 1 t], [Inf yc-1 1]));
SItmp = SItmp(xcF,1:yc-1);
SIfirstM = SIfirstM + SItmp;
end    

SIfirstM = SIfirstM/(wkend-wkstart+1);
 
wkstart= 5*73+1; wkend=  438;
SIlastM = zeros(bla(2), yc-1);
% all I'm doing is just taking the mean for 
for t = wkstart:wkend
SItmp = squeeze(ncread(strcat(budgetfolder,'SeaIceArea.nc'), 'SIarea', [1 1 t], [Inf yc-1 1]));
SItmp = SItmp(xcF,1:yc-1);
SIlastM = SIlastM + SItmp;
end    

SIlastM = SIlastM/(wkend-wkstart+1);
 
save('SImean.mat','SIlastM', 'SIfirstM')
%%
