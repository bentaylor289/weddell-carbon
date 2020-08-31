%%% Averaging T S and Neutral D over a length of time in the Weddell 

%%% for the Iter133 2013-2018 5day carbon data

%%% I'll average over a specified length of time, and a specified region

clear all

load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%
depth = zeros(1,52);
for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
 depth(k) = sum(DRF(1:k)) - DRF(k)/2;
end

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

Gfirst = zeros(bla(2), yc-1, 3);
Tfirst = zeros(bla(2), yc-1, 3);
Sfirst = zeros(bla(2), yc-1, 3);
cells = [13 20 25]; % 100 200 300 m approx
hF = hFacC(xcF,1:yc-1,:);
hF(hF>0)=1;
% all I'm doing is just taking the mean for 
    
for k= 1:3
% this will be much harder now! - each one needs three calculations!
% north boundary of the box
Gtmp = squeeze(ncread(strcat(budgetfolder,'GAMMA.nc'), 'gamma', [1 1 cells(k) wkstart], [Inf yc-1 1 wkend-wkstart+1]));
Stmp = squeeze(ncread(strcat(budgetfolder,'Salt.nc'), 'SALT', [1 1 cells(k) wkstart], [Inf yc-1 1 wkend-wkstart+1]));
Ttmp = squeeze(ncread(strcat(budgetfolder,'Theta.nc'), 'THETA', [1 1 cells(k) wkstart], [Inf yc-1 1 wkend-wkstart+1]));
Gfirst(:,:,k) = mean(Gtmp(xcF, :,:),3).*hF(:,:,cells(k));
Sfirst(:,:,k) = mean(Stmp(xcF, :,:),3).*hF(:,:,cells(k));
Tfirst(:,:,k) = mean(Ttmp(xcF, :,:),3).*hF(:,:,cells(k));
end


wkstart = 366;
wkend = 438;
%bla = size(xcF);
%load('/local/projects/bSOSE_carbon_Ben/Iter133/BlockyLayerBudget1to438.mat', 'isINS'); 

Glast = zeros(bla(2), yc-1, 3);
Tlast = zeros(bla(2), yc-1, 3);
Slast = zeros(bla(2), yc-1, 3);
cells = [13 20 25]; % 100 200 300 m approx
%hF = hFacC(xcF,1:yc-1,:);
%hF(hF>0)=1;
% all I'm doing is just taking the mean for 
    
for k= 1:3
Gtmp = squeeze(ncread(strcat(budgetfolder,'GAMMA.nc'), 'gamma', [1 1 cells(k) wkstart], [Inf yc-1 1 wkend-wkstart+1]));
Stmp = squeeze(ncread(strcat(budgetfolder,'Salt.nc'), 'SALT', [1 1 cells(k) wkstart], [Inf yc-1 1 wkend-wkstart+1]));
Ttmp = squeeze(ncread(strcat(budgetfolder,'Theta.nc'), 'THETA', [1 1 cells(k) wkstart], [Inf yc-1 1 wkend-wkstart+1]));
Glast(:,:,k) = mean(Gtmp(xcF, :,:),3).*hF(:,:,cells(k));
Slast(:,:,k) = mean(Stmp(xcF, :,:),3).*hF(:,:,cells(k));
Tlast(:,:,k) = mean(Ttmp(xcF, :,:),3).*hF(:,:,cells(k));
end


save('TSGamLastvsFirst.mat','cells', 'Tfirst', 'Tlast', 'Slast', 'Sfirst','Glast','Gfirst')
%%
