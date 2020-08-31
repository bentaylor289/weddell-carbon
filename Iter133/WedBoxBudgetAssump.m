%%% Difficult terms of the Weddell carbon budget 
%%% for the Iter133 2013-2018 carbon budget 
%%% testing assumptions and claims made by the tracer_budget document. 
%%%
%%% looping this one for multiple 5-day periods

clear all
load('../Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%
toph = DRF(1);
for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
topA = squeeze(volume(:,:,1)/toph);
wkstart = 1;
wkend = 438;
monstart =wkstart; % lazy coding
monend = wkend;

budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
snaps = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5daySnap';
Nlat = 60;
lats = ncread('../../../data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 302;
lonE = 30;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcFull = [xc1:2160 1:xc2-1];

WeddellFix = true
ycSW = 277;

n1 = monend-monstart +1; % number of iterations.

toptend = zeros(n1,1);
topcont = zeros(n1,1);
partetaC = zeros(n1,1);
partCeta = zeros(n1,1);
bioeta = zeros(n1,1);
adveta = zeros(n1,1);
mixeta = zeros(n1,1);
%tendTot = zeros(n1,1);
%VTrans = zeros(n1,1);
%Cres = zeros(n1,1);

for t = monstart:monend

    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    

load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'tend', 'bio');

tend(isnan(tend)) = 0;
%dilut(isnan(dilut)) = 0;
%surf(isnan(surf)) = 0;
bio(isnan(bio)) = 0;
%res(isnan(res)) = 0;
%corr(isnan(corr)) = 0;
%div(isnan(div)) = 0;

% this will be much harder now! - each one needs three calculations!
% north boundary of the box
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVy_DIC.nc'), 'ADVyTr01', [1 yc 1 t_real], [Inf 1 1 1]));
adv_tosum= double(adv_tosum);
myEta = squeeze(ncread(strcat(budgetfolder, 'SSH.nc'), 'ETAN', [1 yc-1 t_real], [Inf 2 1]));
myEta = mean(myEta,2);
adveta(l) = -sum(sum(adv_tosum(xcFull).*myEta(xcFull))); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 1 1]));
mixeta(l) = -sum(sum(mix_tosum(xcFull).*myEta(xcFull)));

% west boundary
if(WeddellFix)
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 ycSW 1 t_real], [1 yc-ycSW 1 1]));
adv_tosum= double(adv_tosum);
myEta = squeeze(ncread(strcat(budgetfolder, 'SSH.nc'), 'ETAN', [xc1-1 ycSW t_real], [2 yc-ycSW 1]));
myEta = squeeze(mean(myEta,1));
adveta(l) = adveta(l)+sum(sum(adv_tosum.*myEta)); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 ycSW 1 t_real], [1 yc-ycSW 1 1]));
mixeta(l) = mixeta(l)+sum(sum(mix_tosum.*myEta));
end
if(WeddellFix ~= true)% not changed yet!
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = advTot(l)+sum(sum(adv_tosum)); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
mixTot(l) = mixTot(l)+sum(sum(mix_tosum));
end

% east boundary
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc2 1 1 t_real], [1 yc-1 1 1]));
adv_tosum= double(adv_tosum);
myEta = squeeze(ncread(strcat(budgetfolder, 'SSH.nc'), 'ETAN', [xc2-1 1 t_real],[2 yc-1 1]));
myEta = squeeze(mean(myEta,1));
adveta(l) = adveta(l)-sum(sum(adv_tosum.*myEta)); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc2 1 1 t_real], [1 yc-1 1 1]));
mixeta(l) = mixeta(l)-sum(sum(mix_tosum.*myEta));

adveta(l) = adveta(l)/toph;
mixeta(l) = mixeta(l)/toph;
%%
myEta = squeeze(ncread(strcat(budgetfolder, 'SSH.nc'), 'ETAN', [1 1 t_real], [Inf yc-1 1]));
etaSn = squeeze(ncread(strcat(snaps, 's_SSH.nc'), 'ETAN', [1 1 t_real], [Inf yc-1 2])); 
DICSn = squeeze(ncread(strcat(snaps, 'Shots_DIC.nc'), 'TRAC01', [1 1 1 t_real],[Inf yc-1 1 2]));

topcont(l) = 1/(120*3600)*sum(sum(squeeze((toph+ etaSn(xcFull, 1:yc-1, 2)).*topA(xcFull,1:yc-1).*DICSn(xcFull,1:yc-1,2) - (toph+etaSn(xcFull, 1:yc-1,1)).*topA(xcFull,1:yc-1).*DICSn(xcFull,1:yc-1,1))));
partCeta(l) = sum(sum(double(tend(xcFull,1:yc-1, 1, t_array)).*topA(xcFull,1:yc-1,1).*etaSn(xcFull, 1:yc-1, 2)));
partetaC(l) = 1/(120*3600)* sum(sum(squeeze(diff(etaSn(xcFull, :,:),1,3)).*topA(xcFull, 1:yc-1, 1).*DICSn(xcFull,:,1)));
toptend(l) = sum(sum(double(tend(xcFull,1:yc-1,1, t_array)).*volume(xcFull,1:yc-1,1),2),1);
%resTot = sum(sum(sum(double(res(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
bioeta(l) = sum(sum(double(bio(xcFull,1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1).*myEta(xcFull,:)/toph,2),1);

% we need to get that section around the peninsula fixed!
if (WeddellFix)
xcfix = 1788:1812;
ycfix = 1:260;
topcont(l) = topcont(l) + 1/(120*3600)*sum(sum(squeeze((toph+ etaSn(xcfix, ycfix, 2)).*topA(xcfix,ycfix).*DICSn(xcfix,ycfix,2) - (toph+etaSn(xcfix, ycfix,1)).*topA(xcfix,ycfix).*DICSn(xcfix,ycfix,1))));
partCeta(l) = partCeta(l) + sum(sum(double(tend(xcfix, ycfix, 1, t_array)).*topA(xcfix, ycfix,1).*etaSn(xcfix, ycfix, 2)));
partetaC(l) = partetaC(l) + 1/(120*3600)* sum(sum(squeeze(diff(etaSn(xcfix,ycfix,:),1,3)).*topA(xcfix,ycfix, 1).*DICSn(xcfix,ycfix,1)));
toptend(l) = toptend(l) + sum(sum(double(tend(xcfix,ycfix,1, t_array)).*volume(xcfix, ycfix,1),2),1);
bioeta(l) = bioeta(l)+sum(sum(double(bio(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1).*myEta(xcfix, ycfix),2),1);

xcfix = 1796:1812;
ycfix = 261:270;

topcont(l) = topcont(l) + 1/(120*3600)*sum(sum(squeeze((toph+ etaSn(xcfix, ycfix, 2)).*topA(xcfix,ycfix).*DICSn(xcfix,ycfix,2) - (toph+etaSn(xcfix, ycfix,1)).*topA(xcfix,ycfix).*DICSn(xcfix,ycfix,1))));
partCeta(l) = partCeta(l) + sum(sum(double(tend(xcfix, ycfix, 1, t_array)).*topA(xcfix, ycfix,1).*etaSn(xcfix, ycfix, 2)));
partetaC(l) = partetaC(l) + 1/(120*3600)* sum(sum(squeeze(diff(etaSn(xcfix,ycfix,:),1,3)).*topA(xcfix,ycfix, 1).*DICSn(xcfix,ycfix,1)));
toptend(l) = toptend(l) + sum(sum(double(tend(xcfix,ycfix,1, t_array)).*volume(xcfix, ycfix,1),2),1);
bioeta(l) = bioeta(l)+sum(sum(double(bio(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1).*myEta(xcfix,ycfix),2),1);

xcfix = 1808:1812;
ycfix = 271:277;

topcont(l) = topcont(l) + 1/(120*3600)*sum(sum(squeeze((toph+ etaSn(xcfix, ycfix, 2)).*topA(xcfix,ycfix).*DICSn(xcfix,ycfix,2) - (toph+etaSn(xcfix, ycfix,1)).*topA(xcfix,ycfix).*DICSn(xcfix,ycfix,1))));
partCeta(l) = partCeta(l) + sum(sum(double(tend(xcfix, ycfix, 1, t_array)).*topA(xcfix, ycfix,1).*etaSn(xcfix, ycfix, 2)));
partetaC(l) = partetaC(l) + 1/(120*3600)* sum(sum(squeeze(diff(etaSn(xcfix,ycfix,:),1,3)).*topA(xcfix,ycfix, 1).*DICSn(xcfix,ycfix,1)));
toptend(l) = toptend(l) + sum(sum(double(tend(xcfix,ycfix,1, t_array)).*volume(xcfix, ycfix,1),2),1);
bioeta(l) = bioeta(l)+ sum(sum(double(bio(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1).*myEta(xcfix,ycfix),2),1);
end
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

%Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 
clear bio tend
end    
%%


%%
save('Weddell'+string(Nlat)+'Wk'+string(monstart)+'to'+string(monend)+'ASSUMP.mat', 'adveta', 'mixeta', 'partetaC','toptend', 'bioeta', 'partCeta', 'topcont');
