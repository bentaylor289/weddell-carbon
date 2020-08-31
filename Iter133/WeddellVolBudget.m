%%% Weddell box Volume Budget 
%%% for the Iter133 2013-2018 volume budget 
%%% looping this one for multiple 5-day periods

clear all
load('../Iter129/grid.mat', 'hFacC', 'hFacW' , 'hFacS', 'DXG', 'DYG', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
sur = volume(:,:,1)/DRF(1);
wkstart = 1;
wkend = 438;
monstart =wkstart; % lazy coding
monend = wkend;

budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
folder = '/local/data/bSOSE/iter133NEW/1day/bsose_i133_2013to2018_1dy_';
folder2 = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5daySnaps_';
Nlat = 60;
lats = ncread('../../../data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 302;
lonE = 30;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcFull = [xc1:2160 1:xc2-1];
tmp = size(xcFull);
FullLen = tmp(2);
WeddellFix = true
ycSW = 277;

n1 = monend-monstart +1; % number of iterations.

advTot = zeros(n1,1);
SSHTot = zeros(n1,1);
FWTot = zeros(n1,1);
Volres = zeros(n1,1);
%bioTot = zeros(n1,1);

for t = monstart:monend

    %t_array = 1;
    t_real = t;
    l = t- monstart+1;
    

%load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'dilut', 'surf', 'tend', 'bio', 'res','corr');


% this will be much harder now! - each one needs three calculations!
% north boundary of the box
adv_tosum = squeeze(ncread(strcat(budgetfolder,'Vvel.nc'), 'VVEL', [1 yc 1 t_real], [Inf 1 Inf 1]));
adv_tosum= double(adv_tosum);
hFac = squeeze(hFacS(xcFull, yc,:)); 
advTot(l) = -sum(sum(adv_tosum(xcFull,:).*hFac.*DXG(xcFull,yc).*repmat(DRF,[1 FullLen])')); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
%mixTot(l) = -sum(sum(mix_tosum(xcFull,:)));

% west boundary
if(WeddellFix)
adv_tosum = squeeze(ncread(strcat(budgetfolder,'Uvel.nc'), 'UVEL', [xc1 ycSW 1 t_real], [1 yc-ycSW Inf 1]));
adv_tosum= double(adv_tosum);
hFac = squeeze(hFacW(xc1, ycSW:yc-1,:));
advTot(l) = advTot(l)+sum(sum(adv_tosum.*hFac.*DYG(xc1, ycSW:yc-1)'.*repmat(DRF,[1 yc-ycSW])')); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 ycSW 1 t_real], [1 yc-ycSW Inf 1]));
%mixTot(l) = mixTot(l)+sum(sum(mix_tosum));
end
%if(WeddellFix ~= true)
%adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
%adv_tosum= double(adv_tosum);
%advTot(l) = advTot(l)+sum(sum(adv_tosum)); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
%mixTot(l) = mixTot(l)+sum(sum(mix_tosum));
%end

% east boundary
adv_tosum = squeeze(ncread(strcat(budgetfolder,'Uvel.nc'), 'UVEL', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
hFac = squeeze(hFacW(xc2, 1:yc-1,:));
advTot(l) = advTot(l)-sum(sum(adv_tosum.*hFac.*DYG(xc2, 1:yc-1)'.*repmat(DRF, [1 yc-1])')); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
%mixTot(l) = mixTot(l)-sum(sum(mix_tosum));

%%

SSHin = ncread(strcat(folder2, 'SSH.nc'), 'ETAN', [1 1 t], [Inf Inf 2] );
SSH = diff(SSHin,1,3)/(3600*120);

FWin = ncread(strcat(folder, 'oceFWflx.nc'), 'oceFWflx', [1 1 5*t-4], [Inf Inf 5]);
FW = mean(FWin, 3)/1035; 

FW2in = ncread(strcat(folder, 'SIatmFW.nc'), 'SIatmFW', [1 1 5*t-4], [Inf Inf 5]);
FW2 = mean(FW2in, 3)/1035; 

SSHTot(l) = sum(sum(double(SSH(xcFull, 1:yc-1)).*sur(xcFull,1:yc-1)));
FWTot(l) = sum(sum(double(FW(xcFull,1:yc-1)).*sur(xcFull,1:yc-1)));
%corrTot(l) = sum(sum(double(corr(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull, 1:yc-1, 1)));
FWTot2(l) = sum(sum(double(FW2(xcFull,1:yc-1)).*sur(xcFull,1:yc-1)));
%tendTot(l) = sum(sum(sum(double(tend(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
%resTot(l) = sum(sum(sum(double(res(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
%bioTot(l) = sum(sum(sum(double(bio(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));

% we need to get that section around the peninsula fixed!
if (WeddellFix)
xcfix = 1788:1812;
ycfix = 1:260;
SSHTot(l) = SSHTot(l) + sum(sum(double(SSH(xcfix, ycfix)).*sur(xcfix,ycfix)));
FWTot(l) = FWTot(l)+ sum(sum(double(FW(xcfix, ycfix)).*sur(xcfix, ycfix)));
FWTot2(l) = FWTot2(l)+ sum(sum(double(FW2(xcfix, ycfix)).*sur(xcfix, ycfix)));

xcfix = 1796:1812;
ycfix = 261:270;

SSHTot(l) = SSHTot(l) + sum(sum(double(SSH(xcfix, ycfix)).*sur(xcfix,ycfix)));
FWTot(l) = FWTot(l)+ sum(sum(double(FW(xcfix, ycfix)).*sur(xcfix, ycfix)));
FWTot2(l) = FWTot2(l)+ sum(sum(double(FW2(xcfix, ycfix)).*sur(xcfix, ycfix)));

xcfix = 1808:1812;
ycfix = 271:277;
SSHTot(l) = SSHTot(l) + sum(sum(double(SSH(xcfix, ycfix)).*sur(xcfix,ycfix)));
FWTot(l) = FWTot(l)+ sum(sum(double(FW(xcfix, ycfix)).*sur(xcfix, ycfix)));
FWTot2(l) = FWTot2(l)+ sum(sum(double(FW2(xcfix, ycfix)).*sur(xcfix, ycfix)));
end
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

Volres(l) = SSHTot(l)-FWTot(l)-advTot(l); 
%clear dilut surf bio tend res corr
end    
%%
Voladv = advTot;
% figure(7)
% x = 1:8;
% barh(x,[advTot mixTot dilutTot surfTot VTrans bioTot tendTot tendTot-dilutTot-surfTot-advTot-VTrans-bioTot])
% 
% set( gca,'yticklabel',{'Advec', 'Mixing','Dilution' ,'Air-sea', 'Volume', 'Bio','Tend','Resid'})
% title('Circumpolar Box, Lat '+string(lat)+', month ' + string(t_real))
%%

%%
save('Wed60Vol'+string(monstart)+'to'+string(monend)+'.mat', 'SSHTot', 'FWTot','FWTot2', 'Voladv', 'Volres')
