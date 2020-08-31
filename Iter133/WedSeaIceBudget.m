%%% Weddell box Sea Ice volume Budget 
%%% for the Iter133 2013-2018 run 
%%% looping this one for multiple 5-day periods

clear all
load('../Iter129/grid.mat', 'hFacC','DXG','DYG', 'RAC', 'DRF');
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

%budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
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

WeddellFix = true
ycSW = 277;

n1 = monend-monstart +1; % number of iterations.

SIexTot = zeros(n1,1);
SImTot = zeros(n1,1);
SIhTot = zeros(n1,1);
SIres = zeros(n1,1);
%bioTot = zeros(n1,1);

for t = monstart:monend

    %t_array = 1;
    t_real = t;
    l = t- monstart+1;
    

%load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'dilut', 'surf', 'tend', 'bio', 'res','corr');


% this will be much harder now! - each one needs three calculations!
% north boundary of the box
adv_tosum = squeeze(ncread(strcat(folder,'SIvheff.nc'), 'SIvheff', [1 yc 5*t_real-4], [Inf 1 5]));
adv_tosum= squeeze(mean(adv_tosum,2));
%hFac = squeeze(hFacS(xcFull, yc,:)); 
advTot(l) = -sum(adv_tosum(xcFull).*DXG(xcFull,yc)); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
%mixTot(l) = -sum(sum(mix_tosum(xcFull,:)));

% west boundary
if(WeddellFix)
adv_tosum = squeeze(ncread(strcat(folder,'SIuheff.nc'), 'SIuheff', [xc1 ycSW  5*t_real-4], [1 yc-ycSW 5]));
adv_tosum= squeeze(mean(adv_tosum,2));
%hFac = squeeze(hFacW(xc1, ycSW:yc-1,:));
advTot(l) = advTot(l)+sum(adv_tosum.*DYG(xc1, ycSW:yc-1)'); % because this is advection to the north!
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
adv_tosum = squeeze(ncread(strcat(folder,'SIuheff.nc'), 'SIuheff', [xc2 1 5*t_real-4], [1 yc-1 5]));
adv_tosum= squeeze(mean(adv_tosum,2));
%hFac = squeeze(hFacW(xc2, 1:yc-1,:));
advTot(l) = advTot(l)-sum(adv_tosum.*DYG(xc2, 1:yc-1)');% because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
%mixTot(l) = mixTot(l)-sum(sum(mix_tosum));
SIexTot(l) = - advTot(l); % advection is an import! 

%%
%if (t==1)
SIhin = ncread(strcat(folder2, 'SeaIceHeff.nc'), 'SIheff', [1 1 t], [Inf yc 2] );
SIh = squeeze(diff(SIhin,1,3))/(3600*120);
%end
%if (t==438)
%SIhin = ncread(strcat(folder, 'SeaIceHeff.nc'), 'SIheff', [1 1 437*5], [Inf yc 6] );
%SIh = (SIhin(:,:,6)- 0.5*SIhin(:,:,2) - 0.5*SIhin(:,:,1))/(3600*108);
%end
%if (t>1 & t<438)
%SIhin = ncread(strcat(folder, 'SeaIceHeff.nc'), 'SIheff', [1 1 5*(t-1)], [Inf yc 7] );
%SIh = (0.5*SIhin(:,:,7)+ 0.5*SIhin(:,:,6) - 0.5*SIhin(:,:,2)-0.5*SIhin(:,:,1))/(3600*120);
%end

FWin = ncread(strcat(folder, 'oceFWflx.nc'), 'oceFWflx', [1 1 5*t-4], [Inf yc 5]);
FW = mean(FWin, 3)/1035; 

FWin2 = ncread(strcat(folder, 'SIatmFW.nc'), 'SIatmFW', [1 1 5*t-4], [Inf yc 5]);
FW2 = mean(FWin2, 3)/1035; 

SIhTot(l) = sum(sum(double(SIh(xcFull, 1:yc-1)).*sur(xcFull,1:yc-1)));
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
SIhTot(l) = SIhTot(l) + sum(sum(double(SIh(xcfix, ycfix)).*sur(xcfix,ycfix)));
FWTot(l) = FWTot(l)+ sum(sum(double(FW(xcfix, ycfix)).*sur(xcfix, ycfix)));
FWTot2(l) = FWTot2(l)+ sum(sum(double(FW2(xcfix, ycfix)).*sur(xcfix, ycfix)));

xcfix = 1796:1812;
ycfix = 261:270;

SIhTot(l) = SIhTot(l) + sum(sum(double(SIh(xcfix, ycfix)).*sur(xcfix,ycfix)));
FWTot(l) = FWTot(l)+ sum(sum(double(FW(xcfix, ycfix)).*sur(xcfix, ycfix)));
FWTot2(l) = FWTot2(l)+ sum(sum(double(FW2(xcfix, ycfix)).*sur(xcfix, ycfix)));

xcfix = 1808:1812;
ycfix = 271:277;
SIhTot(l) = SIhTot(l) + sum(sum(double(SIh(xcfix, ycfix)).*sur(xcfix,ycfix)));
FWTot(l) = FWTot(l)+ sum(sum(double(FW(xcfix, ycfix)).*sur(xcfix, ycfix)));
FWTot2(l) = FWTot2(l)+ sum(sum(double(FW2(xcfix, ycfix)).*sur(xcfix, ycfix)));
end
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 
SImTot(l) = FWTot(l) - FWTot2(l);
SIres(l) = 0.88*SIhTot(l)+SImTot(l)+SIexTot(l); 
%clear dilut surf bio tend res corr
end    
%%
% figure(7)
% x = 1:8;
% barh(x,[advTot mixTot dilutTot surfTot VTrans bioTot tendTot tendTot-dilutTot-surfTot-advTot-VTrans-bioTot])
% 
% set( gca,'yticklabel',{'Advec', 'Mixing','Dilution' ,'Air-sea', 'Volume', 'Bio','Tend','Resid'})
% title('Circumpolar Box, Lat '+string(lat)+', month ' + string(t_real))
%%

%%
save('WeddellSeaIce60fix'+string(monstart)+'to'+string(monend)+'.mat', 'SIhTot', 'FWTot','FWTot2', 'SIexTot', 'SIres')
