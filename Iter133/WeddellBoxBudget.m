%%% Weddell box carbon budget 
%%% for the Iter133 2013-2018 carbon budget 
%%% looping this one for multiple 5-day periods

clear all
load('../Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
wkstart = 1;
wkend = 438;
monstart =wkstart; % lazy coding
monend = wkend;

budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';

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

advTot = zeros(n1,1);
mixTot = zeros(n1,1);
resTot = zeros(n1,1);
surfTot = zeros(n1,1);
bioTot = zeros(n1,1);
dilutTot = zeros(n1,1);
corrTot = zeros(n1,1);
tendTot = zeros(n1,1);
VTrans = zeros(n1,1);
Cres = zeros(n1,1);

for t = monstart:monend

    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    

load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'dilut', 'surf', 'tend', 'bio', 'res','corr');

tend(isnan(tend)) = 0;
dilut(isnan(dilut)) = 0;
surf(isnan(surf)) = 0;
bio(isnan(bio)) = 0;
res(isnan(res)) = 0;
corr(isnan(corr)) = 0;
%div(isnan(div)) = 0;

% this will be much harder now! - each one needs three calculations!
% north boundary of the box
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVy_DIC.nc'), 'ADVyTr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = -sum(sum(adv_tosum(xcFull,:))); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
mixTot(l) = -sum(sum(mix_tosum(xcFull,:)));

% west boundary
if(WeddellFix)
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 ycSW 1 t_real], [1 yc-ycSW Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = advTot(l)+sum(sum(adv_tosum)); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 ycSW 1 t_real], [1 yc-ycSW Inf 1]));
mixTot(l) = mixTot(l)+sum(sum(mix_tosum));
end
if(WeddellFix ~= true)
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = advTot(l)+sum(sum(adv_tosum)); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
mixTot(l) = mixTot(l)+sum(sum(mix_tosum));
end

% east boundary
adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = advTot(l)-sum(sum(adv_tosum)); % because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
mixTot(l) = mixTot(l)-sum(sum(mix_tosum));

%%
dilutTot(l) = sum(sum(double(dilut(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1)));
surfTot(l) = sum(sum(double(surf(xcFull,1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1)));
corrTot(l) = sum(sum(double(corr(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull, 1:yc-1, 1)));
tendTot(l) = sum(sum(sum(double(tend(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
resTot(l) = sum(sum(sum(double(res(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
bioTot(l) = sum(sum(sum(double(bio(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));

% we need to get that section around the peninsula fixed!
if (WeddellFix)
xcfix = 1788:1812;
ycfix = 1:260;
dilutTot(l) = dilutTot(l) + sum(sum(double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1)));
surfTot(l) = surfTot(l)+ sum(sum(double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1)));
corrTot(l) = corrTot(l)+ sum(sum(double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1)));
tendTot(l) = tendTot(l) + sum(sum(sum(double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
resTot(l) = resTot(l)+sum(sum(sum(double(res(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
bioTot(l) = bioTot(l)+sum(sum(sum(double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));

xcfix = 1796:1812;
ycfix = 261:270;

dilutTot(l) = dilutTot(l) + sum(sum(double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1)));
surfTot(l) = surfTot(l)+ sum(sum(double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1)));
corrTot(l) = corrTot(l)+ sum(sum(double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1)));
tendTot(l) = tendTot(l) + sum(sum(sum(double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
resTot(l) = resTot(l)+sum(sum(sum(double(res(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
bioTot(l) = bioTot(l)+sum(sum(sum(double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));

xcfix = 1808:1812;
ycfix = 271:277;
dilutTot(l) = dilutTot(l) + sum(sum(double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1)));
surfTot(l) = surfTot(l)+ sum(sum(double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1)));
corrTot(l) = corrTot(l)+ sum(sum(double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1)));
tendTot(l) = tendTot(l) + sum(sum(sum(double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
resTot(l) = resTot(l)+sum(sum(sum(double(res(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
bioTot(l) = bioTot(l)+sum(sum(sum(double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
end
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 
clear dilut surf bio tend res corr
end    
%%
% figure(7)
% x = 1:8;
% barh(x,[advTot mixTot dilutTot surfTot VTrans bioTot tendTot tendTot-dilutTot-surfTot-advTot-VTrans-bioTot])
% 
% set( gca,'yticklabel',{'Advec', 'Mixing','Dilution' ,'Air-sea', 'Volume', 'Bio','Tend','Resid'})
% title('Circumpolar Box, Lat '+string(lat)+', month ' + string(t_real))
%%


cv = 3.7843e-4; % units of Tg/yr
%Cres = Cres * cv;
tendTot = tendTot * cv; 
%resTot = resTot *cv;
advTot=  advTot*cv;
surfTot = surfTot*cv;
bioTot = bioTot*cv;
corrTot = corrTot*cv;
mixTot = mixTot*cv;
resTot = (resTot+Cres)*cv;
dilutTot = dilutTot*cv;
%corrTot= corrTot';
budlats = monstart:monend;% sloppy name
close all
figure(61)
%plot(budlats, 100*Cres)
%caxis([-200, 200])
hold on
%plot(budlats, tendTot);
plot(budlats, tendTot+corrTot-dilutTot)
%plot(budlats, dilutTot, 'y')
plot(budlats, advTot)
plot(budlats, surfTot)
plot(budlats, bioTot)
plot(budlats, 100*(resTot))
%%
legend('Storage Tend', 'Advection', 'AirSea','Bio','100*Resid')
%%
%legend('10^3*Resid', 'Tend', 'Modif Tend', 'Dilution', 'Advec', 'VolCorr', 'AirSea','Bio')

title('Weddell DIC budget terms')%+string(lat))
ylabel('Tg/yr')
xlabel('5-day periods')


saveas(gcf, 'weddelltest2.png')
%%
% plot(budlats, advTot-corrTot','b')
% plot(budlats, tendTot-dilutTot,'r')
% plot(budlats, ones(1,24)*(sum(tendTot)-sum(dilutTot))/24, '--k');
% ylim([-8e6, 4e6])
%%
%legend('10^3*Resid','Surf', 'Bio','Cor+Adv', 'Tend', 'Tend Avg');

%%



%%
save('Weddell'+string(Nlat)+'budgetWk'+string(monstart)+'to'+string(monend)+'WEDFIX.mat', 'tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot','resTot');
