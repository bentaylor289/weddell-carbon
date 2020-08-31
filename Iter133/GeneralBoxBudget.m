%%% General box carbon budget;  
%%% for the Iter133 2013-2018 carbon budget 
%%% looping this one for multiple 5-day periods

clear all
%% inputs to change

wkstart = 1;
wkend = 1;

budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
% also dependencies on cellwise budget folder, grid file. 
Nlat = 60;
lonW = 160;
lonE = 230; 

WeddellFix = false; % fix for Antarctic peninsula (we want to stop the Western boundary before 78 S).
% Weddell fix also includes hardcoded locations for Weddell. 
ycSW = 277; % only used if Weddell fix is on

%% Code begins here

load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat', 'hFacC', 'RAC', 'DRF','YC');
lats = YC(1,:);
[min1,yc] = min(abs(lats+Nlat));

xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
if (lonW>lonE)
xcFull = [xc1:2160 1:xc2-1];
end
if (lonE> lonW)
	xcFull = [xc1: xc2-1];
end


monstart =wkstart; % lazy coding
monend = wkend;
n1 = monend-monstart +1; % number of iterations.


volume = zeros(size(hFacC));

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end


advTot = zeros(n1,1);
mixTot = zeros(n1,1);
resCell = zeros(n1,1);%residual from the cellwise budget, which should be tiny!
surfTot = zeros(n1,1);
bioTot = zeros(n1,1);
dilutTot = zeros(n1,1);
corrTot = zeros(n1,1);
tendTot = zeros(n1,1);
resTot = zeros(n1,1);

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
resCell(l) = sum(sum(sum(double(res(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
bioTot(l) = sum(sum(sum(double(bio(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));

% we need to get that section around the peninsula fixed!
if (WeddellFix)
xcfix = 1788:1812;
ycfix = 1:260;
dilutTot(l) = dilutTot(l) + sum(sum(double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1)));
surfTot(l) = surfTot(l)+ sum(sum(double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1)));
corrTot(l) = corrTot(l)+ sum(sum(double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1)));
tendTot(l) = tendTot(l) + sum(sum(sum(double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
resCell(l) = resCell(l)+sum(sum(sum(double(res(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
bioTot(l) = bioTot(l)+sum(sum(sum(double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));

xcfix = 1796:1812;
ycfix = 261:270;

dilutTot(l) = dilutTot(l) + sum(sum(double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1)));
surfTot(l) = surfTot(l)+ sum(sum(double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1)));
corrTot(l) = corrTot(l)+ sum(sum(double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1)));
tendTot(l) = tendTot(l) + sum(sum(sum(double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
resCell(l) = resCell(l)+sum(sum(sum(double(res(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
bioTot(l) = bioTot(l)+sum(sum(sum(double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));

xcfix = 1808:1812;
ycfix = 271:277;
dilutTot(l) = dilutTot(l) + sum(sum(double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1)));
surfTot(l) = surfTot(l)+ sum(sum(double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1)));
corrTot(l) = corrTot(l)+ sum(sum(double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1)));
tendTot(l) = tendTot(l) + sum(sum(sum(double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
resCell(l) = resCell(l)+sum(sum(sum(double(res(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
bioTot(l) = bioTot(l)+sum(sum(sum(double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:))));
end
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

resTot(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l) - mixTot(l)+ corrTot(l); 
clear dilut surf bio tend res corr
end    
%%

cv = 3.7843e-4; % units of Tg/yr
tendTot = tendTot * cv; 
resTot = resTot *cv;
advTot=  advTot*cv;
surfTot = surfTot*cv;
bioTot = bioTot*cv;
corrTot = corrTot*cv;
mixTot = mixTot*cv;
resCell = resCell*cv;
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

title('Box DIC budget terms')%+string(lat))
ylabel('Tg/yr')
xlabel('5-day periods')


saveas(gcf, 'weddelltest2.png')



%%
save('BoxbudgetWk'+string(monstart)+'to'+string(monend)+'.mat', 'Nlat', 'lonE','lonW','resCell', 'tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot','resTot');
