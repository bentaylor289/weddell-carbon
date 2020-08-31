%%% Weddell layerwise carbon budget 
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

densityBounds = [23 27.7 28 28.1 28.2 28.27 28.31 28.35 28.37 28.4 29];
%densityBounds = [ 23 27.55 27.7 27.85 28 28.1 28.2 28.27 28.35 28.4 29]; 
mncurref= densityBounds;
lays = length(mncurref)-1;

advTot = zeros(n1,lays);
mixTot = zeros(n1,lays);
resTot = zeros(n1,lays);
surfTot = zeros(n1,lays);
bioTot = zeros(n1,lays);
dilutTot = zeros(n1,lays);
corrTot = zeros(n1,lays);
tendTot = zeros(n1,lays);
VTrans = zeros(n1,lays);
Cres = zeros(n1,lays);
LayVol = zeros(n1, lays); 

for t = monstart:monend

    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    

load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'mix','dilut', 'surf', 'tend', 'bio', 'res','corr');

tend(isnan(tend)) = 0;
mix(isnan(mix)) = 0;
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
% we need 

G = squeeze(ncread(strcat(budgetfolder, 'GAMMA.nc'), 'gamma', [1 1 1 t_real], [Inf yc-1 Inf 1]));
G1 = G(xcFull,:,:);
binner = discretize(G1, mncurref);
dilutTMP = double(dilut(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1);
surfTMP = double(surf(xcFull,1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1);
corrTMP = double(corr(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull, 1:yc-1, 1);
tendTMP = double(tend(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:);
mixTMP = double(mix(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:);
bioTMP = double(bio(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:);
volTMP = volume(xcFull, 1:yc-1, :);
for i = 1:lays
dilutTot(l,i) = sum(sum(dilutTMP(binner(:,:,1)==i)));
surfTot(l,i) = sum(sum(surfTMP(binner(:,:,1)==i)));
corrTot(l,i) = sum(sum(corrTMP(binner(:,:,1)==i)));
tendTot(l,i) = sum(sum(sum(tendTMP(binner==i))));
bioTot(l,i) = sum(sum(sum(bioTMP(binner==i))));
mixTot(l,i) = sum(sum(sum(mixTMP(binner==i))));
LayVol(l,i) = sum(sum(sum(volTMP(binner==i))));
end

% we need to get that section around the peninsula fixed!
if (WeddellFix)
xcfix = 1788:1812;
ycfix = 1:260;

Gfix = G(xcfix, ycfix,:); 
binner = discretize(Gfix,mncurref);
dilutTMP = double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1);
surfTMP = double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1);
corrTMP = double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1);
tendTMP = double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
mixTMP = double(mix(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
bioTMP = double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
volTMP = volume(xcfix, ycfix,:);
for i = 1:lays
dilutTot(l,i) = dilutTot(l,i)+ sum(sum(dilutTMP(binner(:,:,1)==i)));
surfTot(l,i) = surfTot(l,i)+ sum(sum(surfTMP(binner(:,:,1)==i)));
corrTot(l,i) = corrTot(l,i) + sum(sum(corrTMP(binner(:,:,1)==i)));
tendTot(l,i) = tendTot(l,i)+ sum(sum(sum(tendTMP(binner==i))));
bioTot(l,i) = bioTot(l,i)+ sum(sum(sum(bioTMP(binner==i))));
mixTot(l,i) = mixTot(l,i)+ sum(sum(sum(mixTMP(binner==i))));
LayVol(l,i) = LayVol(l,i)+ sum(sum(sum(volTMP(binner==i))));
end

xcfix = 1796:1812;
ycfix = 261:270;

Gfix = G(xcfix, ycfix,:); 
binner = discretize(Gfix,mncurref);
dilutTMP = double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1);
surfTMP = double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1);
corrTMP = double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1);
tendTMP = double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
mixTMP = double(mix(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
bioTMP = double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
volTMP = volume(xcfix, ycfix,:);

for i = 1:lays
dilutTot(l,i) = dilutTot(l,i)+ sum(sum(dilutTMP(binner(:,:,1)==i)));
surfTot(l,i) = surfTot(l,i)+ sum(sum(surfTMP(binner(:,:,1)==i)));
corrTot(l,i) = corrTot(l,i) + sum(sum(corrTMP(binner(:,:,1)==i)));
tendTot(l,i) = tendTot(l,i)+ sum(sum(sum(tendTMP(binner==i))));
bioTot(l,i) = bioTot(l,i)+ sum(sum(sum(bioTMP(binner==i))));
mixTot(l,i) = mixTot(l,i)+ sum(sum(sum(mixTMP(binner==i))));
LayVol(l,i) = LayVol(l,i)+ sum(sum(sum(volTMP(binner==i))));
end

xcfix = 1808:1812;
ycfix = 271:277;

Gfix = G(xcfix, ycfix,:); 
binner = discretize(Gfix,mncurref);
dilutTMP = double(dilut(xcfix, ycfix, 1, t_array)).*volume(xcfix,ycfix,1);
surfTMP = double(surf(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix,1);
corrTMP = double(corr(xcfix, ycfix, 1, t_array)).*volume(xcfix, ycfix, 1);
tendTMP = double(tend(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
mixTMP = double(mix(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
bioTMP = double(bio(xcfix, ycfix, :, t_array)).*volume(xcfix, ycfix,:);
volTMP = volume(xcfix, ycfix,:);

for i = 1:lays
dilutTot(l,i) = dilutTot(l,i)+ sum(sum(dilutTMP(binner(:,:,1)==i)));
surfTot(l,i) = surfTot(l,i)+ sum(sum(surfTMP(binner(:,:,1)==i)));
corrTot(l,i) = corrTot(l,i) + sum(sum(corrTMP(binner(:,:,1)==i)));
tendTot(l,i) = tendTot(l,i)+ sum(sum(sum(tendTMP(binner==i))));
bioTot(l,i) = bioTot(l,i)+ sum(sum(sum(bioTMP(binner==i))));
mixTot(l,i) = mixTot(l,i)+ sum(sum(sum(mixTMP(binner==i))));
LayVol(l,i) = LayVol(l,i) +sum(sum(sum(volTMP(binner==i))));
end

end
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

%Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 
clear dilut surf bio tend res corr mix
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
%resTot = (resTot+Cres)*cv;
dilutTot = dilutTot*cv;
%corrTot= corrTot';
%budlats = monstart:monend;% sloppy name
%close all
%figure(61)
%plot(budlats, 100*Cres)
%caxis([-200, 200])
%hold on
%plot(budlats, tendTot);
%plot(budlats, tendTot+corrTot-dilutTot)
%plot(budlats, dilutTot, 'y')
%plot(budlats, advTot)
%plot(budlats, surfTot)
%plot(budlats, bioTot)
%plot(budlats, 100*(resTot))
%%
%legend('Storage Tend', 'Advection', 'AirSea','Bio','100*Resid')
%%
%legend('10^3*Resid', 'Tend', 'Modif Tend', 'Dilution', 'Advec', 'VolCorr', 'AirSea','Bio')

%title('Weddell DIC budget terms')%+string(lat))
%ylabel('Tg/yr')
%xlabel('5-day periods')


%saveas(gcf, 'weddelltest2.png')
%%
% plot(budlats, advTot-corrTot','b')
% plot(budlats, tendTot-dilutTot,'r')
% plot(budlats, ones(1,24)*(sum(tendTot)-sum(dilutTot))/24, '--k');
% ylim([-8e6, 4e6])
%%
%legend('10^3*Resid','Surf', 'Bio','Cor+Adv', 'Tend', 'Tend Avg');

%%



%%
save('WeddellLAYERvol'+string(Nlat)+'budgetWk'+string(monstart)+'to'+string(monend)+'Juldens.mat', 'mncurref' ,'LayVol', 'tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot');
