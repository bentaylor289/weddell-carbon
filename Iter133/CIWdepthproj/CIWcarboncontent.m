%%% Summing up the average quantity of carbon for a certain depth 
%%% for the Iter133 2013-2018 5day carbon data

%%% I'll average over a specified length of time, and a specified region and depth

clear all

load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end

wkstart = 322;
wkend = 347;
z1 = 20; %190 m
z2 = 33; %750 m
zd = z2-z1+1;
budgetfolder = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day';

Nlat = 50;
lats = ncread('/local/data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 290;
lonE = 35;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcFull = [xc1:2160 1:xc2-1];

 % number of iterations.
 bla = size(xcFull);
 
Ctot = zeros(bla(2), yc-1); 
for t = wkstart:wkend

    %t_array = 1;
    %t_real = t;
    %l = t- monstart+1;
    



% this will be much harder now! - each one needs three calculations!
% north boundary of the box
dic_tosum = squeeze(ncread(strcat(budgetfolder,'SnapShots_DIC.nc'), 'TRAC01', [1 1 z1 t], [Inf yc-1 zd 1]));
%dic_tosum= double(adv_tosum);
Ctot = Ctot + sum(dic_tosum(xcFull,:,:).*volume(xcFull, 1:yc-1, z1:z2),3); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
%Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 
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

DICcont = Ctot/(wkend-wkstart +1); 
recipV = sum(volume(xcFull,1:yc-1,z1:z2),3); recipV(recipV==0) =inf; recipV = 1./recipV; 
DICconc = DICcont.*recipV; 



save('CIWbiggercarbon'+string(wkstart)+'to'+string(wkend)+'.mat', 'DICcont' ,'DICconc')
%%
% plot(budlats, advTot-corrTot','b')
% plot(budlats, tendTot-dilutTot,'r')
% plot(budlats, ones(1,24)*(sum(tendTot)-sum(dilutTot))/24, '--k');
% ylim([-8e6, 4e6])
%%
%legend('10^3*Resid','Surf', 'Bio','Cor+Adv', 'Tend', 'Tend Avg');

%%



%%
%save('Weddell'+string(Nlat)+'budgetWk'+string(monstart)+'to'+string(monend)+'.mat', 'tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot','resTot');
