%%% Summing up the average quantity of carbon for a certain depth 
%%% for the Iter133 2013-2018 5day carbon data

%%% I'll average over a specified length of time, and a specified region and depth
%%% Also I'm going to plot the concentration and overlay average SSH. 
clear all
close all

load('../Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end

wkstart = 1;
wkend = 73;
z1 = 20; %190 m
z2 = 33; %750 m
zd = z2-z1+1;
budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day';

Nlat = 50;
lats = ncread('../../../data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 290;
lonE = 35;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcFull = [xc1:2160 1:xc2-1];

 % number of iterations.
 bla = size(xcFull);
 
Ctot = zeros(bla(2), yc-1); 
SSHtot = zeros(bla(2), yc-1);

for t = wkstart:wkend

    %t_array = 1;
    %t_real = t;
    %l = t- monstart+1;
    



% this will be much harder now! - each one needs three calculations!
% north boundary of the box
dic_tosum = squeeze(ncread(strcat(budgetfolder,'SnapShots_DIC.nc'), 'TRAC01', [1 1 z1 t], [Inf yc-1 zd 1]));
%dic_tosum= double(adv_tosum);
Ctot = Ctot + sum(dic_tosum(xcFull,:,:).*volume(xcFull, 1:yc-1, z1:z2),3); % because this is advection to the north!
if (t<15)
	SSH_tosum = squeeze(ncread(strcat(budgetfolder,'_SSH.nc'),'ETAN', [1 1 t], [Inf yc-1 1]));
SSHtot = SSHtot + SSH_tosum(xcFull,:,:);
end
if(t>60)

	SSH_tosum = squeeze(ncread(strcat(budgetfolder,'_SSH.nc'),'ETAN', [1 1 t] ,[Inf yc-1 1]));
SSHtot = SSHtot + SSH_tosum(xcFull,:,:);

end
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
%mixTot(l) = -sum(sum(mix_tosum(xcFull,:)));

% west boundary
%adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
%adv_tosum= double(adv_tosum);
%advTot(l) = advTot(l)+sum(sum(adv_tosum)); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
%mixTot(l) = mixTot(l)+sum(sum(mix_tosum));

% east boundary
%adv_tosum = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
%adv_tosum= double(adv_tosum);
%advTot(l) = advTot(l)-sum(sum(adv_tosum)); % because this is advection to the north!
%mix_tosum = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
%mixTot(l) = mixTot(l)-sum(sum(mix_tosum));

%%
%dilutTot(l) = sum(sum(double(dilut(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1)));
%surfTot(l) = sum(sum(double(surf(xcFull,1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1)));
%corrTot(l) = sum(sum(double(corr(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull, 1:yc-1, 1)));
%tendTot(l) = sum(sum(sum(double(tend(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
%resTot(l) = sum(sum(sum(double(res(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
%bioTot(l) = sum(sum(sum(double(bio(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

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
SSHmean = SSHtot/27;% because it's masked for sea ice(wkend-wkstart+1);
DICcont = Ctot/(wkend-wkstart +1); 
recipV = sum(volume(xcFull,1:yc-1,z1:z2),3); recipV(recipV==0) =inf; recipV = 1./recipV; 
DICconc = DICcont.*recipV; 
lines = [-1 -1.2 -1.5 -1.8  -2] ;
c = contour(lonW-360:1/6:lonE-1/6,lats(1:yc-1), SSHtot',lines);
ctr = 1;
K = size(c);


figure(137)
lev = [2.2:0.01:2.31 2.31:0.002:2.34]
contourf(gca,lonW-360:1/6:lonE-1/6, lats(1:yc-1), DICconc', 'Linecolor', 'none')

colorbar()

while ctr < K(2)
	m = c(2,ctr);
	hold on
	if m>10
		plot(gca, c(1,ctr+1:ctr+m), c(2,ctr+1:ctr+m),'b')
end
ctr = ctr+m+1;
end
%cv = 3.7843e-4; % units of Tg/yr
%Cres = Cres * cv;
%tendTot = tendTot * cv; 
%resTot = resTot *cv;
%advTot=  advTot*cv;
%surfTot = surfTot*cv;
%bioTot = bioTot*cv;
%corrTot = corrTot*cv;
%mixTot = mixTot*cv;
%resTot = (resTot+Cres)*cv;
%dilutTot = dilutTot*cv;
%corrTot= corrTot';
%budlats = monstart:monend;% sloppy name
%close(gcf)
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

saveas(gcf, 'CIWSSH1.png');
save('CIWbiggercarbonSSH'+string(wkstart)+'to'+string(wkend)+'.mat', 'DICcont' ,'DICconc', 'SSHmean')
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
