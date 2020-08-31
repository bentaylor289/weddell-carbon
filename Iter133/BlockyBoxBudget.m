%%% Weddell blocky box carbon budget, layerwise 
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
wkend = 146;
monstart = wkstart; % lazy coding
monend = wkend;

%% first find the contour 
% for now, using high C concentration.
budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';

% making a box for the Weddell - only pull data from inside this box
Nlat = 53;
lats = ncread('../../../data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,ycN] = min(abs(lats+Nlat)); 
Slat = 75;
[min1,ycS] = min(abs(lats+Slat)); 
ycF = ycS:ycN-1;
lonW = 290;
lonE = 35;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1; % assuming our area crosses prime meridian
xcF = [xc1:2160 1:xc2-1];
xlen = (lonE+360-lonW)*6;

n1 = monend-monstart +1; % number of iterations

mncurref = [23 27.7 28 28.1 28.2 28.27 28.31 28.35 28.37 28.4 29];
lays = 10; 

%specifics for the contour
zt = 18; zl =15;
cut = 2.322;

% to salvage this approach, we'll need to interpolate the C values to the corners, then do the contour.
quan = ncread(strcat(budgetfolder, 'DIC.nc'),'TRAC01', [1 ycS zt monstart], [Inf ycN-ycS zl n1]);
quan = squeeze(mean(mean(quan,4),3)); % this is a very arbitrary measure... 
quan = quan(xcF,:);
quan(quan==nan)=0;

isINS = zeros(size(quan));
%R = (repmat(1.5:(xlen-0.5), [(ycN-ycS-1) 1]))';
%E = repmat(1.5:(ycN-ycS-0.5), [(xlen-1) 1]);
%cornerquan = interp2(quan, E,R, 'linear'); % interpolates tracer values to corners.
%M = contour(cornerquan, [cut cut]); % produces a contour that goes thru corners. 
isINS(quan>cut)=1;

%ctr=1;
%tmp = size(M);
%lc = tmp(2); 
%maxlen = 0;
%maxctr =0;
%while ctr<lc  
 %   m = M(2, ctr);   
 %   if m > maxlen
%	maxlen = m;
%        maxctr = ctr;
%    end     
%    ctr = ctr + m +1; 
%end 

%xx = M(2, maxctr+1:maxctr+maxlen);% backwards from what I'd expect
%yy = M(1, maxctr+1:maxctr+maxlen); 
%xx= [xx xx(1)];
%yy = [yy yy(1)];
%xx = [100 101 101.5 101  100];
%yy = [101 101.5 101 100.5 101];

%% Budget intro


advTot = zeros(n1,1);
mixTot = zeros(n1,1);
resTot = zeros(n1,1);
surfTot = zeros(n1,1);
bioTot = zeros(n1,1);
dilutTot = zeros(n1,1);
corrTot = zeros(n1,1);
tendTot = zeros(n1,1);
Cres = zeros(n1,1);


advLay = zeros(n1,lays);
mixLay = zeros(n1,lays);
resLay = zeros(n1,lays);
surfLay = zeros(n1,lays);
bioLay = zeros(n1,lays);
dilutLay = zeros(n1,lays);
corrLay = zeros(n1,lays);
tendLay = zeros(n1,lays);
CresLay = zeros(n1,lays);
LayVol = zeros(n1,lays);
LayCont = zeros(n1,lays);

for t = monstart:monend

    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    
load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'adv', 'mix','dilut', 'surf', 'tend', 'bio', 'res','corr');

tend(isnan(tend)) = 0;
dilut(isnan(dilut)) = 0;
surf(isnan(surf)) = 0;
bio(isnan(bio)) = 0;
res(isnan(res)) = 0;
corr(isnan(corr)) = 0;

%Starting advection and diffusion section here. 

%advY = squeeze(ncread(strcat(budgetfolder,'ADVy_DIC.nc'), 'ADVyTr01', [1 ycS 1 t_real], [Inf ycN-ycS+1 Inf 1]));
%advY = advY(xcF,:,:); 
%advX = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
%advX = advX([xcF xc2],:,:); 
%mixY = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 ycS 1 t_real], [Inf ycN-ycS+1 Inf 1]));
%mixY = mixY(xcF,:,:);
%mixX = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
%mixX = mixX([xcF xc2],:,:);

%shoelace = 0;
%x = zeros(k,1);
%adv = zeros(k,52);
%av = zeros(k,52);
%au = zeros(k,52);
%for k=1:maxlen-1
	 
 %   x(k) = floor(min(xx(k+1), xx(k)))+1;
 %   xk = floor(min(xx(k+1), xx(k)))+1;
 %   y(k) = floor(min(yy(k+1), yy(k)))+1;
 %   yk = floor(min(yy(k+1), yy(k)))+1;
    % because this grid is adjusted. 
  %  if (xx(k+1)>= xx(k)) % equals case will go to 0, so no worries.
   %     av(k,:) = advY(xk, yk,:);
%	mvk = mixY(xk, yk,:);

 %   else
  %      av(k,:) = advY(xk, yk+1,:);
   %     mvk = mixY(xk, yk+1,:);	% or tmp = min(yy(k+1), yy(k)); c= frac(c); c*v(xk, yk) + (1-c)*v(xc, yk+1)
    %end
    
%    if (yy(k+1)>= yy(k))
%        au(k,:) = advX(xk+1, yk,:);
%	muk = mixX(xk+1, yk,:);
%    else
%        au(k,:) = advX(xk, yk,:);
%	muk = mixX(xk, yk,:);
%    end
%    adv(k,:) = -(xx(k+1) - xx(k))*av(k) + (yy(k+1)-yy(k))*au(k); % these just describe the fraction of transport
%    mixk = -(xx(k+1) - xx(k))*mvk + (yy(k+1)-yy(k))*muk; % these just describe the fraction of transpotr
%    advTot(l) = advTot(l) + sum(squeeze(adv(k,:)));
%    mixTot(l) = mixTot(l) + sum(squeeze(mixk));
    
    % shoelace theorem calculates area ... just a check here
%    shoek = xx(k)*yy(k+1) - xx(k+1)*yy(k);
%    shoelace = shoelace + shoek;
    
    
%end

%%
% we should only need a few lines. just need a volume that is relevant to our region. 

if (l==1)
	Cvolume = volume(xcF, ycF, :);
	Cvolume = Cvolume.*repmat(isINS, [ 1 1 52]);
end

dilutTot(l) = sum(sum(double(dilut(xcF, ycF, 1, t_array)).*Cvolume(:,:,1)));
surfTot(l) = sum(sum(double(surf(xcF,ycF, 1, t_array)).*Cvolume(:,:,1)));
corrTot(l) = sum(sum(double(corr(xcF, ycF, 1, t_array)).*Cvolume(:, :, 1)));
tendTot(l) = sum(sum(sum(double(tend(xcF,ycF, :, t_array)).*Cvolume)));
resTot(l) = sum(sum(sum(double(res(xcF,ycF, :, t_array)).*Cvolume)));
bioTot(l) = sum(sum(sum(double(bio(xcF,ycF, :, t_array)).*Cvolume)));
advTot(l) = sum(sum(sum(double(adv(xcF,ycF, :, t_array)).*Cvolume)));
mixTot(l) = sum(sum(sum(double(mix(xcF,ycF, :, t_array)).*Cvolume)));

G = squeeze(ncread(strcat(budgetfolder, 'GAMMA.nc'), 'gamma', [1 1 1 t_real], [Inf ycN-1 Inf 1]));
G1 = G(xcF,ycF,:);
binner = discretize(G1, mncurref);
C = squeeze(ncread(strcat(budgetfolder, 'DIC.nc'), 'TRAC01', [1 1 1 t_real], [Inf ycN-1 Inf 1]));

dilutTMP = double(dilut(xcF, ycF, 1, t_array)).*Cvolume(:,:,1);
surfTMP = double(surf(xcF,ycF, 1, t_array)).*Cvolume(:,:,1);
corrTMP = double(corr(xcF, ycF, 1, t_array)).*Cvolume(:, :, 1);
tendTMP = double(tend(xcF,ycF, :, t_array)).*Cvolume;
mixTMP = double(mix(xcF,ycF, :, t_array)).*Cvolume;
bioTMP = double(bio(xcF,ycF, :, t_array)).*Cvolume;
advTMP = double(adv(xcF,ycF, :, t_array)).*Cvolume;
contTMP = C(xcF, ycF, :, t_array).*Cvolume;


for i = 1:lays
dilutLay(l,i) = sum(sum(dilutTMP(binner(:,:,1)==i)));
surfLay(l,i) = sum(sum(surfTMP(binner(:,:,1)==i)));
corrLay(l,i) = sum(sum(corrTMP(binner(:,:,1)==i)));
tendLay(l,i) = sum(sum(sum(tendTMP(binner==i))));
bioLay(l,i) = sum(sum(sum(bioTMP(binner==i))));
mixLay(l,i) = sum(sum(sum(mixTMP(binner==i))));
advLay(l,i) = sum(sum(sum(advTMP(binner==i))));
LayVol(l,i) = sum(sum(sum(Cvolume(binner==i))));
LayCont(l,i) = sum(sum(sum(contTMP(binner==i))));
end

%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 
clear adv mix dilut surf bio tend res corr
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


saveas(gcf, 'weddelltestcont.png')
%%
% plot(budlats, advTot-corrTot','b')
% plot(budlats, tendTot-dilutTot,'r')
% plot(budlats, ones(1,24)*(sum(tendTot)-sum(dilutTot))/24, '--k');
% ylim([-8e6, 4e6])
%%
%legend('10^3*Resid','Surf', 'Bio','Cor+Adv', 'Tend', 'Tend Avg');

%%



%%
save('BlockyLayerBudget'+string(monstart)+'to'+string(monend)+'.mat', 'isINS','tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot','resTot', 'LayVol','LayCont','tendLay', 'dilutLay', 'advLay', 'corrLay', 'surfLay', 'mixLay', 'bioLay', 'advTot', 'mncurref');
