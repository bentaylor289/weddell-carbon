%%% Weddell blocky box carbon budget, layerwise 
%%% for the Iter133 2013-2018 carbon budget 
%%% looping this one for multiple 5-day periods

clear all
load('../Iter129/grid.mat', 'hFacC', 'RAC', 'DRF','hFacW','hFacS', 'DYG','DXG');
volume = zeros(size(hFacC));
%%

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
wkstart = 1;
wkend = 438;
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
ylen = ycN-ycS;
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
cut = 2.323;

quan = ncread(strcat(budgetfolder, 'DIC.nc'),'TRAC01', [1 ycS zt monstart], [Inf ycN-ycS zl n1]);
quan = squeeze(mean(mean(quan,4),3)); % this is a very arbitrary measure... 
quan = quan(xcF,:);
quan(quan==nan)=0;

isINS = zeros(size(quan));
f = 4; % how many dots to fill each box in one dimension. I think 2 is probably fine, 4 is safe.
df = 1/(2*f);
Xh = (repmat(0.5:(xlen-0.5), [ylen 1]));
Yh = repmat(0.5:(ylen-0.5), [xlen 1])';
Xf = repmat(df:2*df:(xlen-df),[ylen*f 1])';
Yf = repmat(df:2*df:(ylen-df),[xlen*f 1]);

fillquan = interp2(Xh, Yh, quan', Xf, Yf, 'nearest');

M = contour(fillquan, [cut cut]); % produces a contour that goes thru corners.
% note that the axes for gridlines now run 0 to xlen and 0 to ylen. 
isINSsimple(quan>cut)=1;

ctr=1;
tmp = size(M);
lc = tmp(2); 
maxlen = 0;
maxctr =0;
while ctr<lc  
    m = M(2, ctr);   
    if m > maxlen
	maxlen = m;
        maxctr = ctr;
    end     
    ctr = ctr + m +1; 
end 

xx = round(2*df*M(2, maxctr+1:maxctr+maxlen));% backwards from what I'd expect
yy = round(2*df*M(1, maxctr+1:maxctr+maxlen));
xxfix = [];
yyfix = [];
Errorct = [];
for j = 2:maxlen
	dist = abs(xx(j)-xx(j-1)) + abs(yy(j)- yy(j-1));
	if (dist==1)
		xxfix = [xxfix xx(j)]; yyfix = [yyfix yy(j)];
	end
	if (dist > 1)
		Errorct = [ERRORct j];
	end
end
if (abs(xxfix(end) - xxfix(1)) + abs(yyfix(end) - yyfix(1)) ==1)
	xxfix = [xxfix xxfix(1)];
	yyfix = [yyfix yyfix(1)];
end 


isINS = inpolygon(Xh,Yh,xxfix,yyfix)';  
xx=xxfix;
yy=yyfix;
tmp = size(xx);
len = max(tmp);%length of new contour

hfC = hFacC(xcF,ycS:ycN-1,:); %???
hfS = hFacS(xcF,ycS:ycN , :);
hfW = hFacW([xcF xc2],ycS:ycN-1,:);
dy = DYG([xcF xc2], ycS:ycN-1);
dx = DXG(xcF, ycS:ycN);
%Budget initialization

ContadvTot = zeros(n1,1);
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

%
NVHtot = zeros(n1, lays);
REStot = zeros(n1, lays);
OVRtot = zeros(n1, lays);
VOLtot = zeros(n1);
VHtot = zeros(n1, lays);
Ntot = zeros(n1, lays);

% Start the computations here

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

advY = squeeze(ncread(strcat(budgetfolder,'ADVy_DIC.nc'), 'ADVyTr01', [1 ycS 1 t_real], [Inf ycN-ycS+1 Inf 1]));
advY = advY(xcF,:,:); 
advX = squeeze(ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
advX = advX([xcF xc2],:,:); 

Vvel = squeeze(ncread(strcat(budgetfolder,'Vvel.nc'), 'VVEL', [1 ycS 1 t_real], [Inf ycN-ycS+1 Inf 1]));
Vvel = Vvel(xcF,:,:); 
Uvel = squeeze(ncread(strcat(budgetfolder,'Uvel.nc'), 'UVEL', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
Uvel = Uvel([xcF xc2],:,:); 

mixY = squeeze(ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 ycS 1 t_real], [Inf ycN-ycS+1 Inf 1]));
mixY = mixY(xcF,:,:);
mixX = squeeze(ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
mixX = mixX([xcF xc2],:,:);

Gam = squeeze(ncread(strcat(budgetfolder,'GAMMA.nc'), 'gamma', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
Gam = Gam(xcF, :,:);
DIC = squeeze(ncread(strcat(budgetfolder,'DIC.nc'), 'TRAC01', [1 ycS 1 t_real], [Inf ycN-ycS Inf 1]));
DIC = DIC(xcF, :,:); 


G = zeros(len-1, 52); 
N = G; HinZ=G; V=G;
DL= zeros(len-1, 1);

advgrid= zeros(len-1,1);
for k=1:len-1
	 
    if (xx(k+1) - xx(k)==1)
	 advgrid(k) = -sum(advY(xx(k+1),yy(k)+1,:),3);
	 mixTot(l) = mixTot(l) - sum(mixY(xx(k+1), yy(k)+1,:),3);
	 Gtmp = squeeze(Gam(xx(k+1), yy(k):(yy(k)+1),:));
	 Ntmp = squeeze(DIC(xx(k+1), yy(k):(yy(k)+1),:));
	 V(k,:) = -Vvel(xx(k+1), yy(k) +1,:);
	 DL(k) = dy(xx(k+1), yy(k)+1,:);
	 HinZ(k,:) = squeeze(hfS(xx(k+1), yy(k)+1,:)).*DRF;
	 mt = squeeze(hfC(xx(k+1), yy(k):yy(k)+1,:));
	 mt2 = sum(mt,1); mt2(mt2==0)=Inf; mt2 = 1./mt2;
	 G(k,:) = sum(Gtmp.*mt,1).*mt2;
	 N(k,:) = sum(Ntmp.*mt,1).*mt2;
	 

 elseif (xx(k+1)- xx(k) ==-1)
	 advgrid(k) = sum(advY(xx(k), yy(k)+1,:),3);
	 mixTot(l) = mixTot(l) + sum(mixY(xx(k), yy(k)+1,:),3);
	
	 Gtmp = squeeze(Gam(xx(k), yy(k):(yy(k)+1),:));
	 Ntmp = squeeze(DIC(xx(k), yy(k):(yy(k)+1),:));
	 V(k,:) = Vvel(xx(k), yy(k) +1,:);
	 DL(k) = dy(xx(k), yy(k)+1,:);
	 HinZ(k,:) = squeeze(hfS(xx(k), yy(k)+1,:)).*DRF;
	 mt = squeeze(hfC(xx(k), yy(k):yy(k)+1,:));
	 mt2 = sum(mt,1); mt2(mt2==0)=Inf; mt2 = 1./mt2;
	 G(k,:) = sum(Gtmp.*mt,1).*mt2;
	 N(k,:) = sum(Ntmp.*mt,1).*mt2;


 elseif (yy(k+1) - yy(k) == 1)
	 advgrid(k) = sum(advX(xx(k)+1, yy(k+1),:),3);
	 mixTot(l) = mixTot(l) + sum(mixX(xx(k)+1, yy(k+1),:),3);

	 Gtmp = squeeze(Gam(xx(k):xx(k)+1, (yy(k)+1),:));
	 Ntmp = squeeze(DIC(xx(k):xx(k)+1, (yy(k)+1),:));
	 V(k,:) = Uvel(xx(k)+1, yy(k)+1,:);
	 DL(k) = dx(xx(k)+1, yy(k)+1,:);
	 HinZ(k,:) = squeeze(hfW(xx(k)+1, yy(k)+1,:)).*DRF;
	 mt = squeeze(hfC(xx(k):xx(k)+1, yy(k)+1,:));
	 mt2 = sum(mt,1); mt2(mt2==0)=Inf; mt2 = 1./mt2;
	 G(k,:) = sum(Gtmp.*mt,1).*mt2;
	 N(k,:) = sum(Ntmp.*mt,1).*mt2;

 elseif (yy(k+1) - yy(k) == -1)
	 advgrid(k) = -sum(advX(xx(k)+1,yy(k),:),3);
	 mixTot(l) = mixTot(l) - sum(mixX(xx(k)+1, yy(k),:),3);
	 
	 Gtmp = squeeze(Gam(xx(k):xx(k)+1, yy(k),:));
	 Ntmp = squeeze(DIC(xx(k):xx(k)+1, yy(k),:));
	 V(k,:) = -Uvel(xx(k)+1, yy(k),:);
	 DL(k) = dx(xx(k)+1, yy(k),:);
	 HinZ(k,:) = squeeze(hfW(xx(k)+1, yy(k),:)).*DRF;
	 mt = squeeze(hfC(xx(k):xx(k)+1, yy(k),:));
	 mt2 = sum(mt,1); mt2(mt2==0)=Inf; mt2 = 1./mt2;
	 G(k,:) = sum(Gtmp.*mt,1).*mt2;
	 N(k,:) = sum(Ntmp.*mt,1).*mt2;

end
 end
 ContadvTot(l) = sum(advgrid);
    
   seclen = len-1; 
   NG = lays;
   VH = zeros(seclen, NG);
   H = VH; NH= VH; NVH = VH;

   if (l==1)
	   Ncell = zeros(n1, seclen, 52);
	   Hlayer = zeros(n1, seclen, NG);
	   vcVHlayer = zeros(n1, seclen, NG);
	   VHlayer = zeros(n1, seclen, NG);
   end
   for kcur=1:NG
     tmp1=V.*HinZ;
     tmp2=HinZ;
     tmp3=V.*HinZ.*N;
     tmp4=HinZ.*N;
     tmp1( (G< mncurref(kcur)) | (G>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     
     
     VH(:,kcur) = (sum(tmp1,2)); %mass transport in layer  (vh)
     H(:,kcur)  = (sum(tmp2,2)); %layer thickness
     NVH(:,kcur)= (sum(tmp3,2)); %tracer transport in layer
     NH(:,kcur) = (sum(tmp4,2));
   end
   
   VH = VH.*DL; NVH = NVH.*DL;
   VH_l = squeeze(sum(VH,2));
   VH_i = squeeze(sum(VH,1));
   Area_i = squeeze(sum(H.*DL,1));
   
   NVH_i = squeeze(sum( NVH,1));

   recipArea = Area_i; recipArea(recipArea==0)=inf; recipArea = 1./recipArea;
   V_i = VH_i.*(recipArea); % weighted average v
   N_i = squeeze(sum(NH.*DL,1)).*recipArea;

NVHtot(l,:) = NVH_i;
VHtot(l,:) = VH_i;
Ntot(l,:) = N_i;
OVRtot(l,:) = VH_i .* N_i;
REStot(l,:) = NVHtot(l,:) - OVRtot(l,:);
VOLtot(l) = sum(VH_i)*sum(sum(NH.*DL))/sum(Area_i);
  
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
save('BlockyNVHBudget'+string(monstart)+'to'+string(monend)+'.mat', 'cut', 'isINS','NVHtot', 'VHtot', 'OVRtot', 'VOLtot', 'Ntot', 'REStot', 'tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot','resTot', 'LayVol','LayCont','tendLay', 'dilutLay', 'advLay', 'corrLay', 'surfLay', 'mixLay', 'bioLay', 'advTot', 'mncurref');
