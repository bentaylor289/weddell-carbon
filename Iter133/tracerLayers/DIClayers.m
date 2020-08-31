%%% Computes tracers along a neutral density layer. 
%%% Goes timestep by timestep.
%% Inputs here

%latN = 61; % latitude of zonal section (S)
%latSW = 65; % latitude of stopping meridional section (once land is reached) on W border.
%latSE = 71; % " " on E border
%lonW = 302; % we'll work on the version that cross the prime meridian later.
%lonE = 25;
t1 = 1; % currently reading from monthly data but that's easy to change.
t2 = 73;

densityBounds = [23 27.7 28 28.1 28.2 28.27 28.31 28.35];
mncurref = densityBounds; 

fold5 = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
foldmon = '/local/data/bSOSE/iter133NEW/monthly/bsose_i133_2013to2018_monthly_';
lats = ncread(strcat(fold5, 'GAMMA.nc'), 'YC');
%[mini, ycN] = min(abs(lats+latN));
%[mini, ycSW] = min(abs(lats+latSW));
%[mini, ycSE] = min(abs(lats+latSE));
%xcW = 6*lonW+1; % not 100% sure if this is quite what I want ... 
%xcE = 6*lonE+1;
%if (lonW > lonE)
%	myxc = [xcW:2160 1 :xcE-1];
%end
%if (lonE > lonW)
%	myxc = xcW:xcE-1;
%end
% for zonal, I want up to but not including xcE  

load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat');
%% Get and consolidate the layer info. 

NG = length(mncurref)-1;
%DONexport = zeros(t2-t1+1,NG);


% Start with the two meridional sections!
NG = length(mncurref)-1;
maskC = hFacC; maskC(maskC>0)=1;
[NX,NY,NZ]=size(hFacS);
%maskS = hFacS; maskS(maskS>0)=1;
%creates a mask so that the averaging with the grid cell to the west
%doesn't get too messy. It's a cool tool. 
%dnmS = maskC(:,:,:) + maskC([1 1:NX-1],:,:);dnmS(dnmS==0)=inf; dnmS=1./dnmS;
lon = XC(:,1);lat = YC(1,:);r = squeeze(RC(:));
Hinz = hFacS*0;
for k = 1:NZ
  Hinz(:,:,k) = hFacW(:,:,k)*DRF(k); % using west side here as well.
end
DYinG = zeros(NX,NY,NG,'single');
for k = 1:NG
  DYinG(:,:,k) = DYG;
end

DICavg = zeros(NX,NY,NG);
Savg = zeros(NX,NY,NG);
PVavg = DICavg;
Navg = DICavg;
Havg = DICavg;
BIOavg = DICavg;
isPres = DICavg;

% to calculate f/H
ROT = 7.92e-5; % rotation rate of earth in radians

for (t = t1:t2)
 
H = zeros(NX,NY,NG);
 %NVH = VH;
 
   G = ncread(strcat(fold5, 'GAMMA.nc'), 'gamma', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'GAMMA'],iter,'rec',t)
   %V =  ncread(strcat(fold5,'Uvel.nc'), 'UVEL', [1 1 1 t] , [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   BIO = ncread(strcat(fold5,'BLGBIOC.nc'), 'BLGBIOC', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   C = ncread(strcat(fold5,'DIC.nc'), 'TRAC01', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   %MAP G,N TO W POINTS
   S = ncread(strcat(fold5,'Salt.nc'), 'SALT', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   Strat = ncread(strcat(fold5, 'Strat.nc'), 'DRHODR', [1 1 1 t], [Inf Inf Inf 1]);
   f = permute(repmat(sin(lats*pi/180), [1 2160 NZ]), [2 1 3]);
   fbyH = Strat.*f*ROT./(1000+G);
   N = ncread(strcat(foldmon,'NO3.nc'), 'TRAC04', [1 1 1 ceil(t/6)], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   G = G.*maskC;
   C = C.*maskC;
   S = S.*maskC;
   N = N.*maskC;
   fbyH = fbyH.*maskC;
   BIO = BIO.*maskC;
   %N = (N + N([1 1:NX-1],:,:)).*dnmS;
   CH = H; SH = H; fbyHH = H; NH = H; BIOH = H;

   %NH = H;
   for kcur=1:NG
     tmp1=S.*Hinz;
     tmp2=Hinz;
     tmp3=Hinz.*fbyH;
     tmp4=Hinz.*C;
     tmp5=Hinz.*N;
     tmp6=Hinz.*BIO;
     tmp1((G< mncurref(kcur)) | (G>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     tmp5(tmp1==0)=0;
     tmp6(tmp6==0)=0;
     SH(:,:,kcur) = (sum(tmp1,3)); %mass transport in layer  (vh)
     H(:,:,kcur)  = (sum(tmp2,3)); %layer thickness
     fbyHH(:,:,kcur)= (sum(tmp3,3)); %tracer transport in layer
     CH(:,:,kcur) = (sum(tmp4,3));
     NH(:,:,kcur) = (sum(tmp5,3));
     BIOH(:,:,kcur) = (sum(tmp6,3));
   end
   
   %VH = VH.*DYinG; NVH = NVH.*DYinG;
   %%
   
   recipH = H; recipH(recipH==0)=inf; recipH = 1./recipH;
   St = SH.*recipH;
   Ct = CH.*recipH; 
   fbyHt = fbyHH.*recipH;
   Nt = NH.*recipH;
   BIOt = BIOH.*recipH;
   isPrest = H; isPrest(isPrest~=0) =1;
   
Savg = Savg + St;
DICavg = DICavg + Ct;
Havg = Havg + H;
PVavg = PVavg + fbyHt;
Navg = Navg + Nt;
BIOavg = BIOavg + BIOt;
isPres = isPres + isPrest;
   end

Savg = Savg/(t2-t1+1);
DICavg = DICavg/(t2-t1+1);
PVavg = PVavg/(t2-t1+1);
Navg = Navg/(t2-t1+1);
BIOavg = BIOavg/(t2-t1+1);
Havg = Havg/(t2-t1+1);
isPres = isPres/(t2-t1+1);
Depavg = Havg; % will give depth at center of layer..
Depavg(:,:,1) = Depavg(:,:,1)/2;
for l=2:NG
	Depavg(:,:,l) = Depavg(:,:,l-1) + 0.5*Havg(:,:,l-1)+0.5*Havg(:,:,l);
end

save('BIOLayers'+string(t1)+ 'to'+string(t2) + '.mat', 'densityBounds','isPres','Navg','PVavg','BIOavg','Savg', 'DICavg', 'Havg','Depavg')


