%%% General carbon content budget;  
%%% Ben Taylor 09/02/2020
%%% for the B-SOSE Iteration 133 2013-2018 carbon budget 
%%% looping this one for multiple 5-day periods

%%% 

clear all
%% inputs to change

wkstart = 1; % best to start with week (really, 5day step) #1 
wkend = 1; % 438 is the final week (6 years, 73 5-day periods in a year)

% folders of bSOSE data
fold5day = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
fold1dy = '../../../data/bSOSE/iter133NEW/1day/bsose_i133_2013to2018_1dy_';
fold5daySnaps = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5daySnaps_';
fold5daySnapShots = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5daySnapShots_';
% Other dependencies: cellwise budget folder (iter133/cellwise5day/), grid file. 

Nlat = 60; % S is assumed, so enter positive latitudes
lonW = 160; 
lonE = 230; % using 0 thru 360, not -180 to 180
% if crossing the prime meridian, fine to leave lonW > lonE (this case is covered).

ycSW = 75; % to avoid reading data from areas of land, you can set the point where meridional sections hit the coast here
ycSE = 75; % not necessary - can set 78 (or higher value) as a default and the code should run fine.  

% Using the sea ice budget adds significant error to the budget
% Not using the budget assumes that sea ice production/melt within the region is exactly balanced by export and change in total sea ice mass; that is, data for sea ice production/melt is combined with residual (i.e. not explicitlyused).  
UseSeaIceBudget = true; 

% code for calculating horizontal and overturning is commented out, but is easily included. 
%densityBounds = [23 27.7 28 28.1 28.2 28.27 28.31 28.35 28.37 28.4 29]; % uncomment this for density bounds

% code for "Weddellfix" not included here; individual budget codes contain it, so code can be copy-pasted from there if needed
% "Weddellfix" is used to accomodate regions where Antarctic coastline crosses a line of longitude multiple times
% for example, 59 W crosses the Antartic peninsula but also cuts off some of the Weddell coastal shelf west of 59 W% "Weddellfix" was implemented to inlucde that marginal region of coastal shelf in the Weddell region.

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


volume = zeros(size(hFacC)); % volume of ocean in each grid cell
for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
sur = volume(:,:,1)/DRF(1); % surface area of ocean in each grid cell

% Tracer box budget total terms
advTot = zeros(n1,1);
mixTot = zeros(n1,1); % parametrized mixing (small away from coast and mixed layer, integrates to ~ 0)
resCell = zeros(n1,1);%integrated residual from the cellwise budget, which should be tiny!
surfTot = zeros(n1,1); % airsea fluxes
bioTot = zeros(n1,1);
dilutTot = zeros(n1,1);
corrTot = zeros(n1,1); % correction for advection
tendTot = zeros(n1,1); % change in concentration*volume, integrated thru the box
resBox = zeros(n1,1); %residual of the integrated volume*concentration budget for the whole box

% Volume Budget total terms
advVol = zeros(n1,1); % advection of volume into the region
SSHVol = zeros(n1,1); % change in total volume of region from beginning snapshot to end snapshot
FWVol = zeros(n1,1); % total FW input into that region
resVol = zeros(n1,1); % residual of volume budget

% Sea Ice Budget total terms
expSI = zeros(n1,1); % total sea ice export
meltSI = zeros(n1,1); % total sea ice production
heffSI = zeros(n1,1); % change in sea ice effective height from beginning to end
resSI  = zeros(n1,1); % residual of the sea ice budget 

phi0 = zeros(n1,1); % average carbon concentration on boundary (full-depth)
resSub5day = zeros(n1,1); % difference between online and offline advection. If it's large, there's a problem. 

for t = monstart:monend

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
adv_tosum = squeeze(ncread(strcat(fold5day,'ADVy_DIC.nc'), 'ADVyTr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = -sum(sum(adv_tosum(xcFull,:))); % negative because this is advection to the north!
mix_tosum = squeeze(ncread(strcat(fold5day,'DFyE_DIC.nc'), 'DFyETr01', [1 yc 1 t_real], [Inf 1 Inf 1]));
mixTot(l) = -sum(sum(mix_tosum(xcFull,:)));

% west boundary
adv_tosum = squeeze(ncread(strcat(fold5day,'ADVx_DIC.nc'), 'ADVxTr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = advTot(l)+sum(sum(adv_tosum));
mix_tosum = squeeze(ncread(strcat(fold5day,'DFxE_DIC.nc'), 'DFxETr01', [xc1 1 1 t_real], [1 yc-1 Inf 1]));
mixTot(l) = mixTot(l)+sum(sum(mix_tosum));

% east boundary
adv_tosum = squeeze(ncread(strcat(fold5day,'ADVx_DIC.nc'), 'ADVxTr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
advTot(l) = advTot(l)-sum(sum(adv_tosum)); 
mix_tosum = squeeze(ncread(strcat(fold5day,'DFxE_DIC.nc'), 'DFxETr01', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
mixTot(l) = mixTot(l)-sum(sum(mix_tosum));

%%
% sum concentration change * cell volume over full region. 
dilutTot(l) = sum(sum(double(dilut(xcFull, 1:yc-1, 1, 1)).*volume(xcFull,1:yc-1,1)));
surfTot(l) = sum(sum(double(surf(xcFull,1:yc-1, 1, 1)).*volume(xcFull,1:yc-1,1)));
corrTot(l) = sum(sum(double(corr(xcFull, 1:yc-1, 1, 1)).*volume(xcFull, 1:yc-1, 1)));
tendTot(l) = sum(sum(sum(double(tend(xcFull,1:yc-1, :, 1)).*volume(xcFull,1:yc-1,:))));
resCell(l) = sum(sum(sum(double(res(xcFull,1:yc-1, :, 1)).*volume(xcFull,1:yc-1,:))));
bioTot(l) = sum(sum(sum(double(bio(xcFull,1:yc-1, :, 1)).*volume(xcFull,1:yc-1,:))));

resBox(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l) - mixTot(l)+ corrTot(l); 
clear dilut surf bio tend res corr
    
%% volume budget begins here. 

adv_tosum = squeeze(ncread(strcat(fold5day,'Vvel.nc'), 'VVEL', [1 yc 1 t_real], [Inf 1 Inf 1]));
adv_tosum= double(adv_tosum);
hFac = squeeze(hFacS(xcFull, yc,:)); 
advVol(l) = -sum(sum(adv_tosum(xcFull,:).*hFac.*DXG(xcFull,yc).*repmat(DRF,[1 FullLen])')); % negative because model output is advection to the north, and we count into box as positive

% west boundary
adv_tosum = squeeze(ncread(strcat(fold5day,'Uvel.nc'), 'UVEL', [xc1 ycSW 1 t_real], [1 yc-ycSW Inf 1]));
adv_tosum= double(adv_tosum);
hFac = squeeze(hFacW(xc1, ycSW:yc-1,:));
advVol(l) = advVol(l)+sum(sum(adv_tosum.*hFac.*DYG(xc1, ycSW:yc-1)'.*repmat(DRF,[1 yc-ycSW])')); 

% east boundary
adv_tosum = squeeze(ncread(strcat(fold5day,'Uvel.nc'), 'UVEL', [xc2 1 1 t_real], [1 yc-1 Inf 1]));
adv_tosum= double(adv_tosum);
hFac = squeeze(hFacW(xc2, 1:yc-1,:));
advVol(l) = advVol(l)-sum(sum(adv_tosum.*hFac.*DYG(xc2, 1:yc-1)'.*repmat(DRF, [1 yc-1])')); 

SSHin = ncread(strcat(fold5daySnaps, 'SSH.nc'), 'ETAN', [1 1 t_real], [Inf Inf 2] ); % Eta snapshots at beginning and end
SSH = diff(SSHin,1,3)/(3600*120); % difference between adjacent snapshots

FWin = ncread(strcat(fold1day, 'oceFWflx.nc'), 'oceFWflx', [1 1 5*t_real-4], [Inf Inf 5]); % total freshwater flux into the region, daily output for the 5 days of the period
FW = mean(FWin, 3)/1035; % averaging the daily values across the 5-day period

SSHVol(l) = sum(sum(double(SSH(xcFull, 1:yc-1)).*sur(xcFull,1:yc-1,1))); % multiply by cellwise sea area
FWVol(l) = sum(sum(double(FW(xcFull,1:yc-1)).*sur(xcFull,1:yc-1,1)));

resVol(l) = SSHVol(l) - FWVol(l) - advVol(l); 

%% Sea Ice Budget begins here
% makes use of FW input file loaded in volume budget. 
    
adv_tosum = squeeze(ncread(strcat(fold1day,'SIvheff.nc'), 'SIvheff', [1 yc 5*t_real-4], [Inf 1 5]));
adv_tosum= squeeze(mean(adv_tosum,2));
expSI(l) = +sum(adv_tosum(xcFull).*DXG(xcFull,yc)); % export is northward (signs opposite from vol budget)

% west boundary
adv_tosum = squeeze(ncread(strcat(fold1day,'SIuheff.nc'), 'SIuheff', [xc1 ycSW  5*t_real-4], [1 yc-ycSW 5]));
adv_tosum= squeeze(mean(adv_tosum,2));
expSI(l) = expSI(l)-sum(adv_tosum.*DYG(xc1, ycSW:yc-1)');

% east boundary
adv_tosum = squeeze(ncread(strcat(fold1day,'SIuheff.nc'), 'SIuheff', [xc2 1 5*t_real-4], [1 yc-1 5]));
adv_tosum= squeeze(mean(adv_tosum,2));
expSI(l) = expSI(l)+sum(adv_tosum.*DYG(xc2, 1:yc-1)');


hSIin = ncread(strcat(fold5daySnaps, 'SeaIceHeff.nc'), 'SIheff', [1 1 t_real], [Inf yc 2] ); % snapshots of sea ice thickness
SIh = squeeze(diff(hSIin,1,3))/(3600*120); % divide by number of seconds; 3 in diff refers to third dimension.

FW2in = ncread(strcat(fold1day, 'SIatmFW.nc'), 'SIatmFW', [1 1 5*t_real-4], [Inf Inf 5]); % freshwater input without sea ice. This is not clearly communicated in the ncinfo, and led to some confusion early on. Feel free to convince yourself that this is what this data is. 
FW2 = mean(FW2in, 3)/1035; 

heffSI(l) = sum(sum(double(SIh(xcFull, 1:yc-1)).*sur(xcFull,1:yc-1)));
FWTot2(l) = sum(sum(double(FW2(xcFull,1:yc-1)).*sur(xcFull,1:yc-1))); % freshwater fluxes, not including sea ice.

meltSI(l) = FMTot(l) - FWTot2(l); % sea ice melt is difference of FW input with and without sea ice. 
resSI(l) = heffSI(l) + meltSI(l) + expSI(l); 

%% calculate phi0, as well as horizontal and overturning as desired

% This is somewhat inefficient, as you could just read these small datasets for all values of t (and avoid repeated ncreads; earlier budget sections do better with one ncread per 5day period). 

% West section
   G1 = squeeze(ncread(strcat(folder, 'GAMMA.nc'), 'gamma', [xcW-1 ycSW 1 t_real] , [2 ycN-ycSW+1 Inf 1]));% two values in order to get the average
   V1 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'UVEL' , [xcW ycSW 1 t_real], [1 ycN-ycSW+1  Inf 1])); 
   N1 = squeeze(ncread(strcat(folder, 'DIC.nc'), 'TRAC01', [xcW-1 ycSW 1 t_real], [2 ycN-ycSW+1 Inf 1])); 
   mask = squeeze(ncread(strcat(folder, 'DIC.nc'), 'hFacC', [xcW-1 ycSW 1], [2 ycN-ycSW+1 Inf ]));
   Hinz1 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'hFacW', [xcW ycSW 1], [1 ycN-ycSW+1 Inf]));
   DYinG1 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'dyG', [xcW ycSW], [1 ycN-ycSW+1]));
   Hinz1 = Hinz1.*DRF';
   % Averaging to get tracer values along cell boundary as opposed to the center of the grid cell
   % mask1 is used to deal with the case where one or both cells isn't fully ocean (ie hFacC~= 1). 
   mask1 = sum(mask,1); mask1(mask1==0) = Inf; mask1 = 1./mask1;
   G1 = squeeze(sum(mask.*G1,1).*mask1);
   N1 = squeeze(sum(mask.*N1,1).*mask1);
   
   % East section
   G3 = squeeze(ncread(strcat(folder, 'GAMMA.nc'), 'gamma', [xcE-1 ycSE 1 t_real], [2 ycN-ycSE+1 Inf 1]));
   V3 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'UVEL' , [xcE ycSE 1 t_real], [1 ycN-ycSE+1 Inf 1])); 
   N3 = squeeze(ncread(strcat(folder, 'DIC.nc'), 'TRAC01', [xcE-1 ycSE 1 t_real], [2 ycN-ycSE+1 Inf 1])); 
   mask = squeeze(ncread(strcat(folder, 'DIC.nc'), 'hFacC', [xcE-1 ycSE 1], [2 ycN-ycSE+1 Inf ])); 
   Hinz3 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'hFacW', [xcE ycSE 1], [1 ycN-ycSE+1 Inf]));
   DYinG3 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'dyG', [xcE ycSE], [1 ycN-ycSE+1]));
   Hinz3 = Hinz3.*DRF';
   mask3 = sum(mask,1); mask3(mask3==0) = Inf; mask3 = 1./mask3;
   G3 = squeeze(sum(mask.*G3,1).*mask3);
   N3 = squeeze(sum(mask.*N3,1).*mask3);
   
   % North section begins here
   G = squeeze(ncread(strcat(folder, 'GAMMA.nc'),'gamma', [1 ycN 1 t_real], [Inf 2 Inf 1])); % grabbing full circumpolar data to deal with the case where the section crosses the prime meridian.    
   V = squeeze(ncread(strcat(folder, 'Vvel.nc'), 'VVEL' , [1 ycN+1 1 t_real], [Inf 1 Inf 1])); 
   N = squeeze(ncread(strcat(folder, 'DIC.nc'), 'TRAC01', [1 ycN 1 t_real], [Inf 2 Inf 1])); 
   mask = squeeze(ncread(strcat(folder, 'DIC.nc'), 'hFacC', [1 ycN 1 ], [Inf 2 Inf])); 
   Hinz2 = squeeze(ncread(strcat(folder, 'Vvel.nc'), 'hFacS', [1 ycN+1 1], [Inf 1 Inf]));
   DXinG2 = squeeze(ncread(strcat(folder, 'Vvel.nc'), 'dxG', [1 ycN+1], [Inf 1]));       
   Hinz2 = Hinz2.*DRF';

   mask2 = sum(mask,2); mask2(mask2==0) = Inf; mask2 = 1./mask2; 
   G2 = squeeze(sum(mask.*G,2).*mask2);
   N2 = squeeze(sum(mask.*N,2).*mask2);
   norbdx = xcFull; % taking only the part within our region 
   V2 = V(norbdx,:);
   N2 = N2(norbdx,:);
   G2 = G2(norbdx,:);
   Hinz2 = Hinz2(norbdx,:);
   DXinG2 = DXinG2(norbdx);

% combining the sections into a single array   
   G = [G1' G2' fliplr(G3)']';
   N = [N1' N2' fliplr(N3)']';
   V = [V1' -V2' -fliplr(V3)']';

   Hinz = [Hinz1' Hinz2' Hinz3']';
   DL = [DYinG1 DXinG2' DYinG3]';

   seclen = size(G); 
   seclen = seclen(1);
   VH = zeros(seclen, NG);
   H = VH; NH= VH; NVH = VH;

   NVHTot(l) = sum(sum(V.*HinZ.*N,2).*DL); % N*V*Area for whole boundary
   resSub5day(l) = advTot(l) - NVHTot(l); 
   phi0(l) = sum(sum(HinZ.*N,2).*DL)/sum(sum(HinZ,2).*DL); % Area-weighted average tracer concentration
   
   %% Overturning calculation 
   % mncurref = densityBounds

   %for kcur=1:NG
   %  tmp1=V.*Hinz;
   %  tmp2=Hinz;
   %  tmp3=V.*Hinz.*N;
   %  tmp4=Hinz.*N;
   %  tmp1( (G< mncurref(kcur)) | (G>=mncurref(kcur+1)) )=0;
   %  tmp2(tmp1==0)=0;
   %  tmp3(tmp1==0)=0;
   %  tmp4(tmp1==0)=0;
     
     
   %  VH(:,kcur) = (sum(tmp1,2)); %mass transport in layer  (vh)
   %  H(:,kcur)  = (sum(tmp2,2)); %layer thickness
   %  NVH(:,kcur)= (sum(tmp3,2)); %tracer transport in layer
   %  NH(:,kcur) = (sum(tmp4,2));
   %end
   
   %VH = VH.*DL; NVH = NVH.*DL;
   %VH_l = squeeze(sum(VH,2));
   %VH_i = squeeze(sum(VH,1));
   %Area_i = squeeze(sum(H.*DL,1));
   
   %NVH_i = squeeze(sum( NVH,1));
   %recipArea = Area_i; recipArea(recipArea==0)=inf; recipArea = 1./recipArea;
   %V_i = VH_i.*(recipArea); % weighted average v
   %N_i = squeeze(sum(NH.*DL,1)).*recipArea;

   %NVHtot(t,:) = NVH_i;
   %VHtot(t,:) = VH_i;
   %Ntot(t,:) = N_i;
   %OVRtot(t,:) = VH_i .* N_i;

%% directly evaluate change in surface content  

etaSn = squeeze(ncread(strcat(fold5daySnaps, 'SSH.nc'), 'ETAN', [1 1 t_real], [Inf yc-1 2])); 
DICsur = squeeze(ncread(strcat(fold5day, 'DIC.nc'), 'TRAC01', [1 1 1 t_real],[Inf yc-1 1 1]));

partetaC(l) = 1/(120*3600)* sum(sum(squeeze(diff(etaSn(xcFull, :,:),1,3)).*sur(xcFull, 1:yc-1, 1).*DICsur(xcFull,:,1))); % d(eta)/dt * area * average DIC concentration gives change in surface carbon content, exactly.

end 

%% Now we combine these to form our content budget. 
% I'll refer to the Carbon Inventory Budget Latex here to connect code to equations there. 
% for now the units are mol/s

eps_damp = resVol; % equation 24 (volume budget). 
eps_vol = partetaC - corrTot + dilutTot - resVol.*phi0; % equation 26
eps_SI = resSI; % equation 28
v_cons_adv = advTot - phi0.*advVol; % volume conserved advection, which appears in equation 29 as the difference between total advection and phi0*v0*A_c. 
% we could also evaluate this using NVHtot - VOLtot, or by creating a volume-conserved velocity field. 
% (See GeneralBoxHorizontal.m code) These are equal up to the sub5day residual

storage = tendTot + partetaC - SSHVol.*phi0 - heffSI.*phi0; % equation 29, LHS (not including residuals)
% SSHVol.*phi0 is the d(eta)/dt * phi_0

if (UseSeaIceBudget) % proceed as in the LaTeX, without hiding the residual
PER = FWVol - meltSI; % equivalent to total FW2. Appears in equation 29, multipled by phi0 as freshwater outflow term

resExtra = storage - eps_vol - eps_SI.*phi0 - surfTot - bioTot -(v_cons_adv) - expSI.*phi0 + PER.*phi0 - mixTot- resBox;
% from equation 29, the extra residual not attributable to the carbon concentration, volume, or sea ice budgets   
resTotal = eps_vol + eps_SI.*phi0 + resBox + resExtra;
end 

if (~UseSeaIceBudget) % let's assume that sea ice export and height change are well captured, and the residual has to do with precipitation (i.e. let's include the eps_SI in meltSI).
meltSI = meltSI - resSI; % include the residual in meltSI, so that now meltSI + expSI + heffSI = 0	
PER = FWVol - meltSI; % no longer equivalent to FW2. 

% I suspect there's a sign backwards here
resExtra = storage - eps_vol - surfTot - bioTot -(v_cons_adv) - expSI.*phi0 + PER.*phi0 - mixTot- resBox;
% from equation 29, the extra residual not attributable to the carbon concentration, volume, or sea ice budgets   
resTotal = eps_vol  + resBox + resExtra;
end 

% At the end, convert everything into more familiar units! 
cv = 3.7843e-4; % units of Tg/yr, instead of mol/s . 
storage = storage * cv; 
resBox = resBox *cv;
resTotal = resTotal*cv;
resExtra = resExtra*cv;
surfTot = surfTot*cv;
bioTot = bioTot*cv;
v_cons_adv = v_cons_adv*cv;
mixTot = mixTot*cv;
eps_vol = eps_vol*cv;
% to compare eps_SI with the other residuals, multiply by phi0 (and by cv as needed)

FW_driven_export = phi0.*PER*cv; 
SI_driven_import = phi0.*expSI*cv;


save('GContentBudgetWk'+string(monstart)+'to'+string(monend)+'.mat', 'Nlat', 'lonE','lonW','resTotal','resExtra','FW_driven_export','SI_driven_import','storage','v_cons_adv', 'tendTot', 'dilutTot', 'advTot', 'corrTot', 'surfTot','bioTot', 'mixTot','resBox');
