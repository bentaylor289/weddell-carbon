%%% Linear decomposition of advection
%%% modified May 7 2020 to work for iter133 bSOSE DON output
%%% The idea here is to replicate the calculations that Graeme made, so
%%% it's easy to direct compare the results. 
%%% It works for Weddell Gyre right now, I hope!

%%% This version combines the sections into one, which will help!
%%% I also aspire to have two ways to calculate each relevant quantity. 


 latN = 61; % latitude of zonal section (S)
 latSW = 63.56; % latitude of stopping meridional section (once land is reached) on W border.
 latSE = 71; % " " on E border
 lonW = 302; % we'll work on the version that cross the prime meridian later.
 lonE = 25; %if either is 0 it's gets messy
%t = 1; % currently reading from monthly data but that's easy to change.

densityBounds = [23 27.5 27.7 27.9 28 28.1 28.2 28.25 28.3 28.4 28.5 29];
mncurref = densityBounds; 
folder = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';

 lats = ncread(strcat(folder, 'GAMMA.nc'), 'YC');
 [mini, ycN] = min(abs(lats+latN));
 [mini, ycSW] = min(abs(lats+latSW));
 [mini, ycSE] = min(abs(lats+latSE));
 xcW = 6*lonW+1; % not 100% sure if this is quite what I want ... 
 xcE = 6*lonE+1;
 if (lonW > lonE)
	 myxc = [xcW:2160 1:xcE-1];
 end
 if (lonE > lonW)
	 myxc = [xcW:xcE-1];
end
	 %xcW = 1813;
%xcE = 180; 

load('../Iter129/grid.mat');

%% Get and consolidate the layer info. 

% Start with the two meridional sections!
NG = length(mncurref)-1;

tf = 438;
NVHtot = zeros(tf,NG);
OVRtot = zeros(tf,NG);
REStot = zeros(tf,NG);
VOLtot = zeros(tf,1);
VHtot = zeros(tf,NG);
Ntot = zeros(tf,NG);
% volume conserved stuff
vcNVH = zeros(tf,NG);
vcOVR = zeros(tf,NG);
vcRES = zeros(tf,NG);

for t = 1:tf
% West section
   G1 = squeeze(ncread(strcat(folder, 'GAMMA.nc'), 'gamma', [xcW-1 ycSW 1 t] , [2 ycN-ycSW+1 Inf 1]));
   V1 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'UVEL' , [xcW ycSW 1 t], [1 ycN-ycSW+1  Inf 1])); %dmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N1 = squeeze(ncread(strcat(folder, 'DON.nc'), 'TRAC07', [xcW-1 ycSW 1 t], [2 ycN-ycSW+1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   mask = squeeze(ncread(strcat(folder, 'DIC.nc'), 'hFacC', [xcW-1 ycSW 1], [2 ycN-ycSW+1 Inf ]));
   Hinz1 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'hFacW', [xcW ycSW 1], [1 ycN-ycSW+1 Inf]));
   DYinG1 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'dyG', [xcW ycSW], [1 ycN-ycSW+1]));
   Hinz1 = Hinz1.*DRF';
   % Use the mask to get the tracer values on the grid cell line
   mask1 = sum(mask,1); mask1(mask1==0) = Inf; mask1 = 1./mask1;
   G1 = squeeze(sum(mask.*G1,1).*mask1);
   N1 = squeeze(sum(mask.*N1,1).*mask1);
   
   % East section
   G3 = squeeze(ncread(strcat(folder, 'GAMMA.nc'), 'gamma', [xcE-1 ycSE 1 t], [2 ycN-ycSE+1 Inf 1]));
   V3 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'UVEL' , [xcE ycSE 1 t], [1 ycN-ycSE+1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N3 = squeeze(ncread(strcat(folder, 'DON.nc'), 'TRAC07', [xcE-1 ycSE 1 t], [2 ycN-ycSE+1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   mask = squeeze(ncread(strcat(folder, 'DIC.nc'), 'hFacC', [xcE-1 ycSE 1], [2 ycN-ycSE+1 Inf ])); 
   Hinz3 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'hFacW', [xcE ycSE 1], [1 ycN-ycSE+1 Inf]));
   DYinG3 = squeeze(ncread(strcat(folder, 'Uvel.nc'), 'dyG', [xcE ycSE], [1 ycN-ycSE+1]));
   Hinz3 = Hinz3.*DRF';
   mask3 = sum(mask,1); mask3(mask3==0) = Inf; mask3 = 1./mask3;
   G3 = squeeze(sum(mask.*G3,1).*mask3);
   N3 = squeeze(sum(mask.*N3,1).*mask3);
   
   % North section begins here
   G = squeeze(ncread(strcat(folder, 'GAMMA.nc'),'gamma', [1 ycN 1 t], [Inf 2 Inf 1]));  % this is the one we like. 
   V = squeeze(ncread(strcat(folder, 'Vvel.nc'), 'VVEL' , [1 ycN+1 1 t], [Inf 1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N = squeeze(ncread(strcat(folder, 'DON.nc'), 'TRAC07', [1 ycN 1 t], [Inf 2 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4
   mask = squeeze(ncread(strcat(folder, 'DIC.nc'), 'hFacC', [1 ycN 1 ], [Inf 2 Inf])); 
   Hinz2 = squeeze(ncread(strcat(folder, 'Vvel.nc'), 'hFacS', [1 ycN+1 1], [Inf 1 Inf]));
   DXinG2 = squeeze(ncread(strcat(folder, 'Vvel.nc'), 'dxG', [1 ycN+1], [Inf 1]));       
   Hinz2 = Hinz2.*DRF';

   mask2 = sum(mask,2); mask2(mask2==0) = Inf; mask2 = 1./mask2;
   G2 = squeeze(sum(mask.*G,2).*mask2);
   N2 = squeeze(sum(mask.*N,2).*mask2);
   norbdx = myxc; 
   V2 = V(norbdx,:);
   N2 = N2(norbdx,:);
   G2 = G2(norbdx,:);
   Hinz2 = Hinz2(norbdx,:);
   DXinG2 = DXinG2(norbdx);

% combining the sections into a single array   
   G = [G1' G2' G3']';
   N = [N1' N2' N3']';
   V = [V1' -V2' -V3']';

   Hinz = [Hinz1' Hinz2' Hinz3']';
   DL = [DYinG1 DXinG2' DYinG3]';
   
   seclen = size(G); 
   seclen = seclen(1);
   VH = zeros(seclen, NG);
   H = VH; NH= VH; NVH = VH;
   for kcur=1:NG
     tmp1=V.*Hinz;
     tmp2=Hinz;
     tmp3=V.*Hinz.*N;
     tmp4=Hinz.*N;
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

NVHtot(t,:) = NVH_i;
VHtot(t,:) = VH_i;
Ntot(t,:) = N_i;
OVRtot(t,:) = VH_i .* N_i;
REStot(t,:) = NVHtot(t,:) - OVRtot(t,:);
VOLtot(t) = sum(VH_i)*sum(sum(NH.*DL))/sum(Area_i);
  
   % dealing with barotropic sloshing
   Vbar = sum(VH_i)/sum(Area_i); % barotropic v
   V1 = V1 - Vbar; 
   V2 = -V2 - Vbar;
   V3 = -V3 - Vbar;
   V = [V1' V2' V3']'; % new volume-conserved v field

   % for now we'll just repeat all the same steps with a new v
   VH = zeros(seclen, NG);
   H = VH; NH= VH; NVH = VH;
   for kcur=1:NG
     tmp1=V.*Hinz;
     tmp2=Hinz;
     tmp3=V.*Hinz.*N;
     tmp4=Hinz.*N;
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

vcNVH(t,:) = NVH_i;
vcOVR(t,:) = N_i.*VH_i;
vcRES(t,:) = NVH_i - N_i.*VH_i;

% directly calculate the v'C' term with volume-conserved V
RecipHDL = H.*DL; RecipHDL(RecipHDL==0) = inf; RecipHDL = 1./RecipHDL;
VprHL = VH - V_i.*H.*DL; 
NprHL = NH.*DL - N_i.*H.*DL;
NprVprH = VprHL.*NprHL.*RecipHDL;

end


%save('LayerwiseTest'+string(tf)+'NEW.mat', 'VH1_y','VH2_x','VH3_y', 'N1_y', 'N2_x', 'N3_y','Area1','Area2','Area3','NVH1','NVH2','NVH3', 'cont1' , 'cont2' , 'cont3','NVHtot', 'OVRtot', 'REStot', 'VOLtot')

save('LayerCon'+string(tf)+'Wed61DON.mat', 'NVHtot', 'OVRtot','REStot', 'VOLtot','VHtot','Ntot','vcNVH', 'vcOVR', 'vcRES', 'NprVprH')

