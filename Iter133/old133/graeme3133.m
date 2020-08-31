%%% Graeme's figures plot
%%% modified March 19 2020 to work for iter133b bSOSE output
%%% The idea here is to replicate the calculations that Graeme made, so
%%% it's easy to direct compare the results. 
%%% It works for Weddell Gyre right now, I hope!

%%% To make my life easier, at this point it runs only at a snapshot. 

%%% we'll want to make it be able to run multiple times though. 
%%% If I want to see what happens over time, I can run it multiples times lol.
%%% that'll have to change lol! 

%%% Next steps would probably be to make this capable of going across the prime meridian
%%% Also adding in 3day capabilities. Also a bunch of sophisticated stuff 
%%% added by budgeting would be cool
%% Inputs here
% latN = ; % latitude of zonal section (S)
% latSW = 70; % latitude of stopping meridional section (once land is reached) on W border.
% latSE = 77; % " " on E border
% lonW = 165; % we'll work on the version that cross the prime meridian later.
% lonE = 210;
%t = 1; % currently reading from monthly data but that's easy to change.

densityBounds = [23 27.7 27.9 28 28.1 28.13 28.15 28.2 28.23 28.25 29];
mncurref = densityBounds; 

% lats = ncread('../../../data/bSOSE/iter105/monthly/monDIC.nc', 'YC');
% [mini, ycN] = min(abs(lats+latN));
% [mini, ycSW] = min(abs(lats+latSW));
% [mini, ycSE] = min(abs(lats+latSE));
% xcW = 3*lonW+1; % not 100% sure if this is quite what I want ... 
% xcE = 3*lonE+1;

ycN = 334; % not used any more
xcW = 1813;
xcE = 180; 
folder = '../../../data/bSOSE/iter133OLD/sections/';

load grid.mat;

%% Get and consolidate the layer info. 

% Start with the two meridional sections!
NG = length(mncurref)-1;
% maskC = hFacC; maskC(maskC>0)=1;
% [NX,NY,NZ]=size(hFacS);
% maskS = hFacS; maskS(maskS>0)=1;
% %creates a mask so that the averaging with the grid cell to the west
% %doesn't get too messy. It's a cool tool. 
% dnmS = maskC(:,:,:) + maskC([1 1:NX-1],:,:);dnmS(dnmS==0)=inf; dnmS=1./dnmS;
% lon = XC(:,1);lat = YC(1,:);r = squeeze(RC(:));
% Hinz = hFacS*0;
% for k = 1:NZ
%   Hinz(:,:,k) = hFacW(:,:,k)*DRF(k); % using west side here as well.
% end
% DYinG = zeros(NX,NY,NG,'single');
% for k = 1:NG
%   DYinG(:,:,k) = DYG;
% end


% starting with the tiny west end
% I'm doing away with these big VH guys; I want everything to be G1, V1,
% etc. I think this is gonna work! 

 
load lon302eG.mat;
neutralD302 = neutralD;
clear neutralD;
load lon30eG.mat;
neutralD30 = neutralD;
clear neutralD;
load lat59sG.mat;
tf= 73;

NVHtot = zeros(tf,NG);
OVRtot = zeros(tf,NG);
REStot = zeros(tf,NG);
VOLtot = zeros(tf,1);


%for t = 1:tf
t = 15;

ycSW = 20; 
ycNW = 80;
   G1 = neutralD302(1,ycSW:ycNW,:,t);
   V1 = squeeze(ncread(strcat(folder, 'bsose_i133b_lon302e_2014_5day_Uvel.nc'), 'UVEL' , [2 ycSW 1 t], [1 ycNW-ycSW+1  Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N1 = squeeze(ncread(strcat(folder, 'bsose_i133b_lon302e_2014_5day_DIC.nc'), 'TRAC01', [2 ycSW 1 t], [1 ycNW-ycSW+1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   Hinz = squeeze(ncread(strcat(folder, 'bsose_i133b_lon302e_2014_5day_Uvel.nc'), 'hFacW', [2 ycSW 1], [1 ycNW-ycSW+1 Inf]));
   DYinG = squeeze(ncread(strcat(folder, 'bsose_i133b_lon302e_2014_5day_Uvel.nc'), 'dyG', [2 ycSW], [1 ycNW-ycSW+1]));
   Hinz = Hinz.*DRF';
   %MAP G,N TO W POINTS??
   %G = (G + G([1 1:NX-1],:,:)).*dnmS; 
   %N = (N + N([1 1:NX-1],:,:)).*dnmS;
   VH = zeros(ycNW- ycSW+1, NG);
   H = VH; NH= VH; NVH = VH;
   for kcur=1:NG
     tmp1=V1.*Hinz;
     tmp2=Hinz;
     tmp3=V1.*Hinz.*N1;
     tmp4=Hinz.*N1;
     tmp1( (G1< mncurref(kcur)) | (G1>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     
     
     VH(:,kcur) = (sum(tmp1,2)); %mass transport in layer  (vh)
     H(:,kcur)  = (sum(tmp2,2)); %layer thickness
     NVH(:,kcur)= (sum(tmp3,2)); %tracer transport in layer
     NH(:,kcur) = (sum(tmp4,2));
   end
   
   VH = VH.*DYinG'; NVH = NVH.*DYinG';
   
   %VH1_l = squeeze(sum(VH(xcW, ycSW:ycN, :), 3));
   %VH3_l = squeeze(sum(VH(xcE, ycSE:ycN, :), 3));
   
   VH1_y = squeeze(sum(VH,1));
   %VH3_y = squeeze(sum( VH(xcE,ycSE:ycN,:),2));
   Area1 = squeeze(sum(H.*DYinG',1));
   %Area3 = squeeze(sum(H(xcE,ycSE:ycN,:).*DYinG(xcE,ycSE:ycN,:),2));
   cont1 = H; % used for contours
   %cont3 = H(xcE, ycSE:ycN,:); 
   
   NVH1 = squeeze(sum( NVH,1));
   %NVH3 = squeeze(sum( NVH(xcE,ycSE:ycN,:),2));
   
   recipArea1 = Area1; recipArea1(recipArea1==0)=inf; recipArea1 = 1./recipArea1;
   %recipArea3 = Area3; recipArea3(recipArea3==0)=inf; recipArea3 = 1./recipArea3
   
   % this guy has to change to V1_y
   V1_y = VH1_y.*(recipArea1); % weighted average v
   % same for him V3 = VH3_y.*(recipArea3);
   NVH1f = NVH;
   N1_y = squeeze(sum(NH.*DYinG',1)).*recipArea1; % weighted average N
save('gra3133.mat', 'NVH1f', 'V1_y')
   %N3_y = squeeze(sum(NH(xcE,ycSE:ycN,:).*DYinG(xcE,ycSE:ycN,:),2)).*recipArea3;
%%   
%%% now we do this again for the other side.    
 % output from neutral D
ycSE = 120; 
ycNE = 261;
lenE = ycNE-ycSE+1;
   G3 = squeeze(neutralD30(1,ycSE:ycNE,:,t));
   V3 = squeeze(ncread(strcat(folder, 'bsose_i133b_lon30e_2014_5day_Uvel.nc'), 'UVEL' , [2 ycSE 1 t], [1 lenE Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N3 = squeeze(ncread(strcat(folder, 'bsose_i133b_lon30e_2014_5day_DIC.nc'), 'TRAC01', [2 ycSE 1 t], [1 lenE Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   Hinz = squeeze(ncread(strcat(folder, 'bsose_i133b_lon30e_2014_5day_Uvel.nc'), 'hFacW', [2 ycSE 1], [1 lenE Inf]));
   DYinG = squeeze(ncread(strcat(folder, 'bsose_i133b_lon30e_2014_5day_Uvel.nc'), 'dyG', [2 ycSE], [1 lenE]));
   Hinz = Hinz.*DRF';
   %MAP G,N TO W POINTS??
   %G = (G + G([1 1:NX-1],:,:)).*dnmS; 
   %N = (N + N([1 1:NX-1],:,:)).*dnmS;
   VH = zeros(lenE, NG);
   H = VH; NH= VH; NVH = VH;
   for kcur=1:NG
     tmp1=V3.*Hinz;
     tmp2=Hinz;
     tmp3=V3.*Hinz.*N3;
     tmp4=Hinz.*N3;
     tmp1( (G3< mncurref(kcur)) | (G3>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     
     
     VH(:,kcur) = (sum(tmp1,2)); %mass transport in layer  (vh)
     H(:,kcur)  = (sum(tmp2,2)); %layer thickness
     NVH(:,kcur)= (sum(tmp3,2)); %tracer transport in layer
     NH(:,kcur) = (sum(tmp4,2));
   end
   
   VH = VH.*DYinG'; NVH = NVH.*DYinG';
   
   %VH1_l = squeeze(sum(VH(xcW, ycSW:ycN, :), 3));
   %VH3_l = squeeze(sum(VH(xcE, ycSE:ycN, :), 3));
   
   %VH1_y = squeeze(sum(VH,1));
   VH3_y = squeeze(sum(VH,1));
   %Area1 = squeeze(sum(H.*DYinG',1));
   Area3 = squeeze(sum(H.*DYinG',1));
   %cont1 = H; % used for contours
   cont3 = H; 
   
   %NVH1 = squeeze(sum( NVH,1));
   NVH3 = squeeze(sum(NVH,1));
   
   %recipArea1 = Area1; recipArea1(recipArea1==0)=inf; recipArea1 = 1./recipArea1;
   recipArea3 = Area3; recipArea3(recipArea3==0)=inf; recipArea3 = 1./recipArea3;
   
   % this guy has to change to V1_y
   %V1_y = VH1_y.*(recipArea1); % weighted average v
   V3_y = VH3_y.*(recipArea3);
   %N1_y = squeeze(sum(NH.*DYinG',1)).*recipArea1; % weighted average N
   N3_y = squeeze(sum(NH.*DYinG',1)).*recipArea3;
NVH3f = NVH;
save('gra3133.mat', 'NVH3f', 'V3_y', '-append'); 
% Now for the horizontal side 2. 
   
%NG = length(mncurref)-1;
%maskC = hFacC; maskC(maskC>0)=1;
%[NX,NY,NZ]=size(hFacS);
%maskS = hFacS; maskS(maskS>0)=1;
%dnmS = maskC(:,:,:) + maskC(:,[1 1:NY-1],:);dnmS(dnmS==0)=inf; dnmS=1./dnmS;
%lon = XC(:,1);lat = YC(1,:);r = squeeze(RC(:));
%Hinz = hFacS*0;
%for k = 1:NZ
%  Hinz(:,:,k) = hFacS(:,:,k)*DRF(k);
%end
%DXinG = zeros(NX,NY,NG,'single');
%for k = 1:NG
%  DXinG(:,:,k) = DXG;
%end

% Starting over with these variables
% xcW and xcE are defined at the beginning

   G2 = squeeze(neutralD([xcW:2160 1:xcE],1,:,t)); % this is the one we like. 

   ycN = 3;
   norbdx = [xcW:2160 1:xcE];
   V = squeeze(ncread(strcat(folder, 'bsose_i133b_lat59s_2014_5day_Vvel.nc'), 'VVEL' , [1 ycN 1 t], [Inf 1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   V2 = V(norbdx,:);
   N2 = squeeze(ncread(strcat(folder, 'bsose_i133b_lat59s_2014_5day_DIC.nc'), 'TRAC01', [1 ycN 1 t], [Inf 1 Inf 1])); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4
  
   N2 = N2(norbdx,:);
   Hinz = squeeze(ncread(strcat(folder, 'bsose_i133b_lat59s_2014_5day_Vvel.nc'), 'hFacS', [1 ycN 1], [Inf 1 Inf]));
   Hinz = Hinz(norbdx,:);
   DXinG = squeeze(ncread(strcat(folder, 'bsose_i133b_lat59s_2014_5day_Vvel.nc'), 'dxG', [1 ycN], [Inf 1]));       
   DXinG = DXinG(norbdx);
   Hinz = Hinz.*DRF';

%VH = zeros(NX,NY,NG);
VH = zeros(xcE+2161-xcW, NG); 
H = VH;
 NVH = VH;
 NH = VH;
 
   %G = ncread('../../../data/bSOSE/iter105/3day/neutrald.nc', 'GAMMA', [1 1 1 10*t-5], [Inf Inf Inf 1]); %rdmds([fpath 'GAMMA'],iter,'rec',t)
   %V =  ncread('../../../data/bSOSE/iter105/monthly/monvvel.nc', 'VVEL', [1 1 1 t] , [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   %N = ncread('../../../data/bSOSE/iter105/monthly/monDIC.nc', 'TRAC01', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   %MAP G,N TO W POINTS
   %G = (G + G(:,[1 1:NY-1],:)).*dnmS; 
   %N = (N + N(:,[1 1:NY-1],:)).*dnmS;
   
   for kcur=1:NG
     tmp1=V2.*Hinz;
     tmp2=Hinz;
     tmp3=V2.*Hinz.*N2;
     tmp4=Hinz.*N2;
     tmp1( (G2< mncurref(kcur)) | (G2>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     VH(:,kcur) = (sum(tmp1,2)); %mass transport in layer  (vh)
     H(:,kcur)  = (sum(tmp2,2)); %layer thickness
     NVH(:,kcur)= (sum(tmp3,2)); %tracer transport in layer
     NH(:,kcur) = (sum(tmp4,2));
   end
   VH = VH.*DXinG; NVH = NVH.*DXinG;
   %%
   VH2_l = squeeze(sum(VH, 2));
   
   VH2_x = squeeze(sum( VH,1));
   
   Area2 = squeeze(sum(H.*DXinG,1));
   cont2 = H;
   NVH2 = squeeze(sum( NVH,1));
   recipArea2 = Area2; recipArea2(recipArea2==0)=inf; recipArea2 = 1./recipArea2;
   V2_x = VH2_x.*recipArea2; % weighted average v
   NVH2f = NVH;
   N2_x = squeeze(sum(NH.*DXinG,1)).*recipArea2; % weighted average N
   save('gra3133.mat', 'NVH2f', 'V2_x','-append') % USE THIS ONEE!! 


NVHtot(t,:) = NVH1 - NVH2 - NVH3;
OVRtot(t,:) = VH1_y.*N1_y - VH2_x.*N2_x - VH3_y.*N3_y;
REStot(t,:) = NVHtot(t,:) - OVRtot(t,:);
VOLtot(t) = sum(VH1_y-VH2_x-VH3_y)*sum(N1_y.*Area1+N2_x.*Area2+N3_y.*Area3)/sum(Area1+Area2+Area3);
%end

%save('2graeme133out.mat', 'VH1_y','VH2_x','VH3_y', 'N1_y', 'N2_x', 'N3_y','Area1','Area2','Area3','NVH1','NVH2','NVH3', 'N1','N2','N3','cont1', 'cont2','cont3')
%save('graeme133out3.mat', 'NVHtot', 'OVRtot','REStot', 'VOLtot')

%% Make the figures
% We need the DIC of course - that's easy. No need to combine the three
% separate sheets. 
% We need a combined layer summed vh and averaged carbon () ... I have that...
% I actually have all this information from before ... 

figure(1)
levs = mncurref(2:end);
TotalVHbyLayer =  1e-6*(-VH1_y + VH2_x + VH3_y);
barh(flipud(TotalVHbyLayer))
yticklabels(fliplr(levs))
xlabel('Volume (Sverdrups)')
ylabel('Neutral Density Upper limits')
title('Volume Export by Layer')
%%
figure(2) 

TotalCbyLayer = (N1_y.*Area1 + N2_x.*Area2 + N3_y.*Area3)./(Area1 + Area2 + Area3);
barh(flipud(TotalCbyLayer))
yticklabels(fliplr(levs))
ylabel('Neutral Density Upper limits')
xlim([2.2 2.35])
title('Carbon concentration by layer')

%%
figure(3)
totlen = length(VH1_l) + length(VH2_l) + length(VH3_l); 
DepthInt = 1e-6*[-VH1_l VH2_l' fliplr(VH3_l)];
toPlot = zeros(1, totlen);
sumV = 0;

for k = 1:1:totlen
    sumV = sumV + DepthInt(k);
    toPlot(k) = sumV;
end

area(toPlot)
%xticks ( [89, 161, length(VH1_l), length(VH1_l) + length(VH2_l), totlen - 161, totlen-89])
%xticklabels( [ 70, 60, latN, latN, 60, 70])

% how do we generalize this? 
% can we make a list of every multiple of 5 between latSW and latN 
% yeah it's floor(latSW/5):-1:ceil(latN/5)
Ticks = [];
Labels = [];
for k = floor(latSW/5):-1:ceil(latN/5)
    [mink, yck] = min(abs(lats+5*k));
    Ticks = [Ticks yck-ycSW];
    Labels = [Labels string(5*k)+' S'];
end
for k = ceil(lonW/5+0.01):1:floor(lonE/5-0.01)
    Ticks = [Ticks length(VH1_l)+15*k-xcW];
    Labels = [Labels string(5*k)+'E'];
end
for k = ceil(latN/5):1:floor(latSE/5)
    [mink, yck] = min(abs(lats+5*k));
    Ticks = [Ticks totlen-yck+ycSE];
    Labels = [Labels string(5*k)+'S'];
end
ylabel('Cumulative Export (Sv)')
xticks(Ticks)
xticklabels(Labels)
hold on 
xline(length(VH1_l));
xline(totlen - length(VH3_l))
hold off
%%

figure (4)

totalC = [N1' N2' flipud(N3)'];
depths = ncread('../../../data/bSOSE/iter105/monthly/monDIC.nc', 'Z', 1 , Inf);
lev = 2.20:0.008:2.40;
contourf(1:1:totlen, depths, totalC,lev, 'LineColor', 'none')
Ticks = [];
Labels = [];
for k = floor(latSW/5):-1:ceil(latN/5)
    [mink, yck] = min(abs(lats+5*k));
    Ticks = [Ticks yck-ycSW];
    Labels = [Labels string(5*k)+' S'];
end
for k = ceil(lonW/5+0.01):1:floor(lonE/5-0.01)
    Ticks = [Ticks length(VH1_l)+15*k-xcW];
    Labels = [Labels string(5*k)+'E'];
end
for k = ceil(latN/5):1:floor(latSE/5)
    [mink, yck] = min(abs(lats+5*k));
    Ticks = [Ticks totlen-yck+ycSE];
    Labels = [Labels string(5*k)+'S'];
end
ylabel('Depth (m)')
xticks(Ticks)
xticklabels(Labels)
hold on 
for k = 1:1:10
   contk = [squeeze(-sum(cont1(:,:,1:k),3)) squeeze(-sum(cont2(:,:,1:k),3))' fliplr(squeeze(-sum(cont3(:,:,1:k),3))) ]
   if (k==2 || k==5) 
       plot(contk, 'Color', 'w', 'LineWidth' ,2) ;
   else
       plot(contk, 'Color', 'w', 'LineWidth', 1);
   end
end


hold off
% I do still want neutral density contours here... I wonder ... 

%%
Overturn = 3.79e2*TotalCbyLayer.*TotalVHbyLayer;
TotalCTrans = 3.79e-4*(-NVH1+NVH2+NVH3);
Residuals = TotalCTrans - Overturn;
figure(5)
barh(flipud(TotalCTrans))
yticklabels(fliplr(levs))
xlabel('Carbon Export (Tg C/yr)');
ylabel('Neutral Density Upper limits')
title('Total Carbon Export from Box')
%%

figure(6)
barh(flipud(Residuals))
yticklabels(fliplr(levs))
title('Residual Carbon Export')
xlabel('Carbon Export (Tg C/yr)');
ylabel('Neutral Density Upper limits')

