%%% Graeme's figures plot ... for new iter 133
%%% The idea here is to replicate the calculations that Graeme made, so
%%% it's easy to direct compare the results. 
%%% It should be able to take any 3-sided box to the Antarctic Coast.

%%% I'm definitely going to want to be able to loop this... and this code will need a lot of mods ..
%%% To clarify, this one will be able to do multiple timesteps. TBD which variables are saved ...
%%% This needs to be able to cross the prime meridian - that's not hard ... 


%% Inputs here

latN = 61; % latitude of zonal section (S)
latSW = 65; % latitude of stopping meridional section (once land is reached) on W border.
latSE = 71; % " " on E border
lonW = 302; % we'll work on the version that cross the prime meridian later.
lonE = 25;
t1 = 7; % currently reading from monthly data but that's easy to change.
t2 = 8;

densityBounds = [ 23 27.7 28 28.1 28.2 28.27 28.3 28.33 28.36 28.4 28.5 29];
mncurref = densityBounds; 

fold5 = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';

lats = ncread(strcat(fold5, 'GAMMA.nc'), 'YC');
[mini, ycN] = min(abs(lats+latN));
[mini, ycSW] = min(abs(lats+latSW));
[mini, ycSE] = min(abs(lats+latSE));
xcW = 6*lonW+1; % not 100% sure if this is quite what I want ... 
xcE = 6*lonE+1;
if (lonW > lonE)
	myxc = [xcW:2160 1 :xcE-1];
end
if (lonE > lonW)
	myxc = xcW:xcE-1;
end
% for zonal, I want up to but not including xcE  


load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat');
%% Get and consolidate the layer info. 

NG = length(mncurref)-1;
DONexport = zeros(t2-t1+1,NG);

for (t = t1:t2)

% Start with the two meridional sections!
NG = length(mncurref)-1;
maskC = hFacC; maskC(maskC>0)=1;
[NX,NY,NZ]=size(hFacS);
maskS = hFacS; maskS(maskS>0)=1;
%creates a mask so that the averaging with the grid cell to the west
%doesn't get too messy. It's a cool tool. 
dnmS = maskC(:,:,:) + maskC([1 1:NX-1],:,:);dnmS(dnmS==0)=inf; dnmS=1./dnmS;
lon = XC(:,1);lat = YC(1,:);r = squeeze(RC(:));
Hinz = hFacS*0;
for k = 1:NZ
  Hinz(:,:,k) = hFacW(:,:,k)*DRF(k); % using west side here as well.
end
DYinG = zeros(NX,NY,NG,'single');
for k = 1:NG
  DYinG(:,:,k) = DYG;
end

%%% this way of doing things is ridiculously slow computationally; once that becomes an issue,
%%% we need to move back to the other strategy of just pulling out what I need. 
% Let's just do them one at a time.
 VH = zeros(NX,NY,NG);
 H = VH;
 NVH = VH;
 
   G = ncread(strcat(fold5, 'GAMMA.nc'), 'gamma', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'GAMMA'],iter,'rec',t)
   V =  ncread(strcat(fold5,'Uvel.nc'), 'UVEL', [1 1 1 t] , [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N = ncread(strcat(fold5,'DON.nc'), 'TRAC07', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   %MAP G,N TO W POINTS
   G = (G + G([1 1:NX-1],:,:)).*dnmS; 
   N = (N + N([1 1:NX-1],:,:)).*dnmS;
   N1 = squeeze(N(xcW, ycSW:ycN, :));
   N3 = squeeze(N(xcE, ycSE:ycN, :));
   H = VH; NVH = H;
   NH = H;
   for kcur=1:NG
     tmp1=V.*Hinz;
     tmp2=Hinz;
     tmp3=V.*Hinz.*N;
     tmp4=Hinz.*N;
     tmp1( (G< mncurref(kcur)) | (G>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     VH(:,:,kcur) = (sum(tmp1,3)); %mass transport in layer  (vh)
     H(:,:,kcur)  = (sum(tmp2,3)); %layer thickness
     NVH(:,:,kcur)= (sum(tmp3,3)); %tracer transport in layer
     NH(:,:,kcur) = (sum(tmp4,3));
   end
   
   VH = VH.*DYinG; NVH = NVH.*DYinG;
   %%
   VH1_l = squeeze(sum(VH(xcW, ycSW:ycN, :), 3));
   VH3_l = squeeze(sum(VH(xcE, ycSE:ycN, :), 3));
   
   VH1_y = squeeze(sum( VH(xcW,ycSW:ycN,:),2));
   VH3_y = squeeze(sum( VH(xcE,ycSE:ycN,:),2));
   Area1 = squeeze(sum(H(xcW,ycSW:ycN,:).*DYinG(xcW,ycSW:ycN,:),2));
   Area3 = squeeze(sum(H(xcE,ycSE:ycN,:).*DYinG(xcE,ycSE:ycN,:),2));
   cont1 = H(xcW, ycSW:ycN,:);
   cont3 = H(xcE, ycSE:ycN,:); 
   
   NVH1 = squeeze(sum( NVH(xcW,ycSW:ycN,:),2));
   NVH3 = squeeze(sum( NVH(xcE,ycSE:ycN,:),2));
   
   recipArea1 = Area1; recipArea1(recipArea1==0)=inf; recipArea1 = 1./recipArea1;
   recipArea3 = Area3; recipArea3(recipArea3==0)=inf; recipArea3 = 1./recipArea3;
   
   V1 = VH1_y.*(recipArea1); % weighted average v
   V3 = VH3_y.*(recipArea3);
   N1_y = squeeze(sum(NH(xcW,ycSW:ycN,:).*DYinG(xcW,ycSW:ycN,:),2)).*recipArea1; % weighted average N
   N3_y = squeeze(sum(NH(xcE,ycSE:ycN,:).*DYinG(xcE,ycSE:ycN,:),2)).*recipArea3;
   
% Now for the horizontal side 2. 
   
NG = length(mncurref)-1;
maskC = hFacC; maskC(maskC>0)=1;
[NX,NY,NZ]=size(hFacS);
maskS = hFacS; maskS(maskS>0)=1;
dnmS = maskC(:,:,:) + maskC(:,[1 1:NY-1],:);dnmS(dnmS==0)=inf; dnmS=1./dnmS;
lon = XC(:,1);lat = YC(1,:);r = squeeze(RC(:));
Hinz = hFacS*0;
for k = 1:NZ
  Hinz(:,:,k) = hFacS(:,:,k)*DRF(k);
end
DXinG = zeros(NX,NY,NG,'single');
for k = 1:NG
  DXinG(:,:,k) = DXG;
end


% Starting over with these variables
 VH = zeros(NX,NY,NG);
 H = VH;
 NVH = VH;
 NH = VH;
 
   G = ncread(strcat(fold5,'GAMMA.nc'), 'gamma', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'GAMMA'],iter,'rec',t)
   V =  ncread(strcat(fold5,'Vvel.nc'), 'VVEL', [1 1 1 t] , [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_state'],ts(t),'rec',4);
   N = ncread(strcat(fold5,'DON.nc'), 'TRAC07', [1 1 1 t], [Inf Inf Inf 1]); %rdmds([fpath 'OUTPUT_3dy/diag_bgc'],ts(t),'rec',4);
   %MAP G,N TO W POINTS
   G = (G + G(:,[1 1:NY-1],:)).*dnmS; 
   N = (N + N(:,[1 1:NY-1],:)).*dnmS;
   N2 = squeeze(N(xcW:xcE-1, ycN, :));
   
   for kcur=1:NG
     tmp1=V.*Hinz;
     tmp2=Hinz;
     tmp3=V.*Hinz.*N;
     tmp4=Hinz.*N;
     tmp1( (G< mncurref(kcur)) | (G>=mncurref(kcur+1)) )=0;
     tmp2(tmp1==0)=0;
     tmp3(tmp1==0)=0;
     tmp4(tmp1==0)=0;
     VH(:,:,kcur) = (sum(tmp1,3)); %mass transport in layer  (vh)
     H(:,:,kcur)  = (sum(tmp2,3)); %layer thickness
     NVH(:,:,kcur)= (sum(tmp3,3)); %tracer transport in layer
     NH(:,:,kcur) = (sum(tmp4,3));
   end
   
   
   VH = VH.*DXinG; NVH = NVH.*DXinG;
   %%
   VH2_l = squeeze(sum(VH(myxc, ycN, :), 3));
   
   VH2_x = squeeze(sum( VH(myxc, ycN, :),1));
   
   Area2 = squeeze(sum(H(myxc, ycN, :).*DYinG(myxc, ycN, :),1));
   cont2 = H(xcW:xcE-1, ycN, :);
   NVH2 = squeeze(sum( NVH(myxc, ycN, :),1));
   recipArea2 = Area2; recipArea2(recipArea2==0)=inf; recipArea2 = 1./recipArea2;
   V2 = VH2_x./(Area2); % weighted average v
   N2_x = squeeze(sum(NH(myxc, ycN, :).*DYinG(myxc, ycN, :),1))./Area2; % weighted average N

   DONexport(t-t1+1,NG) = (-NVH1+NVH2+NVH3);

   end

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






