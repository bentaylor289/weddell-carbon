%%% Weddell process maps by year and layer 
%%% for the Iter133 2013-2018 carbon budget 
%%% looping this one for multiple 5-day periods

clear all
load('/local/projects/bSOSE_carbon_Ben/Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
volume = zeros(size(hFacC));
%%

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
wkstart = 1;
wkend = 438;
monstart =wkstart; % lazy coding
monend = wkend;

budgetfolder = '/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';

Nlat = 60;
lats = ncread('/local/data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 300;
lonE = 30;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcFull = [xc1:2160 1:xc2-1];
xlen = (lonE+360-lonW)*6;
%WeddellFix = true
%ycSW = 277;

n1 = monend-monstart +1; % number of iterations.

densityBounds = [23 27.7 28 28.1 28.2 28.27 28.31 28.35 28.37 28.4 29];
%densityBounds = [ 23 27.55 27.7 27.85 28 28.1 28.2 28.27 28.35 28.4 29]; 
mncurref= densityBounds;
lays = length(mncurref)-1;

%advTot = zeros(n1,lays);
mixMap = zeros(xlen, yc-1,n1,lays);
%resTot = zeros(n1,lays);
%surfTot = zeros(n1,lays);
bioMap = zeros(xlen, yc-1, n1,lays);
%dilutTot = zeros(n1,lays);
%corrTot = zeros(n1,lays);
%tendTot = zeros(n1,lays);
%VTrans = zeros(n1,lays);
%Cres = zeros(n1,lays);
reminMap = zeros(xlen, yc-1,n1, lays); 

for t = monstart:monend

    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    

load('/local/projects/bSOSE_carbon_Ben/Iter133/cellwise5day/carbBudgetWk'+string(t)+'.mat', 'mix');%,'dilut', 'surf', 'tend', 'bio', 'res','corr');

%tend(isnan(tend)) = 0;
mix(isnan(mix)) = 0;
%dilut(isnan(dilut)) = 0;
%surf(isnan(surf)) = 0;
%bio(isnan(bio)) = 0;
%res(isnan(res)) = 0;
%corr(isnan(corr)) = 0;
%div(isnan(div)) = 0;

bio = squeeze(ncread(strcat(budgetfolder, 'NPP.nc'), 'BLGNPP', [ 1 1 1 t_real], [Inf yc-1 Inf 1])); 
remin = squeeze(ncread(strcat(budgetfolder, 'NCP.nc'), 'BLGNCP', [ 1 1 1 t_real], [Inf yc-1 Inf 1]));
remin = bio-remin; 
G = squeeze(ncread(strcat(budgetfolder, 'GAMMA.nc'), 'gamma', [1 1 1 t_real], [Inf yc-1 Inf 1]));
G1 = G(xcFull,:,:);
binner = discretize(G1, mncurref);
%dilutTMP = double(dilut(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1);
%surfTMP = double(surf(xcFull,1:yc-1, 1, t_array)).*volume(xcFull,1:yc-1,1);
%corrTMP = double(corr(xcFull, 1:yc-1, 1, t_array)).*volume(xcFull, 1:yc-1, 1);
%tendTMP = double(tend(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:);
mixcon = double(mix(xcFull,1:yc-1, :)).*volume(xcFull,1:yc-1,:);
biocon = double(bio(xcFull,1:yc-1, :)).*volume(xcFull,1:yc-1,:);
remincon = double(remin(xcFull, 1:yc-1, :)).*volume(xcFull,1:yc-1,:);

for j = 1:xlen
	for k = 1:(yc-1)
bioTMP = squeeze(biocon(j,k,:));
mixTMP = squeeze(mixcon(j,k,:));
reminTMP = squeeze(remincon(j,k,:));
binnerTMP = squeeze(binner(j,k,:));
%dilutTot(l,i) = sum(sum(dilutTMP(binner(:,:,1)==i)));
%surfTot(l,i) = sum(sum(surfTMP(binner(:,:,1)==i)));
%corrTot(l,i) = sum(sum(corrTMP(binner(:,:,1)==i)));
%tendTot(l,i) = sum(sum(sum(tendTMP(binner==i))));

	for i=1:lays 
	bioMap(j,k,l,i) = sum(bioTMP(binnerTMP==i));
	mixMap(j,k,l,i) = sum(mixTMP(binnerTMP==i));
	reminMap(j,k,l,i) = sum(reminTMP(binnerTMP==i));
end
end
end

clear bio remin  mix
end   

% should probably fix this to group actual winters together... not too hard I think! 
year=1; 
bioOut = zeros(xlen, yc-1, ceil(n1/73), lays);
mixOut = zeros(xlen, yc-1, ceil(n1/73), lays);
reminOut = zeros(xlen, yc-1, ceil(n1/73), lays);

while (year*73 < n1)
	bioOut(:,:,year,:) = mean(bioMap(:,:,year*73-72: year*73,:),3);
	reminOut(:,:,year,:) = mean(reminMap(:,:,year*73-72: year*73,:),3);
	mixOut(:,:,year,:) = mean(mixMap(:,:,year*73-72: year*73,:),3);
year = year+1;
end

bioOut(:,:,year,:) = mean(bioMap(:,:,year*73-72:end,:),3);
reminOut(:,:,year,:) = mean(reminMap(:,:,year*73-72:end,:),3);
mixOut(:,:,year,:) = mean(mixMap(:,:,year*73-72:end,:),3);

%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

%Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 


cv = 3.7843e-4; % units of Tg/yr
%Cres = Cres * cv;
%tendTot = tendTot * cv; 
%resTot = resTot *cv;
%advTot=  advTot*cv;
%surfTot = surfTot*cv;
bioOut = bioOut*cv;
reminOut= reminOut*cv;
mixOut = mixOut*cv;
%resTot = (resTot+Cres)*cv;
%dilutTot = dilutTot*cv;
%corrTot= corrTot';

%%



%%
save('Processmaps'+string(Nlat)+'budgetWk'+string(monstart)+'to'+string(monend)+'Juldens.mat', 'mncurref', 'reminOut','bioOut', 'mixOut');
