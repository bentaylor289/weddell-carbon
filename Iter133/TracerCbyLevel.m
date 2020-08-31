%%% find concentration of nitrate at each level 
%%% for the Iter133 2013-2018 carbon budget 
%%% looping this one for multiple 5-day periods

clear all
%% inputs to change

wkstart = 1;
wkend = 438;

budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_';
% also dependencies on cellwise budget folder, grid file. 
Nlat = 53;
lonW = 290;
lonE = 35; 

WeddellFix = false; % fix for Antarctic peninsula (we want to stop the Western boundary before 78 S).
% Weddell fix also includes hardcoded locations for Weddell. 
ycSW = 277; % only used if Weddell fix is on

zBounds = [1 50 100 150 200 300 500 2000]; % depths
zbounds = [1 8 13 18 21  26 31 40];
%densityBounds = [ 23 27.55 27.7 27.85 28 28.1 28.2 28.27 28.35 28.4 29]; 
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

load('/local/projects/bSOSE_carbon_Ben/Iter133/BlockyNVHBudget1to438.mat', 'isINS');
monstart =wkstart; % lazy coding
monend = wkend;
n1 = monend-monstart +1; % number of iterations.


volume = zeros(size(hFacC));

for k=1:52
 volume(:,:,k) = hFacC(:,:,k).*RAC(:,:)*DRF(k);
end
mncurref= zBounds;
lays = length(mncurref)-1;

LayContN = zeros(n1,lays);
LayVolN = zeros(n1, lays); 

for t = monstart:monend

    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    
conc = squeeze(ncread('/local/data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5daySnapShots_NO3.nc', 'TRAC04',[1 1 1 t_real], [Inf yc-1 Inf 1]));% somehow I don't have the 5days downloaded which is funny!
%load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'adv','mix','dilut', 'surf', 'tend', 'bio', 'res','corr');

%% Calculation is very simple!

isINS2 = [zeros(length(xcFull),77) isINS];
vol1 = volume(xcFull,1:yc-1,:).*repmat(isINS2,[1 1 52]);
%G = squeeze(ncread(strcat(budgetfolder, 'GAMMA.nc'), 'gamma', [1 1 1 t_real], [Inf yc-1 Inf 1]));
%G1 = G(xcFull,:,:); % only the relevant section
%binner = discretize(G1, mncurref); % divides into density classes 
contTMP = double(conc(xcFull,1:yc-1, :)).*vol1;
volTMP = vol1;
% zbounds defined at beginning to work for this. 
for i = 1:lays
LayContN(l,i) = sum(sum(sum(contTMP(:,:,zbounds(i):(zbounds(i+1)-1)))));
LayVolN(l,i) = sum(sum(sum(volTMP(:,:,zbounds(i):(zbounds(i+1)-1)))));
end

% we need to get that section around the peninsula fixed!
% turns out volume transport is also important! 

%Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 
%clear dilut surf bio tend res corr mix
end    
%%
LevContN =LayContN;
LevVolN = LayVolN;
save('NLevFull.mat', 'LevContN','LevVolN','isINS','zbounds')%
