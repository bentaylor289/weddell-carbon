%%% Difficult terms of the Weddelli carbon budget 
%%% for the Iter133 2013-2018 carbon budget 
%%% testing assumptions and claims made by the tracer_budget document. 
%%%
%%% looping this one for multiple 5-day periods

clear all
load('../Iter129/grid.mat', 'hFacC', 'RAC', 'DRF');
surf = squeeze(zeros(size(hFacC(:,:,1))));
%%
toph = DRF(1);
%for k=1:52
 surf = hFacC(:,:,1).*RAC(:,:);
%end
%topA = squeeze(volume(:,:,1)/toph);
wkstart = 1;
wkend = 438;
monstart =wkstart; % lazy coding
monend = wkend;

file = '../../../data/bSOSE/iter133NEW/monthly/bsose_i133_2013to2018_monthly_SeaIceArea.nc';
%snaps = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5daySnap';
Nlat = 61;
lats = ncread('../../../data/bSOSE/iter122/monthly/bsose_i122_2013to2017_monthly_DIC.nc', 'YC');
[min1,yc] = min(abs(lats+Nlat));
lonW = 302;
lonE = 25;
xc1 = 6*lonW+1;
xc2 = 6*lonE+1;
xcFull = [xc1:2160 1:xc2-1];

WeddellFix = true
ycSW = 277;

n1 = monend-monstart +1; % number of iterations.

%toptend = zeros(n1,1);
%topcont = zeros(n1,1);
%partetaC = zeros(n1,1);
%partCeta = zeros(n1,1);
%SIarea = zeros(n1,1);
%adveta = zeros(n1,1);
%mixeta = zeros(n1,1);
%tendTot = zeros(n1,1);
%VTrans = zeros(n1,1);
%Cres = zeros(n1,1);

%for t = monstart:monend
t= 1
    t_array = 1;
    t_real = t;
    l = t- monstart+1;
    
SIa = ncread(file, 'SIarea'); 
%load('cellwise5day/carbBudgetWk'+string(t)+'.mat', 'tend', 'bio');

%tend(isnan(tend)) = 0;
%dilut(isnan(dilut)) = 0;
%surf(isnan(surf)) = 0;
%bio(isnan(bio)) = 0;
%resTot = sum(sum(sum(double(res(xcFull,1:yc-1, :, t_array)).*volume(xcFull,1:yc-1,:))));
totfrac = sum(sum(surf(xcFull,1:yc-1).*SIa(xcFull,1:yc-1,:),2),1)/sum(sum(surf(xcFull,1:yc-1))) ;

%divTot(l) = sum(sum(sum(double(div(:,1:yc-1, :, t_array)).*volume(:,1:yc-1,:))));
%%
% turns out volume transport is also important! 

%Cres(l) = tendTot(l)-dilutTot(l)-surfTot(l)-advTot(l)-bioTot(l)+resTot(l) - mixTot(l)+ corrTot(l); 

%%


%%
save('Weddell'+string(Nlat)+'Wk'+string(monstart)+'to'+string(monend)+'SIfrac.mat', 'totfrac');
