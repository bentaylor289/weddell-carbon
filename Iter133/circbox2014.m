%%% cirumpolar box budget for 2014 in iter 133.


cy = 334;%latitude 59 S
advy=3% make sure the advection lines up with this guy! 4 with 335, 2 w 333, etc. 
ts = 74;
tl = 5; % length of time

dt = 432000; % 5day timestep

folder = '../../../data/bSOSE/iter133OLD/';

x = 1:2160; % all longitudes
ylen = 570;
y = 1:ylen; % everywhere poleward of pretty far north; almost everywhere
% don't let y= 588
z = 1:52;
load grid.mat 
[nx,ny,nz] = size(hFacC);
dz = permute(repmat(DRF(z),[1,nx,ylen]),[2,3,1]).*hFacC(x,y,z); 
volume = zeros(nx,ylen,nz);
areaWest = zeros(nx+1,ylen,nz);
areaSouth = zeros(nx,ylen+1,nz);
areaTop = zeros(nx,ylen,nz+1);
for k=1:nz
 volume(:,y,k) = hFacC(x,y,k).*RAC(x,y)*DRF(k);
 areaTop(:,y,k) = RAC(x,y);
 if x(end)==nx
  areaWest(:,:,k)  = DYG([x 1],y).*DRF(k).*hFacW([x 1],y,k);
 else
  areaWest(:,:,k)  = DYG([x x(end)+1],y).*DRF(k).*hFacW([x x(end)+1],y,k);
 end
 areaSouth(:,:,k) = DXG(x,[y ylen+1]).*DRF(k).*hFacS(x,[y ylen+1],k); 
end
areaTop(:,:,nz+1) = RAC(x,y);
area = RAC(x,y);

nt = tl;
surf = single(zeros(1,nt));
dilut = single(zeros(1,nt));
bio = single(zeros(1,nt));
mix = single(zeros(1,nt));
adv = single(zeros(1,nt));
%adv_h = single(zeros(nx,ylen,nz,nt));
%adv_v = single(zeros(nx,ylen,nz,nt));
corr = single(zeros(1,nt));
%div = single(zeros(nx,ylen,nz,nt));
tend = single(zeros(1,nt));

for t = ts:tl+ts-1
    tr= t-ts+1; %for indexing the results
%%
airseat = ncread(strcat(folder,'bsose_i133_2013to2018_5day_surfCO2flx.nc'), 'BLGCFLX', [1 1 t], [Inf cy 1]);
airseainfo = ncinfo(strcat(folder,'bsose_i133_2013to2018_5day_surfCO2flx.nc'));
airseatimes = ncread(strcat(folder,'bsose_i133_2013to2018_5day_surfCO2flx.nc'), 'time');
surf(tr) = sum(sum(airseat.*volume(:,1:cy,1)));

%%
wvel2014 = ncread(strcat(folder, 'bsose_i133b_surfPlwrd59s_2014_5day_WVELDIC.nc'), 'WTRAC01',[1 1 1 tr], [Inf cy 1 1]);
wvels = squeeze(wvel2014);
corr(tr) = sum(sum(double(wvels).*volume(:,1:cy,1)/4.2));

%%
biot = ncread(strcat(folder,'bsose_i133_2013to2018_5day_BLGBIOC.nc'), 'BLGBIOC', [1 1 1 t], [Inf cy Inf 1]);
bio(tr) = sum(sum(sum(biot.*volume(:,1:cy,:))));

%%
snap = ncread(strcat(folder,'bsose_i133_2013to2018_5daySnapShots_DIC.nc'), 'TRAC01', [1 1 1 t], [Inf cy Inf 2]);
tend(tr) = sum(sum(sum(squeeze(diff(snap,1,4)/dt).*volume(:,1:cy,:))));
%%
snap = ncread(strcat(folder,'bsose_i133_2013to2018_5daySnapShots_DIC.nc'), 'time');

%%
advy = ncread(strcat(folder,'sections/bsose_i133b_lat59s_2014_5day_ADVy_DIC.nc'), 'YG');
advt = ncread(strcat(folder,'sections/bsose_i133b_lat59s_2014_5day_ADVy_DIC.nc'), 'ADVyTr01',[1 3 1 tr],[Inf 1 Inf 1]);
adv(tr) = -sum(sum(squeeze(advt)));

%%
dift = ncread(strcat(folder,'sections/bsose_i133b_lat59s_2014_5day_DFyE_DIC.nc'), 'DFyETr01',[1 3 1 tr],[Inf 1 Inf 1]);
mix(tr) = -sum(sum(squeeze(dift)));

%%
forct = ncread(strcat(folder,'bsose_i133_2013to2018_5day_ForcDIC.nc'), 'ForcTr01', [1 1 1 t], [Inf cy 1 1]);
forct = squeeze(forct);
dilut(tr) = sum(sum(forct.*volume(:,1:cy,1)));
end
%%
residual = tend + corr - adv - dilut - bio -surf - mix;
%%

figure(3)
hold on 
plot(residual)
plot(adv)
plot(surf)
plot(bio)
plot(surf+residual)
legend('res', 'adv', 'surf', 'bio', '???') 

