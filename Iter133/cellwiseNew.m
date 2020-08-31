%%% B-SOSE Cellwise DIC budget for Iter 133 (whole SO) 
%%% Ariane Verdy, May 2018
%%% edited Ben Taylor, Apr 15 2020


%%% Changing filenames for runs other than Iter133 may be necessary (but changing the folder names should take care of most)
%%% Likewise for other tracers - but again, there shouldn't be too many errrors!
%%% Other important info: if run for multiple months, this script will output
%%% separate budget files for each month to the folder JRFWcarbBudmonthly (see save function near the end).


clear all
close all

% file name
%diag_budget_file = 'diag_dic_budget';
budgetfolder = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_5day_'; % includes filename prefix

budgetfolder2 = '../../../data/bSOSE/iter133NEW/5day/bsose_i133_2013to2018_'; % includes filename prefix
wkstart = 1; % starting "week"
wkend = 73; % ending "week"


% select area
x = 1:2160; % all longitudes
ylen = 570;
y = 1:ylen; % everywhere poleward of pretty far north; almost everywhere
% don't let y= 588

z = 1:52; % top 1000 m

monstart = wkstart; monend = wkend; % old variablename
%right now I only want one time step, at month 5

%% CODE STARTS HERE



load /local/projects/bSOSE_carbon_Ben/Iter129/grid.mat 
[nx,ny,nz] = size(hFacC);
dz = permute(repmat(DRF(z),[1,nx,ylen]),[2,3,1]).*hFacC(x,y,z); % what is this?


%save ('grid.mat' ,'hFacS', 'hFacC', 'hFacW', 'XC', 'YC', 'RC', 'RF', 'DRC','RAC','DRF','DXG','DYG')
%%
% cell volume, face areas (for flux calculations)

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


for t=1:monend-monstart+1
% initialize
nt=1;
surf = single(zeros(nx,ylen,nz,nt));
dilut = single(zeros(nx,ylen,nz,nt));
bio = single(zeros(nx,ylen,nz,nt));
mix = single(zeros(nx,ylen,nz,nt));
adv = single(zeros(nx,ylen,nz,nt));
adv_h = single(zeros(nx,ylen,nz,nt));
adv_v = single(zeros(nx,ylen,nz,nt));
corr = single(zeros(nx,ylen,nz,nt));
div = single(zeros(nx,ylen,nz,nt));
tend = single(zeros(nx,ylen,nz,nt));


	%%
% read diagnostics
% calculate tendencies in mol/m3/s


% tendency due to air-sea flux 
% diagnostic: BLGCFLX
% surface flux (mol/m3/s)

% flux = rdmds('diag_airsea',ts(t),'rec',3); 
flux = ncread(strcat(budgetfolder,'surfCO2flx.nc'), 'BLGCFLX', [1 1 t+monstart-1], [Inf ylen 1]);
surf(:,:,1,nt) = flux(x,y)./(DRF(1)*squeeze(hFacC(x,y,1)));
%%
% tendency due to dilution
% diagnostic: FORCTR01
% forcing tendency (mol/m3/s) includes effects of E-P-R and sponge layer contributions

dilut(:,:,1,nt) = ncread(strcat(budgetfolder,'ForcDIC.nc'), 'ForcTr01', [1 1 1 t+monstart-1], [Inf ylen 1 1]);

% end
bio(:,:,:,nt) = ncread(strcat(budgetfolder,'BLGBIOC.nc'), 'BLGBIOC', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);

%%
% advection
% diagnostics: ADVxTr01, ADVyTr01, ADVrTr01
% advective flux (mol/s)
% advflux = rdmds(diag_budget_file,ts(t),'rec',1:3);
advflux_x(x,:,:) = ncread(strcat(budgetfolder,'ADVx_DIC.nc'), 'ADVxTr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
advflux_y(:,:,:) = ncread(strcat(budgetfolder,'ADVy_DIC.nc'), 'ADVyTr01', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
advflux_z(:,:,z) = ncread(strcat(budgetfolder,'ADVr_DIC.nc'), 'ADVrTr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
if x(end)==nx
 advflux_x = advflux_x([x 1],y,z);
else
 advflux_x = advflux_x([x x(end)+1],y,z);
end

% if y(end)==ny
%     advflux_y = advflux_y(x,[y y(end)],z);
% else 
% advflux_y = advflux_y(x,[y y(end)+1],z);
% end
if z(end)==nz
 advflux_z = advflux_z(x,y,[z z(end)]); 
 advflux_z(:,:,end) = 0*advflux_z(:,:,end);
else
 advflux_z = advflux_z(x,y,[z z(end)+1]); 
end
adv_x = diff(advflux_x,1,1)./volume;
adv_y = diff(advflux_y,1,2)./volume;
adv_z = -diff(advflux_z,1,3)./volume;

% minus sign because advective flux is on lhs, moving to rhs
adv(:,:,:,nt) = -(adv_x + adv_y + adv_z);

 

% % THIS STUFF IS WRONG. tryna write a divergence term using the UTRAC01 terms.
% utrans = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_UTRAC01.nc'), 'UTRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% vtrans = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_VTRAC01.nc'), 'VTRAC01', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
% wtrans = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_WTRAC01.nc'), 'WTRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% 
% % calculate div from diff's
% div(1:2159, :, :) = div(1:2159, :, :) + diff(utrans,1);
% div(2160,:, : ) = div(2160, : ,:) + (utrans(2160,:,:) - utrans(1,:,:));
% div(:,:,:) = div(:,:,:) + diff(vtrans,1,2);
% div(:,:,1:51) = div(:,:,1:51) + diff(wtrans,1, 3);
% div(:,:,52) = div(:,:,52) + wtrans(:,:,52);

% divergence, correctly written. 
% diagnostics: UVEL, VVEL, WVEL
%vel = rdmds(diag_state_file,ts(t),'rec',3:5);
% vel = zeros(nx,ylen+1,nz,3);
% 
% vel(:,1:ylen,:,1) = ncread(strcat(budgetfolder,'../bsose_i129_2013to2018_monthly_Uvel.nc'), 'UVEL', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% vel(:,:,:,2) = ncread(strcat(budgetfolder,'../bsose_i129_2013to2018_monthly_Vvel.nc'), 'VVEL', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
% vel(:,1:ylen,:,3) = ncread(strcat(budgetfolder,'../bsose_i129_2013to2018_monthly_Wvel.nc'), 'WVEL', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% if x(end)==nx
%  U = vel([x 1],y,z,1).*areaWest;
% else
%  U = vel([x x(end)+1],y,z,1).*areaWest;
% end
% V= vel(x,[y ylen+1],z,2).*areaSouth;
% 
% if z(end)==nz
%  W = vel(x,y,[z z(end)],3).*areaTop; 
%  W(:,:,end) = 0*W(:,:,end);
% else
%  W = vel(x,y,[z z(end)+1],3).*areaTop; 
% end
% div_x=diff(U,1,1)./volume;
% div_y=diff(V,1,2)./volume;
% div_z=-diff(W,1,3)./volume;
% 
% % tracer field (mol/m3)
% tracer = ncread(strcat(budgetfolder,'../bsose_i122_2013to2017_monthly_DIC.nc'), 'TRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% 
% % tmp = rdmds('diag_bgc',ts(t),'rec',1);
% % tracer = tmp(x,y,z);
% div(:,:,:,t) = tracer.*(div_x + div_y + div_z);
% 
% % advection components
% adv_h(:,:,:,t) = -(adv_x + adv_y)+tracer.*(div_x + div_y);
% adv_v(:,:,:,t) = -(adv_z)+tracer.*(div_z);

% correction to vertical advection at z=0
% diagnostic: WTRAC01
tmp = ncread(strcat(budgetfolder,'WVELDIC.nc'), 'WTRAC01', [1 1 1 t+monstart-1], [Inf ylen 1 1]);
corr(:,:,1,nt) = tmp(x,y,1)./dz(x,y,1); % we may need to divide by the area as well!
%%
% mixing
% diagnostics: DFxETr01, DFyETr01, DFrITr01
% diffusive flux (mol/s)
%diffflux = rdmds(diag_budget_file,ts(t),'rec',4:6);
diffflux = zeros(nx, ylen+1, nz, 3);
diffflux(:,1:ylen,:,1) = ncread(strcat(budgetfolder,'DFxE_DIC.nc'), 'DFxETr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
diffflux(:,:,:,2) = ncread(strcat(budgetfolder,'DFyE_DIC.nc'), 'DFyETr01', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
diffflux(:,1:ylen,:,3) = ncread(strcat(budgetfolder,'DFrI_DIC.nc'), 'DFrITr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
if x(end)==nx
 diffflux_x = diffflux([x 1],y,z,1);
else
 diffflux_x = diffflux([x x(end)+1],y,z,1);
end
diffflux_y = diffflux(x,[y ylen+1],z,2);
if z(end)==nz
 diffflux_z = diffflux(x,y,[z z(end)],3); 
 diffflux_z(:,:,end) = 0*diffflux_z(:,:,end);
else
 diffflux_z = diffflux(x,y,[z z(end)+1],3); 
end
mix_x = diff(diffflux_x,1,1)./volume;
mix_y = diff(diffflux_y,1,2)./volume;
mix_z = -diff(diffflux_z,1,3)./volume;

% minus sign because diffusive flux is on lhs, moving to rhs
mix(:,:,:,nt) = -(mix_x + mix_y + mix_z);
%%
% total tendency
dt = 3600*24*5; % number of seconds in a 5-day period
snap = ncread(strcat(budgetfolder2,'5daySnapShots_DIC.nc'), 'TRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 2]);
tend(:,:,:,nt) = diff(snap(x,y,z,:),1,4)/dt;




clear tmp snap flux diffflux* div_* vel U V W advflux* %mix_* \%



% fix divisions by hFacC=0
surf(isnan(surf)) = 0;
adv(isnan(adv)) = 0;
adv_h(isnan(adv_h)) = 0;
adv_v(isnan(adv_v)) = 0;
mix(isnan(mix)) = 0;
corr(isnan(corr)) = 0;

%%
% remove correction from advection
%adv = adv-corr;
%adv_v = adv_v-corr;
%%

% residual
res = adv+div+mix+bio+surf+dilut-tend-corr;
save('cellwise5day/carbBudgetWk' + string(t+ monstart-1) +'.mat', 'adv', 'mix', 'bio', 'surf', 'dilut', 'tend', 'res', 'corr');
end % for t
%adv=adv+div;
%clear div
%save('carbBudgetmon' + string(monstart) + 'thru'+string(monend)+'.mat', 'adv', 'adv_h', 'adv_v', 'adv_z', 'adv_y', 'adv_x', 'mix', 'mix_x', 'mix_y', 'mix_z', 'bio', 'surf', 'dilut', 'tend', 'res','div');
%What would I do if I got this to work, tho?? This is wild ...


%% 
%save('carbBudgetJRFWmon'+string(monstart)' + 'thru' + string(monend)+'.mat', 'adv', 'mix', 'tend','res', 'dilut', 'surf', 'bio');
% 
% check that the terms balance locally
% plot a single time, single location
%%
t=1; x1=520; y1=470;
%load('grid.mat', 'RC');
figure(5);
zfig = 1:1:36;
RC = RC(zfig);
hold on
plot(squeeze(tend(x1,y1,zfig,t)),RC); hold on
plot(squeeze(surf(x1,y1,zfig,t)),RC);
plot(squeeze(dilut(x1,y1,zfig,t)),RC);
plot(squeeze(mix(x1,y1,zfig,t)),RC);
plot(squeeze(bio(x1,y1,zfig,t)),RC);
plot(squeeze(adv(x1,y1,zfig,t)),RC);
plot(squeeze(res(x1,y1,zfig,t)),RC,'--k');
ylabel('depth'); xlabel('mol/m3/s');
xlim([-5e-9 5e-9])
title('Depth profile at (x1,y1) = '+string(x1) + ', ' +string(y1))
legend('tend','surf','dilut','mix','bio','adv','res') % tend and res removed
saveas(gcf, 'Cellwisetestfig.png')


%%
% title('Carbon Profile at Lat = 100 W, 45 S')
% 
% %% residual
% figure(111)
% dep=1;
% xlon=1830;
% plot(res(xlon, :, dep))
% hold on 
% plot(adv(xlon,: , dep))
% plot(tend(xlon,: , dep)*1)
% plot(surf(xlon,: , dep)*1)
% plot(dilut(xlon,: , dep)*1)
% plot(mix(xlon,:,dep)*1)
% %%
% legend('res', 'adv', 'tend','surf','dilut','mix')
% %%
% title('Surface terms thru Rio outflow, unscaled')
% 
