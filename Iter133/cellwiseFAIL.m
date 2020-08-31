%%% Check cellwise budget works for a random column
% Think I'll need more than I have for that... 


%%% Using data from year 2014, old iter 133. 

cx = 2;
cy = 336;%latitude 59 S
ts = 74;
tl = 73; % length of time

folder = '../../../data/bSOSE/iter133OLD/';

%%
airsea = ncread(strcat(folder,'bsose_i133_2013to2018_5day_surfCO2flx.nc'), 'BLGCFLX', [cx cy ts], [1 1 tl]);
airseainfo = ncinfo(strcat(folder,'bsose_i133_2013to2018_5day_surfCO2flx.nc'));
airseatimes = ncread(strcat(folder,'bsose_i133_2013to2018_5day_surfCO2flx.nc'), 'time');

%%
wvel2014 = ncread(strcat(folder, 'bsose_i133b_surfPlwrd59s_2014_5day_WVELDIC.nc'), 'WTRAC01');
wvels = squeeze(wvel2014(:,:,1,:));

%%
bio = ncread(strcat(folder,'bsose_i133_2013to2018_5day_BLGBIOC.nc'), 'BLGBIOC', [cx cy 1 ts], [1 1 Inf tl]);

%% disappointing that I failed to get the advective terms to balance the cellwise budget. But that's OK!




