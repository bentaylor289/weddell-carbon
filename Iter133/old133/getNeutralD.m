folder = '../../../data/bSOSE/iter133OLD/sections/bsose_i133b_';

% Tinf = ncinfo(strcat(folder, 'lat59s_2014_5day_THETA.nc'))
 T = ncread(strcat(folder, 'lat59s_2014_5day_THETA.nc'), 'THETA', [1 3 1 1], [Inf 1 Inf Inf]);
 S = ncread(strcat(folder ,'lat59s_2014_5day_SALT.nc'), 'SALT', [1 3 1 1], [Inf 1 Inf Inf]);

% % Tinf = ncinfo(strcat(folder, 'lat59s_2014_5day_THETA.nc'))
%  T = ncread(strcat(folder, 'lon302e_2014_5day_THETA.nc'), 'THETA', [2 1 1 1], [1 Inf Inf Inf]);
%  S = ncread(strcat(folder ,'lon302e_2014_5day_SALT.nc'), 'SALT', [2 1 1 1], [1 Inf Inf Inf]);
%  
%  T = ncread(strcat(folder, 'lon30e_2014_5day_THETA.nc'), 'THETA', [2 1 1 1], [1 Inf Inf Inf]);
% S = ncread(strcat(folder ,'lon30e_2014_5day_SALT.nc'), 'SALT', [2 1 1 1], [1 Inf Inf Inf]);
% check that these have the right units.

% Poorly translated to MATLAB by Ben Taylor 2020
% def calc_sigmantr(ptem,psal):
%     '''function sigmantr = calc_sigmantr(ptem,psal)
%         CALCULATE APPROXIMATE NEUTRAL DENSITY FROM POTENTIAL TEMPERATURE AND
%         SALINITY, using functional form of MacDougall and Jackett (2005)
%         INPUT
%            ptem    potential temperature in degrees celsius
%            psal    salinity in psu
%         OUTPUT
%            sigmantr    approximation of neutral density from empirical formula of
%                         MacDougall and Jackett (2005)
%         Copied across from cdfsigntr.f90 from CDFTOOLS package by G.MacGilchrist
%         (gmacgilchrist@gmail.com)'''

    dl_t = T;
    dl_s = S;
    dl_sr= sqrt(abs(dl_s));

    % Numerator
    % T-polynome
    dl_r1=((-4.3159255086706703E-4*dl_t+8.1157118782170051E-2 ) .*dl_t+2.2280832068441331E-1) .*dl_t + 1002.3063688892480E0;
    % S-T Polynome
    dl_r2=(-1.7052298331414675E-7*dl_s-3.1710675488863952E-3*dl_t-1.0304537539692924E-4).*dl_s;
    
    % Denominator
    % T-Polynome
    dl_r3=(((-2.3850178558212048E-9*dl_t-1.6212552470310961E-7).*dl_t+7.8717799560577725E-5).*dl_t+4.3907692647825900E-5).*dl_t+1.0e0;
    % S-T Polynome
    dl_r4=((-2.2744455733317707E-9*dl_t.*dl_t+6.0399864718597388E-6).*dl_t-5.1268124398160734E-4 ).*dl_s;
    % S-T Polynome
    dl_r5=(-1.3409379420216683E-9*dl_t.*dl_t-3.6138532339703262E-5).*dl_s.*dl_sr;
%%
    %%% Neutral density
    neutralD = ( dl_r1 + dl_r2 )./( dl_r3 + dl_r4 + dl_r5 ) - 1000e0;
    
save('lat59sG.mat', 'neutralD');
    
    
