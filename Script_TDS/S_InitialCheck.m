%% Cex version
% Program to perform TDS Analysis
% Sergio Revuelta v1.0
% 

% Cleaning everything before starts
clc, clf, close all, clear all
%%
f_namedata = 'pType';
% Constants
c_c = 299792458; % Light velocity
% Variables
lv_f = 0;
% v_dsu = 1.06E-3; % Substrate thickness
v_dsas = 0.3e-3; % Sample thickness. In this case will be the sample thickness between both subtrates
v_dsar = v_dsas; % Sample thickness. In this case will be air thickness between both subtrates
v_d = v_dsas;
 
% % Parameters
v_fmin = 0.4e12;             % minimum value of interesting frequency domain
v_fmax = 2e12;             % maximum value of interesting frequency domain
%%
% Files
v_filenameRef = 'TDS_Air_20210427.dat'; % Insert file name for reference
v_filenameSam = 'TDS_pType-0.3-0.5Ohmcm_20210427.dat'; % Insert file name for sample
v_directory = strcat(pwd,'\Data_files\');
f_FileEref = strcat(v_directory,v_filenameRef);
f_FileEmod = strcat(v_directory,v_filenameSam);

v_Timestep = 0.05E-12; % Timestep of the file

% Reading the files
a_fref = dlmread(f_FileEref,' ',13,0); 
a_fmod = dlmread(f_FileEmod,' ',13,0); 
a_Eref = a_fref(:,2);
a_Emod = a_fmod(:,2);

v_osm = 1/5*sum(a_Emod(1:5));        % Determine offset modulation
a_Emod = a_Emod - v_osm;               % Remove offset

v_osr = 1/5*sum(a_Eref(1:5));        % Determine offset reference
a_Eref = (a_Eref - v_osr);               % Remove offset
%%
%%%%%%%%%%%%%%%%%%%%%%% FFT
a_TimeData=(0:v_Timestep:v_Timestep*(numel(a_Emod)-1))*1e12;

% Cleaning the pulses
a_Eref = f_hamming(v_Timestep,a_Eref);
a_Emod = f_hamming(v_Timestep,a_Emod);

lv_f = lv_f + 1;
figure(lv_f)
plot(a_TimeData,real(a_Eref),a_TimeData,real(a_Emod))
title('Pulses in TD')

[a_RefFreq, a_RefAmpl, a_RefPhasebis] = my_fft(a_TimeData,a_Eref',v_fmin*1e-12,v_fmax*1e-12,pi);
a_RefPhase = a_RefPhasebis;
a_Eref=a_RefAmpl.*exp(1i*a_RefPhase);
[a_ModFreq, a_ModAmpl, a_ModPhasebis] = my_fft(a_TimeData,a_Emod',v_fmin*1e-12,v_fmax*1e-12,pi);
a_ModPhase = a_ModPhasebis;
a_Emod=a_ModAmpl.*exp(1i.*a_ModPhase);

a_freq = a_RefFreq*1e12;
%%
% Plot Frequency Ampl
lv_f = lv_f + 1;
figure(lv_f)
plot(a_freq,a_RefAmpl,a_freq,a_ModAmpl)
title('Pulses in FD')

% Plot phases
lv_f = lv_f + 1;
figure(lv_f)
plot(a_freq,a_RefPhase,a_freq,a_ModPhase,a_freq,-abs(a_ModPhase)+abs(a_RefPhase))
title('Phases')
legend('Ref','Mod','Delta')

% Retrieve initial refrac
a_nreal = 1+(abs(-abs(a_ModPhase)+abs(a_RefPhase)))./(2*pi*a_freq*v_d/c_c);
a_nimag = -(c_c./(v_d*4*pi.*a_freq)).*log(a_ModAmpl./a_RefAmpl);
a_z0 = a_nreal - 1i*a_nimag;

lv_f = lv_f + 1;
figure(lv_f)
plot(a_freq,real(a_z0),'bo',a_freq,imag(a_z0),'ro')
title('Initial guess')
legend('Real','Imag')