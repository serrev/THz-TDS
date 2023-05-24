%% Version 1.5 
% Program to perform TDS Analysis
% 
% Along time ago, in a code far, far away...

% In this program we are using the scheeme of having a sample sandwiched in two substrates.
% (Air|Substrate|Sample|Substrate|Air). It can be easily extend to a different scheeme 

% In this code: Sample means: Substrate+Sample+Substrate
% In this code: Reference means: Substrate+Air+Substrate

% Developed by: Sergio Revuelta
%%
% Cleaning everything before starts
clc, clf, close all, clear all
%%
f_namedata = 'Sample';
% Constants
c_c = 299792458; % Light velocity
c_e0 = 8.8541878128e-12;
% Variables
lv_f = 0;
v_dsu = 1.06E-3; % Substrate thickness
v_dsas = 0.32e-3; % Sample thickness. In this case will be the sample thickness between both subtrates
v_dsar = 0.25e-3; % Sample thickness. In this case will be air thickness between both subtrates
v_d = v_dsas;
 
% % Parameters
v_fmin = 0.3e12;             % minimum value of interesting frequency domain
v_fmax = 1.8e12;             % maximum value of interesting frequency domain
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
a_Emod = (a_Emod - v_osm);               % Remove offset

v_osr = 1/5*sum(a_Eref(1:5));        % Determine offset reference
a_Eref = (a_Eref - v_osr);               % Remove offset
%%
%%%%%%%%%%%%%%%%%%%%%%% FFT
a_TimeData=(0:v_Timestep:v_Timestep*(numel(a_Emod)-1))*1e12;
a_Ereftime = a_Eref;
a_Emodtime = a_Emod;
% Cleaning the pulses
% a_Eref = f_hamming(v_Timestep,a_Eref);
% a_Emod = f_hamming(v_Timestep,a_Emod);
%
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
plot(a_freq,a_RefPhase,a_freq,a_ModPhase,a_freq,(a_ModPhase)-(a_RefPhase))
title('Phases')
legend('Ref','Mod','Delta')

% Retrieve initial refrac
a_nreal = 1+(abs(abs(a_ModPhase)-abs(a_RefPhase)))./(2*pi*a_freq*v_d/c_c);
a_R = ((1-a_nreal)/(1+a_nreal)).^2;
a_nimag = -(c_c./(v_d*4*pi.*a_freq)).*log((a_ModAmpl*(1-a_R))./a_RefAmpl);	
a_z0 = a_nreal + 1i*a_nimag;

lv_f = lv_f + 1;
figure(lv_f)
plot(a_freq,real(a_z0),'bo',a_freq,abs(imag(a_z0)),'ro')
title('Initial guess')
legend('Real','Imag')

%% Getting matrix of solutions
% Refract matrix (triangular matrix, because yes... In principle the resolution is higher with less datapòints that one rectangular matrix)
% z = x + 1i*y; x:[xb,xe]; y:[yb,ye]; r: step
v_xb = 1;
v_xe = 12;
v_yb = 0.1;
v_ye = 3;

v_r = 0.005;

v_X=v_xe-v_xb;
v_Y=v_ye-v_yb;

v_n=ceil(v_Y/v_r+1); %round to the nearest int
v_dy=v_Y/(v_n-1); 
v_m=ceil(v_X/sqrt(v_r^2-v_dy^2/4) +1); %Round to nearest int following a triangule
v_dx = v_X/(v_m-1);

a_vx=linspace(v_xb,v_xe,v_m); % initial xvector 
a_vy=linspace(v_yb,v_ye,v_n); % initial yvector
[a_x,a_y]=meshgrid(a_vx,a_vy); % Initial mesh

% m_refrac=[x,y]; %Final mesh of elements matrix array of 2 columns (Re|Im)
% m_refrac = [vx,vy];
m_refrac = a_x + 1i.*a_y;
m_diff = zeros(size(m_refrac)); % Solutions per frequency
% m_diff2 = zeros(size(m_refrac)); % Aux for 2 max value
% m_diff3 = zeros(size(m_refrac)); % Aux for 3 max value
% m_diff4 = zeros(size(m_refrac)); % Aur for 4 max value
a_sol = zeros(size(a_freq));

[v_numf,v_numc] = size(m_refrac);
%% Initialization

%%%%%%%%%%%%%%%%%%%%%%%

%	AIR |	AIR	|	MEDIUM3	|	AIR	|	AIR	
% If we have a sample in a substrate sandwich (AIR|SUBS|SAM|SUBS|AIR)
% If we have a sample on top of a substrate (AIR|SUBS|SAM|AIR|AIR)% In any case, first extract the refractive index for substrates (without sample, i.e. AIR|AIR|AIR|AIR|AIR) then, with the data obtained, extract the refractive index of the sample

% lv_n0 = 1.0005;
% lv_n1 = 1.0005;
% % lv_n1 = m_refrac;
% lv_n2 = m_refrac  ;
% % lv_n2 = 1.0005;
% lv_n3 = 1.0005;
% % lv_n3 = m_refrac;
% lv_n4 = 1.0005;

a_nones = ones(size(a_freq));

% Loading refractive index of substrate obtained in program S_TDS_ASASA.m
load('Substrate.mat')
a_refrac = (conj(a_refracinit));

% Refractive index
a_n0 = 1.0005.*a_nones;
% a_n1 = m_refrac;
a_n1 = a_refrac.*a_nones;
% a_n2 = 1.0005.*a_nones;+
a_n2s = m_refrac;
a_n2r = 1.0005.*a_nones;
% a_n3 = m_refrac;
a_n3 = a_refrac.*a_nones;
a_n4 = 1.0005.*a_nones;

v_d0 = 1;
v_d1 = v_dsu;
v_d2s = v_dsas;
v_d2r = v_dsar;
v_d3 = v_dsu;
v_d4 = 1;

v_numberroots = 100;

%% Function to solve
for lv_k = 1:length(a_freq)
	clear m_diff v_val1 v_val2 v_indf v_indc m_diff2 m_diff3 
    clear v_indf2 v_indf3 v_indc2 v_indc3 v_sol1 v_sol2 v_sol3 v_sol4
    clear v_checksol1 v_checksol2 v_checksol3 v_checksol4 
    clear v_indf4 v_indc4 a_checksol v_getmaxsol v_indgetmaxsol
    clear lv_freq
% 	for lv_i = 1:length(vx)
% 		for lv_j = 1:length(vy)
% 			clear m_fun lv_n2
% 			lv_n2 = vx(lv_i) + 1i*vy(lv_j);
% 			m_fun = f_trans(lv_n0,lv_n1)*f_propagation(lv_n1,a_freq(lv_k),c_c,v_d1)*f_trans(lv_n1,lv_n2)*f_propagation(lv_n2,a_freq(lv_k),c_c,v_d2)*f_FP(lv_n0,lv_n1,lv_n2,a_freq(lv_k),c_c,v_d1)*...
% 				f_trans(lv_n2,lv_n3)*f_propagation(lv_n3,a_freq(lv_k),c_c,v_d3)*f_FP(lv_n1,lv_n2,lv_n3,a_freq(lv_k),c_c,v_d2)*f_trans(lv_n3,lv_n4)*f_propagation(lv_n4,a_freq(lv_k),c_c,v_d4)*...
% 				f_FP(lv_n2,lv_n3,lv_n4,a_freq(lv_k),c_c,v_d3);
% 			m_diff(lv_i,lv_j) = abs((m_fun - a_Eref(lv_k)/a_Emod(lv_k)));
% 		end
% 	end

    clear m_fun lv_n0 lv_n1 lv_n2 lv_n3 lv_n4
    
    lv_n0 = a_n0(lv_k);
    lv_n1 = a_n1(lv_k);
    lv_n2s = a_n2s;
	lv_n2r = a_n2r(lv_k);
    lv_n3 = a_n3(lv_k);
    lv_n4 = a_n4(lv_k);
    
    % Obtaining reference and sample transfer functions

    m_refe = f_getcomplexfield(lv_n0,lv_n1,lv_n2r,lv_n3,lv_n4,c_c,v_d0,v_d1,v_d2r,v_d3,v_d4,a_freq(lv_k));
    m_meas = f_getcomplexfield(lv_n0,lv_n1,lv_n2s,lv_n3,lv_n4,c_c,v_d0,v_d1,v_d2s,v_d3,v_d4,a_freq(lv_k)); 

%     m_refe = f_trans(lv_n0,lv_n1).*f_propagation(lv_n1,a_freq(lv_k),c_c,v_d1).*f_trans(lv_n1,lv_n2r).*f_propagation(lv_n2r,a_freq(lv_k),c_c,v_d2r).*...
%         f_trans(lv_n2r,lv_n3).*f_propagation(lv_n3,a_freq(lv_k),c_c,v_d3).*f_trans(lv_n3,lv_n4);
    
%     m_refe = f_trans(lv_n1,lv_n2).*f_propagation(lv_n1,a_freq(lv_k),c_c,v_d2s).*f_trans(lv_n1,lv_n3);
    
    % Obtaining reference transfer function (pulse obtained with sample)
%     m_meas = f_trans(lv_n0,lv_n1).*f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,a_freq(lv_k),c_c,v_d2s).*...
%         f_trans(lv_n2,lv_n3).*f_trans(lv_n3,lv_n4).*f_FP(lv_n1,lv_n2,lv_n3,a_freq(lv_k),c_c,v_d2s);
    
%     m_meas = f_trans(lv_n0,lv_n1).*f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,a_freq(lv_k),c_c,v_d2s).*...
%         f_trans(lv_n2,lv_n3).*f_trans(lv_n3,lv_n4);
    
%     m_meas = f_trans(lv_n0,lv_n1).*f_trans(lv_n1,lv_n2s).*f_propagation(lv_n2s,a_freq(lv_k),c_c,v_d2s).*...
%         f_trans(lv_n2s,lv_n3).*f_propagation(lv_n3,a_freq(lv_k),c_c,v_d3).*f_trans(lv_n3,lv_n4);
    
    
    % Auxiliar matrix
    m_fun = (conj(m_meas./m_refe) - ((a_Emod(lv_k)/a_Eref(lv_k))));

    m_diff = abs(real(m_fun)) + abs(imag(m_fun));

    %% Solving
%     N = 10;
    clear a_inds a_sols a_checksols a_checksolsb
    clear a_tmpm_diff a_rows a_cols a_indsb
    a_inds = zeros(v_numberroots,1);
    a_sols = zeros(v_numberroots,1);
    a_checksols = zeros(v_numberroots,1);
    a_checksolsb = zeros(v_numberroots,2);
    
    a_tmpm_diff = m_diff(:);
    for lv_i=1:v_numberroots
      [~, a_inds(lv_i)] = min(a_tmpm_diff);
      a_tmpm_diff(a_inds(lv_i)) = inf;
    end
    
    [a_rows, a_cols] = ind2sub(size(m_diff), a_inds);
    a_indsb = nonzeros(a_inds');
    a_inds = a_indsb';
    
    if lv_k == 1
        for lv_i = 1:length(a_inds)
            a_sols(lv_i) = a_vx(a_cols(lv_i)) + 1i*a_vy(a_rows(lv_i));
            a_checksols(lv_i) = abs(abs(a_sols(lv_i))/abs(a_z0(lv_k)));
%             a_checksols(lv_i) = real(a_sols(lv_i)/real(3.4)); %7.87 is the
%             initial result @ 1 THz
            a_checksolsb(lv_i,:) = [a_sols(lv_i) a_checksols(lv_i)];
        end
    else
        for lv_i = 1:length(a_inds)
            a_sols(lv_i) = a_vx(a_cols(lv_i)) + 1i*a_vy(a_rows(lv_i));
%             a_checksols(lv_i) = abs(abs(a_sols(lv_i))/(real(a_sol(lv_k-1)))); %Here I put this to choose the value closer to the previous obtained
            a_checksols(lv_i) = abs(abs(a_sols(lv_i))/abs(a_z0(lv_k)));
            a_checksolsb(lv_i,:) = [a_sols(lv_i) a_checksols(lv_i)];
        end
    end
%     a_checksol = [v_sol1,v_sol2,v_sol3,v_sol4;v_checksol1,v_checksol2,...
%                   v_checksol3,v_checksol4];
    clear v_getmaxsol v_indgetmaxsol a_solbis
    [v_getmaxsol,v_indgetmaxsol] = min(abs(log10(a_checksolsb(:,2))));
                  
	a_sol(lv_k) = conj(a_checksolsb(v_indgetmaxsol,1));
    %% Testing purposes
%     lv_f = lv_f + 1;
%     figure(lv_f)
%     contourf(x,y,log10(m_diff))
%     colorbar('southoutside');
%     lv_freq = a_freq(lv_k)/1E12;
% %     if lv_freq > 1
%         lv_freq = a_freq(lv_k)/1E12;
%         lv_f = lv_f + 1;
%         figure(lv_f)
%         surfc(a_vx,a_vy,log10(m_diff),'EdgeColor','none')
%         view(2)
%         colorbar('southoutside');
%         title(sprintf('Frequency %dTHz',lv_freq))
% %     end
    clc
    disp('-------------------------')
    disp('-------------------------')
    str_log1 = sprintf('%0.2f %% completed',round(100*lv_k/length(a_freq),2));
    str_log2 = sprintf('Frequency %0.2f THz processed. Limit in %0.2f THz',a_freq(lv_k)/1E12,a_freq(length(a_freq))/1E12);
    disp(str_log1);
    disp(str_log2);

end

%% Plotting results
% lv_f = lv_f + 1;
% figure(lv_f)
% plot(a_freq,real(a_sol),'bo',a_freq,imag(a_sol),'ro')

% a_ratr = zeros(size(a_freq));
% a_rati = zeros(size(a_freq));
% a_ratr = real(a_sol)./real(a_z0);
% a_rati = imag(a_sol)./imag(a_z0);

a_dielre = (real(a_sol).^2-imag(a_sol).^2);
a_dielim = 2*real(a_sol).*imag(a_sol);
a_diel = a_dielre + 1i.*a_dielim;
% a_condre = real(a_sol).*imag(a_sol).*a_freq;
% a_condim = (1 - (real(a_sol).^2-imag(a_sol).^2)).*(a_freq./2);
% a_cond = a_condre + 1i.*a_condim;
% a_refracback = 0;
% a_dielback = a_refracback^2;
% a_cond = -1i.*((a_diel-a_dielback).*8.85E-12).*2*pi.*a_freq;

v_permittinf = 12; %Silicon
a_recond = 2*real(a_sol).*imag(a_sol).*2*pi.*a_freq*c_e0;
a_imcond = (v_permittinf-real(a_sol).^2+imag(a_sol).^2).*2*pi.*a_freq*c_e0;
a_cond = a_recond + 1i.*a_imcond;

lv_f = lv_f + 1;
figure(lv_f)
subplot(2,1,1)
plot(a_freq,real(a_z0),'b-',a_freq,real(a_sol),'bo')
legend('Initial','Final')
subplot(2,1,2)
plot(a_freq,imag(a_z0),'r-',a_freq,imag(a_sol),'ro')
legend('Initial','Final')

% a_aux1 = ones(size(a_freq));

% lv_f = lv_f + 1;
% figure(lv_f)
% subplot(2,1,1)
% plot(a_freq,a_ratr,'bo',a_freq,a_aux1,'b-')
% subplot(2,1,2)
% plot(a_freq,a_rati,'ro',a_freq,a_aux1,'r-')

lv_f = lv_f + 1;
figure(lv_f)
subplot(2,1,1)
plot(a_freq,real(a_diel),'bo')
legend('Real permittivity')
subplot(2,1,2)
plot(a_freq,imag(a_diel),'ro')
legend('Imag permittivity')

lv_f = lv_f + 1;
figure(lv_f)
subplot(2,1,1)
plot(a_freq,real(a_cond),'bo')
legend('Real conductivity')
subplot(2,1,2)
plot(a_freq,imag(a_cond),'ro')
legend('Imag conductivity')


%% Saving results
% Save file
f_namedatamat = strcat(f_namedata,'.mat');
delete(f_namedatamat)
a_refracinit = a_sol;
a_freqinit = a_freq;
save(f_namedatamat,'a_freqinit','a_refracinit');


m_rreal = real(a_sol);
m_rimag = imag(a_sol);
m_ereal = real(a_diel);
m_eimag = imag(a_diel);
m_creal = real(a_cond);
% m_cimag = imag(a_cond);
m_amplmod = a_ModAmpl;
m_amplref = a_RefAmpl;
m_phasemod = a_ModPhase;
m_phaseref = a_RefPhase;
% m_time = a_TimeData;
% m_sigref = a_Ereftime;
% m_sigmod = a_Emodtime;

m_saver = [a_freq; m_rreal; m_rimag; m_ereal; m_eimag; m_creal; m_amplmod; m_amplref; m_phasemod; m_phaseref];
newfilename = sprintf('%s%s',v_directory,f_namedata,'_TDSresultsFD', '.dat');
fid = fopen(newfilename, 'wt');
%Header
fprintf(fid,'Frequency(Hz) Real_refrac Imag_refrac Real_diel Imag_diel Real_cond AmplitudeSample AmplitudeRef PhaseSample PhaseRef\n');
%Body
fprintf(fid,'%6.2f %12.8f %12.8f  %12.8f  %12.8f %12.8f %12.8f %12.8f %12.8f  %12.8f\r\n', m_saver);
fclose(fid);

% m_saverfft = [m_time; m_sigref; m_sigmod];
% newfilename = sprintf('%s%s',v_directory,f_namedata,'_TDSresultsTD', '.dat');
% fid = fopen(newfilename, 'wt');
% %Header
% fprintf(fid,'Time SignalSample SignalRef\n');
% %Body
% fprintf(fid,'%6.2f %12.8f  %12.8f\r\n', m_saver);
% fclose(fid);
