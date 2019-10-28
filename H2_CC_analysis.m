%% prepare k, E, etc. 
% run C_C_DelayPlot.m
% this will also generate perturbative CC models

% trim NaN's (if starting at E=0)
E = E(2:end); 
k = k(2:end); 
CCP = CCP(2:end); 
CCPA = CCPA(2:end); 
CCPAp = CCPAp(2:end); 

%% proccess data (first load in data)

nstates = 6; 
phase_data = flipud(cat(3, H2_SB12_phase, H2_SB14_phase, H2_SB16_phase)); 
mean_data = flipud(cat(3, mean(H2_SB12_paramout(:,:,:),3), mean(H2_SB14_paramout(:,:,:),3), mean(H2_SB16_paramout(:,:,:),3))); 

SB_phase_data = cat(3, ...
               [mean_data(1:nstates,2,1)'; unwrap(phase_data(1:nstates,1,1))'], ...
               [mean_data(1:nstates,2,2)'; unwrap(phase_data(1:nstates,1,2))'], ...
               [mean_data(1:nstates,2,3)'; unwrap(phase_data(1:nstates,1,3))' + 2*pi]) - pi; 
SB_phase_error = cat(3, ...
               [mean_data(1:nstates,2,1)'; phase_data(1:nstates,2,1)'], ...
               [mean_data(1:nstates,2,2)'; phase_data(1:nstates,2,2)'], ...
               [mean_data(1:nstates,2,3)'; phase_data(1:nstates,2,3)']); 
SB_delay_data = SB_phase_data; SB_delay_data(2,:) = (SB_delay_data(2,:)+glob_phase).*(T_L*1000/2/(2*pi)); %20191009 data processing accidentally added a 120as phase offset between Ar and H2
SB_delay_error = cat(3, ...
               [mean_data(1:nstates,2,1)'; phase_data(1:nstates,2,1)'.*(T_L*1000/2/(2*pi))], ...
               [mean_data(1:nstates,2,2)'; phase_data(1:nstates,2,2)'.*(T_L*1000/2/(2*pi))], ...
               [mean_data(1:nstates,2,3)'; phase_data(1:nstates,2,3)'.*(T_L*1000/2/(2*pi))]); 
           
Ar_phase = [Ar_SB12_phase; Ar_SB14_phase; Ar_SB16_phase; Ar_SB18_phase]; 
Ar_phase(:,1) = Ar_phase(:,1) - pi; 
Ar_phase(3:4,1) = Ar_phase(3:4,1)+2*pi; 
% Ar_phase(:,1) = unwrap(Ar_phase(:,1)); 
Ar_delay = Ar_phase.*(T_L*1000/2/(2*pi)); 
Ar_slope = [Ar_SB12_slope; Ar_SB14_slope; Ar_SB16_slope; Ar_SB18_slope]; 
Ar_energy = [12 14 16 18]*1240/810 - 15.763; 


%% make some Coulomb scattering phase arrays           
[sigma01, delay01] = coulombScatteringPhase(0, 1, E); 
[sigma11, delay11] = coulombScatteringPhase(1, 1, E); 
[sigma21, delay21] = coulombScatteringPhase(2, 1, E); 
[sigma02, delay02] = coulombScatteringPhase(0, 2, E); 
[sigma12, delay12] = coulombScatteringPhase(1, 2, E); 
[sigma22, delay22] = coulombScatteringPhase(2, 2, E); 

%% Import Anatoli's TDSE values for Argon, H2

% import from CSV and convert to table to match H2 calculations
tmp = csvread('/Users/annaliw/Documents/lab/plots/Calculations/H2_TDSE/ArTDSE_810.csv',0,0); 
x_Ee = tmp(:,1)-15.736; 
% t = 2*flipud(tmp(:,2)); 
t = tmp(:,2); 
n = 12:2:24; 
ArTDSE_810 = table(x_Ee, n', t); 

% import H2 as tables
H2TDSE_810_145 = readtable('/Users/annaliw/Documents/lab/plots/Calculations/H2_TDSE/H2-R145.dat');
H2TDSE_810_150 = readtable('/Users/annaliw/Documents/lab/plots/Calculations/H2_TDSE/H2-R150.dat');
H2TDSE_810_140 = readtable('/Users/annaliw/Documents/lab/plots/Calculations/H2_TDSE/H2-omega056.dat');
H2TDSE_785_140 = readtable('/Users/annaliw/Documents/lab/plots/Calculations/H2_TDSE/H2-omega058.dat');
H2TDSE_760_140 = readtable('/Users/annaliw/Documents/lab/plots/Calculations/H2_TDSE/H2-omega06.dat');

%% get and plot XUV phase

Ar_theory = 2*cat(2, [Ar_energy(1); -48], Ar_Dahlstrom(:,1:3)); 
% Ar_theory = ArTDSE_810.t(1:4); 

figure; errorbar(Ar_energy, Ar_delay(:,1)-pi, Ar_delay(:,2), 'ko', 'DisplayName', 'Argon raw measurement');
hold on; plot(ArTDSE_810.x_Ee(1:4), Ar_theory(2,:), 'rv', 'DisplayName', 'Argon TDSE'); 
xlabel('electron kinetic energy'); 
ylabel('delay (as)'); 
xlim([2 12.5]); 
legend; 
goodplot(); 
hold off; 

XUV_delay = Ar_delay(:,1) - pi - Ar_theory(2,:)'; 
figure; errorbar(Ar_energy, XUV_delay, Ar_delay(:,2), 'ko', 'DisplayName', 'XUV delay'); 
xlabel('electron kinetic energy'); 
ylabel('delay (as)'); 
xlim([2 12.5]); 
legend; 
goodplot(); 
hold off; 

%% *** compare H2TDSE to H2 data ***

xdata = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
phase_data = squeeze(reshape(SB_delay_data(2,1:(end-1),:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, 1:(end-1),:), [1 15]));

% % use the mean of the measured H2 values
% xdata = squeeze(reshape(mean(SB_delay_data(1,2:end,:),2), [1 3])); 
% phase_data = squeeze(reshape(mean(SB_delay_data(2,2:end,:),2), [1 3])); 
% phase_error = squeeze(reshape(mean(SB_delay_error(2,2:end,:),2), [1 3])); 

% subtract out XUV contribution by referencing to Ar data and adding in
% provided Ar TDSE values
plot_data = phase_data - reshape(repmat(XUV_delay(1:3), [1 5])', [1 15]); 
% plot_data = phase_data ... 
%                 - Ar_delay(1:3,1)' ...  
%                 + ArTDSE_810.t(1:3)'/2; 
plot_error = phase_error; 


figure; hold on; 

subplot(2, 3, [1 2 3]); hold on; 
errorbar(xdata, plot_data, plot_error, 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t, H2TDSE_810_140.err, ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t, H2TDSE_810_145.err, ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t, H2TDSE_810_150.err, ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
legend; 
goodplot(); 

subplot(2, 3, 4); hold on; 
errorbar(xdata(1:5), plot_data(1:5), plot_error(1:5), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(1), H2TDSE_810_140.t(1), H2TDSE_810_140.err(1), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(1), H2TDSE_810_145.t(1), H2TDSE_810_145.err(1), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(1), H2TDSE_810_150.t(1), H2TDSE_810_150.err(1), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
goodplot(); 

subplot(2, 3, 5); hold on; 
errorbar(xdata(6:10), plot_data(6:10), plot_error(6:10), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(2), H2TDSE_810_140.t(2), H2TDSE_810_140.err(2), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(2), H2TDSE_810_145.t(2), H2TDSE_810_145.err(2), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(2), H2TDSE_810_150.t(2), H2TDSE_810_150.err(2), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
goodplot(); 

subplot(2, 3, 6); hold on; 
errorbar(xdata(11:15), plot_data(11:15), plot_error(11:15), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(3), H2TDSE_810_140.t(3), H2TDSE_810_140.err(3), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(3), H2TDSE_810_145.t(3), H2TDSE_810_145.err(3), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(3), H2TDSE_810_150.t(3), H2TDSE_810_150.err(3), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
goodplot(); 

%% *** compare ArTDSE to Ar data ***

xdata = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
phase_data = squeeze(reshape(SB_delay_data(2,2:end,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, 2:end,:), [1 15]));

% subtract out XUV contribution by referencing to Ar data and adding in
% provided Ar TDSE values
plot_data = Ar_delay(1:3,1) ...  
                - squeeze(SB_delay_data(2,5,:)) ...
                + H2TDSE_810_150.t(1:3); 

plot_error = Ar_delay(1:3,2); 
% plot_data = tmp + reshape(repmat(ArTDSE_810.t(1:3), [1 5])', [1 15]); 
% xdata = reshape( (repmat(12:2:16, [5 1])'*1240/810 - repmat(IP(2:end), [3 1]))', [1 15]);  

figure; hold on; 
errorbar((12:2:16)*1240/810 - 15.736, plot_data, plot_error, 'ko', 'DisplayName', 'Ar measurement'); 
plot(ArTDSE_810.x_Ee(1:3), ArTDSE_810.t(1:3), ...
    'rv', 'DisplayName', 'Ar TDSE 810nm'); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
legend; 
hold off; 
goodplot()

%% *** Argon CC plot ***
%% calculate CC delays
CCP_tmp = (fliplr(unwrap(fliplr(mod(CCP,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 
CCPA_tmp = (fliplr(unwrap(fliplr(mod(CCPA,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 
CCPAp_tmp = (fliplr(unwrap(fliplr(mod(CCPAp,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 

[Serov_CC_plus, Serov_mult_plus] = Serov_curve(1, E + 1240/810);
Serov_CC_plus = fliplr(cumtrapz(fliplr(E+1240/810), fliplr(Serov_CC_plus))); 
[Serov_CC_minus, Serov_mult_minus] = Serov_curve(1, E - 1240/810); 
Serov_CC_minus = fliplr(cumtrapz(fliplr(E-1240/810), fliplr(Serov_CC_minus))); 
Serov_CC = (Serov_CC_plus - Serov_CC_minus)*24.2/(2*1240/810); 
clear('Serov_CC_plus', 'Serov_CC_minus'); 
Ivanov_CC_plus = fliplr(cumtrapz(fliplr(E+1240/810), fliplr(Ivanov_curve(0,1, E + 1240/810))));
Ivanov_CC_minus = fliplr(cumtrapz(fliplr(E-1240/810), fliplr(Ivanov_curve(0,1, E - 1240/810))));
Ivanov_CC = (Ivanov_CC_plus - Ivanov_CC_minus)*24.2/(2*1240/810); 
clear('Ivanov_CC_plus', 'Ivanov_CC_minus'); 

%% process Argon single photon scattering phases
% actually double photon :(
% % Mauritsson Argon single photon values
% Ar_Wigner = [...
% 2.6997245179063363, -45.35276635590114; 
% 5.785123966942148, -29.13460625059369; 
% 8.980716253443525, -16.36607906200379]; 
% Ar_XUV = Ar_delay(:,1) - Ar_Wigner(:,2); 

% import from CSV and convert to table to match H2 calculations
% energy is in Rydberg, phase in in units of pi. Needs to be converted. 
tmp = csvread('/Users/annaliw/Documents/lab/plots/Calculations/Ar_SinglePhoton/Ar_3StoP.csv',0,0); 
x_Ee = 13.6057*tmp(:,1); 
phi = pi*tmp(:,2); 
phi = interp1(x_Ee, phi, E, 'spline', 'extrap'); 
Ar_3StoP_interp = table(E, phi); 

% find correct energy values for each peak
ArSBind = 11:1:19; 
ArSBind(:) = 0; 
for ii=1:numel(ArSBind)
    n = ii+10; 
    ArSBind(ii) = find(abs(E-n*1240/810+15.736)<0.02, 1); 
end

% get delta phi and convert to tau
% also need to add the p-wave coulomb phase shift (sigma01
t = 12:2:16; 
t(1) = Ar_3StoP_interp.phi(ArSBind(3)) - Ar_3StoP_interp.phi(ArSBind(1)) ...
        + sigma01(ArSBind(3)) - sigma01(ArSBind(1)); 
t(2) = Ar_3StoP_interp.phi(ArSBind(5)) - Ar_3StoP_interp.phi(ArSBind(3)) ...
        + sigma01(ArSBind(5)) - sigma01(ArSBind(3));
t(3) = Ar_3StoP_interp.phi(ArSBind(7)) - Ar_3StoP_interp.phi(ArSBind(5)) ...
        + sigma01(ArSBind(7)) - sigma01(ArSBind(5)); 
t = t' * (T_L*1000/2/(2*pi)); 
x_Ee = ((12:2:16)*1240/810 - 15.736)'; 
Ar_singlephoton = table(x_Ee, t); 

% special Serov version
t = 12:2:16; 
t(1) = (Ar_3StoP_interp.phi(ArSBind(3)) + sigma01(ArSBind(3)))*Serov_mult_plus(ArSBind(3)) ...
    - (Ar_3StoP_interp.phi(ArSBind(1)) + sigma01(ArSBind(1)))*Serov_mult_minus(ArSBind(1)); 
t(2) = (Ar_3StoP_interp.phi(ArSBind(5)) + sigma01(ArSBind(5)))*Serov_mult_plus(ArSBind(5)) ...
    - (Ar_3StoP_interp.phi(ArSBind(3)) + sigma01(ArSBind(3)))*Serov_mult_minus(ArSBind(3)); 
t(3) = (Ar_3StoP_interp.phi(ArSBind(7)) + sigma01(ArSBind(7)))*Serov_mult_plus(ArSBind(7)) ...
    - (Ar_3StoP_interp.phi(ArSBind(5)) + sigma01(ArSBind(5)))*Serov_mult_minus(ArSBind(5)); 
t = t' * (T_L*1000/2/(2*pi)); 
x_Ee = ((12:2:16)*1240/810 - 15.736)'; 
Ar_singlephoton_Serov = table(x_Ee, t); 

%% plot ArCC through subtraction, but this is kind of circular
% xdata = reshape(SB_delay_data(1,2:4,:), [1 9]); 
xdata = (12:2:16)*1240/810 - IP(2); 
phase_data = squeeze(SB_delay_data(2,5,:)); 
phase_error = squeeze(SB_delay_error(2,5,:)); 

% % prepare H2 values
% tmp = [H2TDSE_810_140.t, H2TDSE_810_145.t, H2TDSE_810_150.t];  
% tmp = reshape(tmp(1:3,:)', [1 9]); 
% % tmp = squeeze(mean(tmp(1:3,:),2)'); 

plot_data = phase_data - Ar_delay(:,1) + Ar_singlephoton.t(1:3) - H2TDSE_810_150.t(1:3); % should give Ar CC values
plot_theory = ArTDSE_810.t(1:3)-Ar_singlephoton.t(1:3); 

figure; hold on; 
plot(xdata, -plot_data, 'ko', 'DisplayName', 'CC measurement');
plot(Ar_singlephoton.x_Ee, plot_theory, 'kv', 'DisplayName', 'CC TDSE'); 
plot(E, CCP_tmp, 'Color', [0,0,1], 'DisplayName', 'Dahlstrom CC l=0'); 
plot(E, CCPA_tmp, 'Color',[0,0.3,1], 'DisplayName', 'Dahlstrom CC l=0 (A)'); 
plot(E, CCPAp_tmp, 'Color',[0,0.5,1], 'DisplayName', 'Dahlstrom CC l=0 (Ap)'); 
plot(E, Ivanov_CC, 'b--', 'DisplayName', 'Ivanov CLC'); 
% plot(E, Serov_CC, 'b-.', 'DisplayName', 'Serov CLC'); 
legend; 
xlim([1.6 14]); 
ylim([-800 400]); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
goodplot()

%% compare Argon measurement to separable model

xdata = (12:2:16)*1240/810 - IP(2); 
phase_data = reshape(SB_delay_data(2,5,:), [1 3]); 
phase_error = reshape(SB_delay_error(2,5,:), [1 3]); 

% prepare H2 values
% tmp = [H2TDSE_810_140.t, H2TDSE_810_145.t, H2TDSE_810_150.t];  
% tmp = reshape(tmp(1:3,:)', [1 3]); 

plot_data = phase_data - Ar_delay(:,1)' - H2TDSE_810_150.t(1:3)'; % referenced argon measurement
plot_theory = ArTDSE_810.t(1:3); 

% find Argon energy values in E
E_Ar = 1:3; 
for ii=1:3
    E_Ar(ii) = find(abs(E-Ar_singlephoton.x_Ee(ii)) < 0.02, 1); 
end

figure; hold on; 
plot(xdata, -plot_data, 'o-', 'DisplayName', 'Argon measurement');
plot(Ar_singlephoton.x_Ee, plot_theory, '*-', 'DisplayName', 'Argon TDSE'); 
plot(E(E_Ar), CCP_tmp(E_Ar)+Ar_singlephoton.t(1:3)', '+-', 'Color', [0,0,1], 'DisplayName', 'Dahlstrom CC l=0'); 
plot(E(E_Ar), CCPA_tmp(E_Ar)+Ar_singlephoton.t(1:3)', '+-', 'Color',[0,0.3,1], 'DisplayName', 'Dahlstrom CC l=0 (A)'); 
plot(E(E_Ar), CCPAp_tmp(E_Ar)+Ar_singlephoton.t(1:3)', '+-', 'Color',[0,0.5,1], 'DisplayName', 'Dahlstrom CC l=0 (Ap)'); 
plot(E(E_Ar), Ivanov_CC(E_Ar)+Ar_singlephoton.t(1:3)', 'bs-', 'DisplayName', 'Ivanov CLC'); 
plot(E(E_Ar), Serov_CC(E_Ar)+Ar_singlephoton_Serov.t(1:3)', 'bv-', 'DisplayName', 'Serov CLC'); 
legend; 
% xlim([1.6 14]); 
% ylim([-800 200]); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
goodplot()

%% *** compare measurement + CC model to Wigner delay ***

% CCP_tmp = (fliplr(unwrap(fliplr(mod(CCP,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 
% CCPA_tmp = (fliplr(unwrap(fliplr(mod(CCPA,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 
% CCPAp_tmp = (fliplr(unwrap(fliplr(mod(CCPAp,2*pi))))-2*pi).*(T_L*1000/2/(2*pi));

% TDSE_H2a = [2.2689, 5.3361, 8.3823; -116.1756, -91.1418, -62.2917]; 
% TDSE_diff = [mean(SB_delay_data(2,2:end,1))-TDSE_H2a(2,1), ...
%              mean(SB_delay_data(2,2:end,2))-TDSE_H2a(2,2), ...
%              mean(SB_delay_data(2,2:end,3))-TDSE_H2a(2,3)];  
         
xdata = squeeze(reshape(SB_delay_data(1,2:end,:), [1 15])); 
phase_data = squeeze(reshape(SB_delay_data(2,2:end,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, 2:end,:), [1 15]));

% assume Serov and Ivanov CC delays previously calculated in above cell
tolerance = 0.02; 
H2ind = xdata; 
for ii=1:numel(xdata)
    H2ind(ii) = find(abs(E-xdata(ii))<tolerance, 1); 
end
Arind = Ar_energy; 
for ii=1:numel(Ar_energy)
    Arind(ii) = find(abs(E-Ar_energy(ii))<tolerance, 1); 
end
TDSEind = TDSE_H2a(1,:); 
for ii=1:numel(TDSEind)
    TDSEind(ii) = find(abs(E-TDSE_H2a(1,ii))<tolerance, 1); 
end
Serov_CCsub_data = Serov_CC(H2ind) - reshape(repmat(Serov_CC(Arind), [5 1]), [1 15]);
Serov_CCsub_TDSE = Serov_CC(TDSEind); 
Ivanov_CCsub_data = Ivanov_CC(H2ind)- reshape(repmat(Ivanov_CC(Arind), [5 1]), [1 15]); 
Ivanov_CCsub_TDSE = Ivanov_CC(TDSEind); 

plot_data = phase_data - reshape(repmat(Ar_XUV', [5 1]), [1 15]) ...
                + 0.7741.*(T_L*1000/4/pi); 
plot_error = sqrt(...
                  reshape(repmat(Ar_phase(:,2)'.^2, [5 1]), [1 15]) + ...
                  phase_error.^2); 

set(groot,'defaultLineLineWidth',2.0)
figure; hold on; 
errorbar(xdata, plot_data - Serov_CCsub_data, plot_error, 'bo', ...
         'DisplayName', 'data with Serov CC subtraction');
errorbar(xdata, plot_data - Ivanov_CCsub_data, plot_error, 'co', ...
         'DisplayName', 'data with Ivanov CC subtraction');
plot(TDSE_H2a(1,:), TDSE_H2a(2,:) - Serov_CCsub_TDSE, 'bv-', ...
    'DisplayName', 'TDSE with Serov CC subtraction');
plot(TDSE_H2a(1,:), TDSE_H2a(2,:) - Ivanov_CCsub_TDSE, 'cv-', ...
    'DisplayName', 'TDSE with Ivanov CC subtraction');
plot(E(2:end), delay11, 'k-', 'DisplayName', 'Wigner Delay l=1'); 
legend;  
xlim([1.5 10]); 
ylim([-100 400]); 
title('H2 photoionization delay calculated with different CC models')
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 

%% *** compare measurement and TDSE ***
% can only compare relative phases within vib. peak because of XUV phase
% uncomment TDSE parts when Vlad is done calculating them

TDSE_H2a = [2.2689, 5.3361, 8.3823; -116.1756, -91.1418, -62.2917]; 
% TDSE_H2a_vibrelphase = TDSE_H2a; 
% for ii=1:size(TDSE_H2a_vibrelphase,3)
%     TDSE_H2a_vibrelphase(2,2:end,ii) = TDSE_H2a(2,2:end,ii) - TDSE_H2a(2,end,ii); 
% end
SB_vibrelphase = SB_phase_data; 
for ii=1:size(SB_delay_data,3)
    SB_vibrelphase(2,2:end,ii) = SB_phase_data(2,2:end,ii) - SB_phase_data(2,end,ii); 
end

xdata = squeeze(reshape(SB_vibrelphase(1,2:end,:), [1 15])); 
phase_data = squeeze(reshape(SB_vibrelphase(2,2:end,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, 2:end,:), [1 15]));
% xcalc = squeeze(reshape(TDSE_H2_vibrelphase(1,2:end,:), [1 numel(TDSE_H2_vibrephase(1,2:end,:))])); 
% phase_calc = squeeze(reshape(TDSE_H2_vibrelphase(2,2:end,:), [1 numel(TDSE_H2_vibrelphase(2,2:end,:))])); 

figure; hold on; 
subplot(2,1,1); 
% errorbar(xdata, phase_data, phase_error, 'o', 'DisplayName', 'relative phase w/in v-state'); 
plot(xdata, phase_data, 'o', 'DisplayName', 'measured'); 
ylabel('phase (radians)'); 
title('relative phase w/in v-state'); 
ylim([-0.5 0.5]); 
legend; 
% subplot(2,1,2); 
% plot(xcalc, phase_calc, 'o', 'DisplayName', 'relative phase w/in v-state'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('phase (radians)'); 
% ylim([-0.5 0.5]); 
% legend; 

%% *** RABBITT spectrogram ***
h1 = plotfun_rabbitspectrum(9:2:19, 810, Ebins, sum(sum(E_SpectraArray,2),3), 'average'); 
ax1 = gca; 
h2 = plotfun_rabbitspectrum(10:2:18, 810, Ebins, twoOmega_signal, 'twoOmega'); 
ax2 = gca; 
%%
h3 = figure; hold on; 
tmp = fft(E_SpectraArray,[],2) .* exp(1j*permute(repmat(peak_phase, [900 1 223]), [1 3 2])); 
tmp = ifft(fftshift(tmp,2),[],2);
carpet = abs(sum(tmp, 3));% - median(median(abs(sum(tmp,3)),1),2); 
% carpet = carpet - repmat(mean(carpet,1),[900 1]);
% carpet = carpet - repmat(XUV_only./sum(XUV_only), [1 223]); 
% carpet = carpet - repmat(mean(carpet,2),[1 223]); 
map = zeros(101 , 3);
map(1:50,1) = 1; map(1:50,2) = linspace(0,1,50)'; map(1:50,3) = linspace(0,1,50)'; 
map(51:101,1) = linspace(1,0,51)'; map(51:101,2) = linspace(1,0,51)'; map(51:101,3) = 1; 
imagesc(Ebins, stageTimes, carpet'); 
% colormap(map); 
colorbar; 
% xlim([1.5 15]); 
% ylim([-1 5]*10^(-15))
ax3 = gca; 

%%
tmp = fft(E_SpectraArray,[],2) .* exp(1j*permute(repmat(peak_phase, [900 1 223]), [1 3 2])); 
tmp = ifft(fftshift(tmp,2),[],2); 