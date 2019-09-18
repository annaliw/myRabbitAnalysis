%% prepare k, E, etc. 
% run C_C_DelayPlot.m
% this will also generate perturbative CC models

%% proccess data (first load in data)

nstates = 6; 
phase_data = cat(3, H2_SB12_phase, H2_SB14_phase, H2_SB16_phase); 
mean_data = cat(3, mean(H2_SB12_paramout(:,:,:),3), mean(H2_SB14_paramout(:,:,:),3), mean(H2_SB16_paramout(:,:,:),3)); 

% measured_phases = [unwrap(phase_data(1:nstates,1,1))', ...
%                    unwrap(phase_data(1:nstates,1,2))', ...
%                    unwrap(phase_data(1:nstates,1,3))']; 
% errors_phase =    [phase_data(1:nstates,2,1)', ...
%                    phase_data(1:nstates,2,2)', ...
%                    phase_data(1:nstates,2,3)'];

SB_phase_data = cat(3, ...
               [mean_data(1:nstates,2,1)'; unwrap(phase_data(1:nstates,1,1))'], ...
               [mean_data(1:nstates,2,2)'; unwrap(phase_data(1:nstates,1,2))'], ...
               [mean_data(1:nstates,2,3)'; unwrap(phase_data(1:nstates,1,3))']); 
SB_phase_error = cat(3, ...
               [mean_data(1:nstates,2,1)'; phase_data(1:nstates,2,1)'], ...
               [mean_data(1:nstates,2,2)'; phase_data(1:nstates,2,2)'], ...
               [mean_data(1:nstates,2,3)'; phase_data(1:nstates,2,3)']); 
SB_delay_data = cat(3, ...
               [mean_data(1:nstates,2,1)'; unwrap(phase_data(1:nstates,1,1))'.*(T_L*1000/2/(2*pi))], ...
               [mean_data(1:nstates,2,2)'; unwrap(phase_data(1:nstates,1,2))'.*(T_L*1000/2/(2*pi))], ...
               [mean_data(1:nstates,2,3)'; unwrap(phase_data(1:nstates,1,3))'.*(T_L*1000/2/(2*pi))]); 
SB_delay_error = cat(3, ...
               [mean_data(1:nstates,2,1)'; phase_data(1:nstates,2,1)'.*(T_L*1000/2/(2*pi))], ...
               [mean_data(1:nstates,2,2)'; phase_data(1:nstates,2,2)'.*(T_L*1000/2/(2*pi))], ...
               [mean_data(1:nstates,2,3)'; phase_data(1:nstates,2,3)'.*(T_L*1000/2/(2*pi))]); 

%% make some Coulomb scattering phase arrays           
[sigma01, delay01] = coulombScatteringPhase(0, 1, E); 
[sigma11, delay11] = coulombScatteringPhase(1, 1, E); 
[sigma21, delay21] = coulombScatteringPhase(2, 1, E); 
[sigma02, delay02] = coulombScatteringPhase(0, 2, E); 
[sigma12, delay12] = coulombScatteringPhase(1, 2, E); 
[sigma22, delay22] = coulombScatteringPhase(2, 2, E); 
%% create Ivanov CLC plots
% single Coulomb center, definite l, atomic Hydrogen model

Ivanov_CC0to1 = Ivanov_curve(1, 1, E(2:end))*24.2 + delay01; 
Ivanov_CC1to0 = Ivanov_curve(0, 1, E(2:end))*24.2 + delay11; 
Ivanov_CC1to2 = Ivanov_curve(2, 1, E(2:end))*24.2 + delay11; 

Ivanov_CC0 = Ivanov_curve(0, 1, E)*24.2; 
Ivanov_CC2 = Ivanov_curve(2, 1, E)*24.2; 

% plotfun_compareToTheory(SB_phase_data(:,2:end,:), SB_phase_error(:,2:end,:), E(2:end), ...
%     [Ivanov_CC0to1; Ivanov_CC1to0; Ivanov_CC1to2], ...
%      ["Ivanov CC l=0 to 1", "Ivanov CC l=1 to 0", "Ivanov CC l=1 to 2"], ...
%      200); 


%% create Serov CLC plots
% single Coulomb center, definite l (in tW), Z=2 molecular Hydrogen or Z=1? 

[Serov_Z1, wignerFactor_Z1] = Serov_curve(1, E);
[Serov_Z2, wignerFactor_Z2] = Serov_curve(2, E);

Serov_Z1l0 = Serov_Z1(2:end).*24.2+wignerFactor_Z1(2:end).*delay01; 
Serov_Z1l1 = Serov_Z1(2:end).*24.2+wignerFactor_Z1(2:end).*delay11; 
Serov_Z2l0 = Serov_Z2(2:end).*24.2+wignerFactor_Z2(2:end).*delay01; 
Serov_Z2l1 = Serov_Z2(2:end).*24.2+wignerFactor_Z2(2:end).*delay11; 
  
% plotfun_compareToTheory(SB_phase_data(:,2:end,:), SB_phase_error(:,2:end,:), E(2:end), ...
%     [Serov_Z1l0; Serov_Z1l1; Serov_Z2l0; Serov_Z2l1], ...
%      ["Serov Z=1 l=0", "Serov Z=1 l=1", "Serov Z=2 l=0", "Serov Z=2 l=1"], ...
%      300);   

%% Discretize
% integrate to convert CLC to CC (and discretize wigner delay)
int_lim_E = repmat([11 13 15 17], [5 1])*(1240/wavelength) - repmat(IP(2:end)', [1 4]); 
refsize = size(int_lim_E); 
int_lim_E = sort(reshape(int_lim_E, [1 numel(int_lim_E)])); 
% remove = find(int_lim_E < E(1)); 
% int_lim_E(remove) = []; 
int_lim_ind = int_lim_E; 
for ii=1:numel(int_lim_E) 
    int_lim_ind(ii) = find(abs(E-int_lim_E(ii))<0.1, 1); 
    tol = E(int_lim_ind(ii)) - E(int_lim_ind(ii)-1); 
    int_lim_ind(ii) = find(abs(E-int_lim_E(ii))<tol, 1); 
end
int_lim_ind = reshape(int_lim_ind, [refsize(1) refsize(2)]); 

CCP_s_tmp = (fliplr(unwrap(fliplr(mod(CCP,2*pi))))-2*pi); 
% CCP_d_tmp = (fliplr(unwrap(fliplr(mod(CCP_d,2*pi))))-2*pi); 

Ivanov_CC1to0_disc = zeros([refsize(1), refsize(2)-1]); 
Ivanov_CC1to2_disc = zeros([refsize(1), refsize(2)-1]); 
Ivanov_CC0_disc = zeros([refsize(1), refsize(2)-1]); 
Serov_Z1l1_disc = zeros([refsize(1), refsize(2)-1]);
Serov_CC_disc = zeros([refsize(1), refsize(2)-1]); 
delay11_disc = zeros([refsize(1), refsize(2)-1]); 
CCP_s_disc = zeros([refsize(1), refsize(2)-1]); 
CCP_d_disc = zeros([refsize(1), refsize(2)-1]); 
E_disc = zeros([refsize(1), refsize(2)-1]); 
for ii = 1:refsize(1) % enumerate over vibrational states
    for jj = 2:refsize(2) % enumerate over harmonics
        start_ind = int_lim_ind(ii, jj-1); 
        end_ind = int_lim_ind(ii, jj); 
        Ivanov_CC1to0_disc(ii, jj-1) = trapz(E(start_ind:end_ind), Ivanov_CC1to0(start_ind:end_ind));  
        Ivanov_CC0_disc(ii, jj-1) = trapz(E(start_ind:end_ind), Ivanov_CC0(start_ind:end_ind)); 
        Ivanov_CC1to2_disc(ii, jj-1) = trapz(E(start_ind:end_ind), Ivanov_CC1to2(start_ind:end_ind));  
        Serov_Z1l1_disc(ii, jj-1) = trapz(E(start_ind:end_ind), Serov_Z1l1(start_ind:end_ind));
        Serov_CC_disc(ii, jj-1) = trapz(E(start_ind:end_ind), Serov_Z1(start_ind:end_ind)*24.2);
        CCP_s_disc(ii, jj-1) = trapz(E(start_ind:end_ind), CCP_s_tmp(start_ind:end_ind));
%         CCP_d_disc(ii, jj-1) = trapz(E(start_ind:end_ind), CCP_d_tmp(start_ind:end_ind));
        delay11_disc(ii, jj-1) = (sigma11(end_ind) - sigma11(start_ind)).*(T_L*1000/2/(2*pi)); 
        E_disc(ii, jj-1) = (E(end_ind) + E(start_ind))/2; 
    end
end  

Ivanov_CC1to0_disc = reshape(Ivanov_CC1to0_disc, [1 numel(Ivanov_CC1to0_disc)]); 
Ivanov_CC1to2_disc = reshape(Ivanov_CC1to2_disc, [1 numel(Ivanov_CC1to2_disc)]); 
Ivanov_CC0_disc = reshape(Ivanov_CC0_disc, [1 numel(Ivanov_CC0_disc)]); 
Serov_Z1l1_disc = reshape(Serov_Z1l1_disc, [1 numel(Serov_Z1l1_disc)]); 
Serov_CC_disc = reshape(Serov_CC_disc, [1 numel(Serov_CC_disc)]); 
delay11_disc = reshape(delay11_disc, [1 numel(delay11_disc)]); 
CCP_s_disc = reshape(CCP_s_disc, [1 numel(CCP_s_disc)]); 
% CCP_d_disc = reshape(CCP_d_disc, [1 numel(CCP_d_disc)]); 
E_disc = reshape(E_disc, [1 numel(E_disc)]); 

%% Lifted Ar Wigner delays from Mauritsson
% need to have processed Argon RABBITT data

Ar_Wigner = [...
2.6997245179063363, -45.35276635590114; 
5.785123966942148, -29.13460625059369; 
8.980716253443525, -16.36607906200379]; 

Ar_phase = [Ar_SB12_phase; Ar_SB14_phase; Ar_SB16_phase]; 
Ar_phase(:,1) = unwrap(Ar_phase(:,1)); 
Ar_delay = Ar_phase.*(T_L*1000/2/(2*pi)); 
Ar_slope = [Ar_SB12_slope; Ar_SB14_slope; Ar_SB16_slope]; 
Ar_energy = [12 14 16]*1240/810 - 15.763; 

Ar_XUV = Ar_delay(:,1) - Ar_Wigner(:,2); 

%%
% CCP_s_tmp = CCP_s_disc+2*pi; 
% CCP_d_tmp = CCP_d_disc+2*pi; 
% CCP_s_tmp = (fliplr(unwrap(fliplr(mod(CCP_s,2*pi))))-2*pi); 
% CCP_s_tmp = -CCP_s; 

set(groot,'defaultLineLineWidth',2.0)
figure; hold on; 
errorbar(squeeze(reshape(SB_phase_data(1,2:end,:), [1 15])), ...
         squeeze(reshape(SB_phase_data(2,2:end,:), [1 15])), ...
         squeeze(reshape(SB_phase_error(2, 2:end,:), [1 15])), 'o', ...
         'DisplayName', 'Measurement');
% plot(E, CCP_s_tmp, 'k-', 'DisplayName', 'CC phase l=0'); 
% plot(E_disc, -delay11_disc./(T_L*1000/2/(2*pi)), 'r-', 'DisplayName', 'Coulomb phase l=1'); 
% plot(E_disc, -delay11_disc./(T_L*1000/2/(2*pi)) + CCP_s_tmp(SB_ind), 'b-', 'DisplayName', 'Atomic phase l=0'); 
xlim([1.5 10]); 
legend; 
xlabel('electron kinetic energy (eV)'); 
ylabel('phase (radians)'); 

%%

CCP_s_tmp = (fliplr(unwrap(fliplr(mod(CCP_s,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 

TDSE_H2a = [2.2689, 5.3361, 8.3823; -116.1756, -91.1418, -62.2917]; 
TDSE_diff = [mean(SB_delay_data(2,2:end,1))-TDSE_H2a(2,1), ...
             mean(SB_delay_data(2,2:end,2))-TDSE_H2a(2,2), ...
             mean(SB_delay_data(2,2:end,3))-TDSE_H2a(2,3)];  

set(groot,'defaultLineLineWidth',2.0)
figure; hold on; 
errorbar(squeeze(reshape(SB_delay_data(1,2:end,:), [1 15])), ...
         squeeze(reshape(SB_delay_data(2,2:end,:), [1 15]))-reshape(repmat(TDSE_diff, [5 1]), [1 15]), ...
         squeeze(reshape(SB_delay_error(2, 2:end,:), [1 15])), 'o', ...
         'DisplayName', 'shifted measurement');
plot(TDSE_H2a(1,:), TDSE_H2a(2,:), 'ko-', 'DisplayName', 'TDSE');
plot(E_disc, delay11_disc + CCP_s_tmp(SB_ind), 'b-', 'DisplayName', 'Atomic phase l=0'); 
plot(E_disc, Serov_Z1l1_disc, 'bo-.', 'DisplayName', 'Serov CLC + tW'); 
legend; 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 

%%

CCP_tmp = (fliplr(unwrap(fliplr(mod(CCP,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 
CCPA_tmp = (fliplr(unwrap(fliplr(mod(CCPA,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 
CCPAp_tmp = (fliplr(unwrap(fliplr(mod(CCPAp,2*pi))))-2*pi).*(T_L*1000/2/(2*pi)); 

[Serov_CC_plus, ~] = Serov_curve(1, E + 1240/810);
[Serov_CC_minus, ~] = Serov_curve(1, E - 1240/810);
Serov_CC = (Serov_CC_plus - Serov_CC_minus)*24.2/4/pi; 
clear('Serov_CC_plus', 'Serov_CC_minus'); 
Ivanov_CC_plus = Ivanov_curve(0,1, E + 1240/810);
Ivanov_CC_minus = Ivanov_curve(0,1, E - 1240/810);
Ivanov_CC = (Ivanov_CC_plus - Ivanov_CC_minus)*24.2/4/pi; 
clear('Ivanov_CC_plus', 'Ivanov_CC_minus'); 


TDSE_H2a = [2.2689, 5.3361, 8.3823; -116.5337, -90.1702, -61.2589]; 
TDSE_diff = [mean(SB_delay_data(2,2:end,1))-TDSE_H2a(2,1), ...
             mean(SB_delay_data(2,2:end,2))-TDSE_H2a(2,2), ...
             mean(SB_delay_data(2,2:end,3))-TDSE_H2a(2,3)];  
         
xdata = squeeze(reshape(SB_delay_data(1,2:end,:), [1 15])); 
phase_data = squeeze(reshape(SB_delay_data(2,2:end,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, 2:end,:), [1 15]));

% plot_data = phase_data - ...
%             reshape(repmat(Ar_XUV', [5 1]), [1 15]) - ...
%             reshape(repmat(TDSE_diff, [5 1]), [1 15]); 
plot_data = phase_data - reshape(repmat(Ar_XUV', [5 1]), [1 15]) - ...
                reshape(repmat(TDSE_H2a(2,:), [5 1]), [1 15]) + 0.7741.*(T_L*1000/4/pi); % - Ar_SB18_phase(1).*(T_L*1000/2/(2*pi)); 
plot_error = sqrt(...
                  reshape(repmat(Ar_phase(:,2)'.^2, [5 1]), [1 15]) + ...
                  phase_error.^2); 

% CCP_shift = CCP_s_tmp(SB_ind(7)) - mean(plot_data(6:10)); 
% Ivanov_shift = Ivanov_CC0_disc(7) - mean(plot_data(6:10));
% Serov_shift = Serov_CC_disc(7) - mean(plot_data(6:10));

set(groot,'defaultLineLineWidth',2.0)
figure; hold on; 
errorbar(xdata, -plot_data, plot_error, 'o', 'DisplayName', 'CC measurement');
plot(E, CCP_tmp, 'Color', [0,0,1], 'DisplayName', 'Atomic phase l=0'); 
plot(E, CCPA_tmp, 'Color',[0,0.3,1], 'DisplayName', 'Atomic phase l=0 (A)'); 
plot(E, CCPAp_tmp, 'Color',[0,0.5,1], 'DisplayName', 'Atomic phase l=0 (Ap)'); 
% plot(E_disc, Ivanov_CC0_disc.*(1/2/(2*pi)), 'b--', 'DisplayName', 'Ivanov CLC'); 
% plot(E_disc, Serov_CC_disc.*(1/2/(2*pi)), 'b-.', 'DisplayName', 'Serov CLC'); 
% plot(E, CCP_s_tmp, 'b-', 'DisplayName', 'Atomic phase l=0'); 
plot(E, Ivanov_CC, 'b--', 'DisplayName', 'Ivanov CLC'); 
plot(E, Serov_CC, 'b-.', 'DisplayName', 'Serov CLC'); 
legend; 
xlim([1.6 14]); 
% ylim([-500 500]); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 