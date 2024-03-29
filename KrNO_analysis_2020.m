% import data from most recent KrNO .mat file
% read in theory Kr PI values
% T = readtable('/Users/annawang/Documents/data/Kr_theory_Kheifets/Kr_4p_phase.csv'); 
T = readtable('/Users/annawang/Documents/data/Kr_theory_Kheifets/Kr_theory_HJW_kheifets.csv'); 
Kr_theory_photonE = T.Var1; % note these are the 14eV state values (P 3/2)
% Kr_theory_phase = T.Var2*pi; 
Kr_theory_phase = T.Var2 / (T_L*1000/2/(2*pi)); % convert to phase from as

% choose correct values
sidebands = [8 10 12 14 16 18 20]; 
Kr_T = zeros([2, numel(sidebands)]); 
for ii=1:1:numel(sidebands)
    [~,Kr_4ptheory_Eind] = min(abs(Kr_theory_photonE-sidebands(ii)*1240/wavelength)); % here energy is in photon energy, not photoelectron energy
    Kr_T(:,ii) = [Kr_theory_photonE(Kr_4ptheory_Eind); Kr_theory_phase(Kr_4ptheory_Eind)]; 
end

%% get CC delays
% run C_C_DelayPlot.m
CCmodel = CCPAp; 

Kr_CC = zeros([2, numel(sidebands)]); 
NO_CC = zeros([2 numel(sidebands)]); 
for ii=1:1:numel(sidebands)
    [~,Kr_14_Eind] = min(abs(E-(sidebands(ii)*1240/wavelength-IP(1)))); % only using Ip=14eV
    Kr_CC(:,ii) = [E(Kr_14_Eind); CCmodel(Kr_14_Eind)]; 
    
    [~,NO_X_Eind] = min(abs(E-(sidebands(ii)*1240/wavelength-9.55))); 
    NO_CC(:,ii) = [E(NO_X_Eind); CCmodel(NO_X_Eind)]; 
end

Kr_CC(2,:) = unwrap(mod(Kr_CC(2,:),2*pi)); 
NO_CC(2,:) = unwrap(NO_CC(2,:)); 

%% get XUV phase from Kr 12-20
Kr_phase = [Kr_SB12_phase(:,1) Kr_SB14_phase(:,1) Kr_SB16_phase(:,1) Kr_SB18_phase(:,1) Kr_SB20_phase(:,1)]; 
Kr_phase(2,:) = unwrap(Kr_phase(2,:)); 
Kr_phase(2,5) = Kr_phase(2,5)+pi;
XUV = unwrap(mod(Kr_phase(2,:) - Kr_T(2,3:7) - Kr_CC(2,3:7), 2*pi)); 
% XUV = XUV./2; 
% XUV = mod(Kr_phase(2,:) - Kr_CC(2,3:6), 2*pi); 

% fit line to SB12-16 XUV data
p = polyfit(12:2:16, XUV(1:3), 1); 
XUV_fit = polyval(p, 8:2:20); 
% XUV_fit(end) = XUV(end); 

figure; hold on; 
plot(12:2:20, XUV, 'ko', 'DisplayName', 'XUV data'); 
plot([8 10 20], [XUV_fit(1:2) XUV_fit(7)], 'bo', 'DisplayName', 'XUV extrapolation'); 
plot(7:0.1:21, polyval(p, 7:0.1:21), 'r-', 'DisplayName', 'fit'); 
legend; 
xlabel('sideband'); ylabel('phase'); 
goodplot(24); 

%% read in Taran's experiment
NO_LCLS = readtable('/Users/annawang/Box/data/PIdelays_Anna20210406.txt', 'FileType', 'text', 'Delimiter', ','); 
% photonEnergy_eV_, delay_as_, stdErr_delay__as_
% change names to more useful names
NO_LCLS.Properties.VariableNames{'photonEnergy_eV_'} = 'photonE_eV';
NO_LCLS.Properties.VariableNames{'delay_as_'} = 'delay_as';
NO_LCLS.Properties.VariableNames{'stdErr_delay__as_'} = 'delay_err_as';
%% put it all together!

NO_exp = [NO_SB8_phase(1,1) NO_SB10_phase(1,1) NO_SB12_phase(1,1) NO_SB14_phase(1,1) NO_SB16_phase(1,1) NO_SB18_phase(1,1) NO_SB20_phase(1,1)+pi; ...
          NO_SB8_phase(1,2) NO_SB10_phase(1,2) NO_SB12_phase(1,2) NO_SB14_phase(1,2) NO_SB16_phase(1,2) NO_SB18_phase(1,2) NO_SB20_phase(1,2)]; 
NO_exp = NO_exp'; 
NO_exp(:,1) = unwrap(NO_exp(:,1)); % + pi; 
% NO_exp(1,1) = NO_exp(1,1)+2*pi; 
% NO_exp = NO_RABBITT_phase; 

% Kr_exp = Kr_phase(2,:); % only take Ip=14eV

% NO_PI = (NO_exp(:,1)'-Kr_exp) - NO_CC(2,:) + Kr_CC(2,:) + Kr_T(2,:); 
% NO_PI = (NO_exp(:,1)'-Kr_exp)
NO_streak = NO_exp(:,1)' - XUV_fit; % + 800/(T_L*1000/4/pi); 
% NO_streak = NO_exp(:,1)' - XUV_fit(1:6); 
% NO_streak(3:6) = NO_exp(3:6,1)' - XUV; 

% figure; hold on; 
% errorbar(sidebands*1240/810-9.55, NO_streak, NO_exp(:,2)', 'kd', 'DisplayName', 'NO "streak" phase'); 
% xlabel('photoelectron energy (eV)'); ylabel('phase (radians)'); 
% legend; goodplot(24); 
% 
% figure; hold on; 
% plot(Kr_theory_photonE, Kr_theory_phase, 'rd-', 'DisplayName', 'Kr PI theory'); 
% xlabel('photon energy (eV)'); ylabel('phase (radians)'); 
% legend; goodplot(24); 
% 
% figure; hold on; 
% plot(Kr_CC(1,:), Kr_CC(2,:), 'ro-', 'DisplayName', 'Kr CC'); 
% plot(NO_CC(1,:), NO_CC(2,:), 'ko-', 'DisplayName', 'NO CC'); 
% xlabel('photoelectron energy (eV)'); ylabel('phase (radians)'); 
% legend; goodplot(24); 
% 
% figure; hold on; 
% errorbar(sidebands*1240/810-9.55, NO_exp(:,1)', NO_exp(:,2)', 'ko', 'DisplayName', 'NO RABBIT');
% xlabel('photoelectron energy (eV)'); ylabel('phase (radians)'); 
% legend; goodplot(24); 

% NO_PI = NO_exp(:,1)' - XUV_fit(1:6) - NO_CC(2,:); 
NO_PI = NO_streak - NO_CC(2,:); 

% figure; hold on; 
% errorbar(sidebands*1240/810-9.55, NO_PI*(T_L*1000/2/(2*pi)), NO_exp(:,2)'*(T_L*1000/2/(2*pi)), 'kd', 'DisplayName', 'NO PI'); 
% plot(E(2:end), delay01, 'b--', 'DisplayName', 'Coulomb Delay l=0, Z=1'); 
% plot(E(2:end), delay02, 'r--', 'DisplayName', 'Coulomb Delay l=0, Z=2'); 
% % plot(E(2:end), delay03, 'g--', 'DisplayName', 'Coulomb Delay l=0, Z=3'); 
% xlabel('photoelectron energy (eV)'); ylabel('delay (as)'); 
% legend; goodplot(24); 
% % ylim([-300 2500])
% xlim([1 20])

figure; hold on; 
errorbar(sidebands*1240/810-9.55, NO_streak*(T_L*1000/2/(2*pi)), NO_exp(:,2)'*(T_L*1000/2/(2*pi)), 'kd', 'DisplayName', 'valence ionization delay'); 
errorbar(NO_LCLS.photonE_eV-543.65, NO_LCLS.delay_as, NO_LCLS.delay_err_as, 'rd', 'DisplayName', 'core ionization delay'); 
% plot(E(2:end), delay01 + CCmodel(2:end)*(T_L*1000/2/(2*pi)), 'b--', 'DisplayName', 'Coulomb Delay Z=1 + CCPAp'); 
% plot(E(2:end), delay02 + CCmodel(2:end)*(T_L*1000/2/(2*pi)), 'r--', 'DisplayName', 'Coulomb Delay Z=2 + CCPAp'); 
% plot(E(2:end), delay03 + CCmodel(2:end)*(T_L*1000/2/(2*pi)), 'g--', 'DisplayName', 'Coulomb Delay Z=3 + CCPAp'); 
plot(E, delay01 + Serov_CC_Z1, 'b-.', 'DisplayName', 'Coulomb field'); 
% plot(E(2:end), delay01 + Ivanov_CC_Z1(2:end), 'b:', 'DisplayName', 'Coulomb Delay Z=1 + Ivanov Z=1'); 
% plot(E(2:end), delay02 + Serov_CC_Z2(2:end), 'r-.', 'DisplayName', 'Coulomb Delay Z=2 + Serov Z=2'); 
% plot(E(2:end), delay02 + Ivanov_CC_Z2(2:end), 'r:', 'DisplayName', 'Coulomb Delay Z=2 + Ivanov Z=2'); 
% plot(E(2:end), delay03 + Serov_CC_Z3(2:end), 'g-.', 'DisplayName', 'Coulomb Delay Z=3 + Serov Z=3'); 
% plot(E(2:end), delay03 + Ivanov_CC_Z3(2:end), 'g:', 'DisplayName', 'Coulomb Delay Z=3 + Ivanov Z=3'); 

xlabel('photoelectron energy (eV)'); ylabel('delay (as)'); 
legend; goodplot(24); 
ylim([-300 2500])
xlim([1 20])


%% add in 2.4um CLC and compare to 800 CLC
% run C_C_DelayPlot.m with 2.4um
CCmodel = CCPAp; 

Kr_CC_2400 = zeros([2, numel(sidebands)]); 
NO_CC_2400 = zeros([2 numel(sidebands)]); 
% calculate CC delay at experiment energy points
for ii=1:1:numel(sidebands)
    [~,Kr_14_Eind] = min(abs(E-(sidebands(ii)*1240/wavelength-IP(1)))); % only using Ip=14eV
    Kr_CC_2400(:,ii) = [E(Kr_14_Eind); CCmodel(Kr_14_Eind)]; 
    
    [~,NO_X_Eind] = min(abs(E-(sidebands(ii)*1240/wavelength-9.55))); 
    NO_CC_2400(:,ii) = [E(NO_X_Eind); CCmodel(NO_X_Eind)]; 
end

% Kr_CC(2,:) = unwrap(Kr_CC(2,:)); 
NO_CC_2400(2,:) = unwrap(NO_CC_2400(2,:)); 

% figure; hold on; 
% errorbar(sidebands*1240/810-9.55, (NO_PI + NO_CC_2400(2,:))*(T_L*1000/2/(2*pi)), NO_exp(:,2)'*(T_L*1000/2/(2*pi)), 'kd', 'DisplayName', 'NO PI + 2.4um CLC'); 
% errorbar(sidebands*1240/810-9.55, NO_streak*(T_L*1000/2/(2*pi)), NO_exp(:,2)'*(T_L*1000/2/(2*pi)), 'ko', 'DisplayName', 'NO PI + 800nm CLC'); 
% plot(E(2:end), delay01, 'b--', 'DisplayName', 'Coulomb Delay l=0, Z=1'); 
% plot(E(2:end), delay02, 'r--', 'DisplayName', 'Coulomb Delay l=0, Z=2'); 
% % plot(E(2:end), delay03, 'g--', 'DisplayName', 'Coulomb Delay l=0, Z=3'); 
% xlabel('photoelectron energy (eV)'); ylabel('delay (as)'); 
% legend; goodplot(24); 
% % ylim([-300 2500])
% xlim([1 20])

figure; hold on; 
errorbar(sidebands*1240/810-9.55, NO_PI*((2*pi/(1240/800/E_AU) * T_AU)*1000/2/(2*pi)) + NO_CC_2400(2,:)*((2*pi/(1240/2400/E_AU) * T_AU)*1000/2/(2*pi)), NO_exp(:,2)'*((2*pi/(1240/800/E_AU) * T_AU)*1000/2/(2*pi)), 'ro', 'DisplayName', 'NO PI + 2.4um CLC'); 
errorbar(sidebands*1240/810-9.55, NO_streak*((2*pi/(1240/800/E_AU) * T_AU)*1000/2/(2*pi)), NO_exp(:,2)'*((2*pi/(1240/800/E_AU) * T_AU)*1000/2/(2*pi)), 'bo', 'DisplayName', 'NO PI + 800nm CLC'); 
% plot(E(2:end), delay01, 'b--', 'DisplayName', 'Coulomb Delay l=0, Z=1'); 
% plot(E(2:end), delay02, 'r--', 'DisplayName', 'Coulomb Delay l=0, Z=2'); 
% plot(E(2:end), delay03, 'g--', 'DisplayName', 'Coulomb Delay l=0, Z=3'); 
xlabel('photoelectron energy (eV)'); ylabel('delay (as)'); 
legend; goodplot(24); 
% ylim([-300 2500])
xlim([1 20])
