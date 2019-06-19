%% plot Coulomb scattering phase
% first run Jame's C_C_DelayPlot.m

hw = 1240/wavelength; 

% k = sqrt(2*((11:1:21)*hw-IP(1))/27.211); 
% k = logspace(-1, 0, 1000); 
% E = k.^2/2; 
% l = 1; 
sigma = unwrap(angle(igamma(l + 1 - 1j./k, 0))); 
delay = 24.2 * diff(sigma)./diff(E/E_AU); 

figure; 
plot(E, sigma); 
figure; 
plot(E(2:end), delay, 'o-'); 

%% get sigma_n
ref_n = 5; 
nstates = 6; 
hw = 1240/wavelength; 
sblist = [12, 14, 16]; 
% [sigmaE, sigmak, sigma_n] = coulombScatteringPhase(sblist, 810, IP); 
% [tmp1, tmp2, sigma_0] = coulombScatteringPhase(sblist, 810, IP(ref_n)); 
% sigma_sub = sigma_n - reshape(squeeze(repmat(sigma_0, [1, size(IP)]))',size(sigma_n)); 
enlist = cat(3, ...
             repmat((sblist-1)*hw, [nstates,1]) - repmat(IP(1:nstates), [numel(sblist),1])', ...
             repmat((sblist+1)*hw, [nstates,1]) - repmat(IP(1:nstates), [numel(sblist),1])', ...
             repmat((sblist)*hw,   [nstates,1]) - repmat(IP(1:nstates), [numel(sblist),1])');  
saveshape = size(enlist); 
enlist = enlist(:); 
enind  = zeros([saveshape]); 
enind = enind(:); 
for ii=1:numel(enlist)
    [tmp, enind(ii)] = min(abs(E-enlist(ii)/E_AU)); 
end
enind = reshape(enind, saveshape); 
minusind = enind(:,:,1); minusind = minusind(:); 
plusind  = enind(:,:,2); plusind  = plusind(:); 
nind     = enind(:,:,3); nind = nind(:); 
sigma_n = sigma(plusind) - sigma(minusind); 
sigmaE = E(nind); 

figure; plot(sigmaE, sigma_n, 'o-'); 
xlabel('photoelectron energy (a.u.)'); ylabel('Coulomb scattering phase'); 


%% reference RABBITT phase to v=ref_n
ref_n = 3; 
nstates = 6; 
hw = 1240/wavelength; 
sblist = [12, 14, 16]; 

measured_phases = [unwrap(sb12_phase(1:nstates,1))', ...
                   unwrap(sb14_phase(1:nstates,1))', ...
                   unwrap(sb16_phase(1:nstates,1))']; 
errors_phase =    [sb12_phase(1:nstates,2)', ...
                   sb14_phase(1:nstates,2)', ...
                   sb16_phase(1:nstates,2)']; 
% sub_phases = [sb12_param(:,2)'-sb12_param(ref_n,2), sb14_param(1:nstates,2)'-sb14_param(ref_n,2), sb16_param(1:nstates,2)'-sb16_param(ref_n,2)]; 
% measured_phases = [sb12_jackknife_phase(1:nstates)', sb14_jackknife_phase(1:nstates)', sb16_jackknife_phase(1:nstates)']; 
% sub_phases = [sb12_jackknife_phase(1:nstates)'-sb12_jackknife_phase(ref_n), sb14_jackknife_phase(1:nstates)'-sb14_jackknife_phase(ref_n), sb16_jackknife_phase(1:nstates)'-sb16_jackknife_phase(ref_n)]; 
measured_energy = [mean(sb12_param(:,2),3)', mean(sb14_param(:,2),3)', mean(sb16_param(:,2),3)']; 

% cc_phase = sub_phases - sigma_n; 

% figure; plot(sigmaE, sigma_n, 'o-'); 
% xlabel('photoelectron energy (a.u.)'); ylabel('Coulomb scattering phase'); 
figure; errorbar(measured_energy, measured_phases, errors_phase, 'o-'); 
xlabel('photoelectron energy (a.u.)'); ylabel('RABBITT phase');
% figure; plot(sigmaE, cc_phase, 'o-'); 
% xlabel('photoelectron energy (a.u.)'); ylabel('Continuum-continuum phase');


%% Plot CC phases
% % first run Jame's C_C_DelayPlot.m
% cc_time = cc_phase.*(T_L*1000/2/(2*pi)); 
% 
% figure; hold on; 
% plot(sigmaE, cc_time - 0.25*cc_time(end), 'o',...
%     'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
% % errorbar(sigmaE, cc_phase + 2*cc_phase(end), tmperror, 'o',...
% %     'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); 
% plot(fliplr(E/E_AU), CCP.*(T_L*1000/2/(2*pi)), 'LineWidth', 2); 
% plot(fliplr(E/E_AU), CCPA.*(T_L*1000/2/(2*pi)), 'LineWidth', 2); 
% plot(fliplr(E/E_AU), CCPAp.*(T_L*1000/2/(2*pi)), 'LineWidth', 2); 
% xlim([min(sigmaE) max(sigmaE)]); 
% xlabel('photoelectron energy (a.u.)'); ylabel('time delay'); 
% legend('H2 data cc phase', 'CCP', 'CCPA', 'CCPAp'); 

%% subtract CC models and Wigner Delays
% sigma and CC models should have been calculated with the same energy axis 

curve_CCP = fliplr(CCP(2:end)).*wavelength/2/0.3/(2*pi); 
curve_CCPA = fliplr(CCPA(2:end)).*wavelength/2/0.3/(2*pi); 
curve_CCPAp = fliplr(unwrap(CCPAp(2:end))).*wavelength/2/0.3/(2*pi); 
subcurve_CCP = delay + fliplr(CCP(2:end)).*wavelength/2/0.3/(2*pi); 
subcurve_CCPA = delay + fliplr(CCPA(2:end)).*wavelength/2/0.3/(2*pi); 
subcurve_CCPAp = delay + fliplr(unwrap(CCPAp(2:end))).*wavelength/2/0.3/(2*pi); 

% CLC delay formulas
tmpE = E/E_AU; 
tmpl = wavelength/0.0529; 
tmpw = 2*pi*137/tmpl; 
ge = 0.5772; 
alpha = 0; 
% tmpE = E; 
% w = 2*pi*0.29979*24.2/(wavelength*0.001)/1000; 

TCLC = -(1./(2*tmpE).^(3/2)).*(log(4*tmpE*tmpl)); 
TCLC_11 = -(1./(2*tmpE).^(3/2)).*(log(4*tmpE/tmpw) - ge + pi*tmpw./(8*tmpE)); 
TCLC_15 = -(1./(2*tmpE).^(3/2)).*(log(0.37*2*pi*tmpE/tmpw) - 1); 
TCLC_18 = -(1./(1 + alpha*tmpw/((2*tmpE).^(3/2)))).*(1./(2*tmpE).^(3/2)).*(log(4*tmpE/tmpw) - ge - 1 + 3*pi*tmpw./(4*(2*tmpE).^(3/2))); 
% TCLC_18 = -(1./(2*tmpE).^(3/2)).*(log(4*tmpE/tmpw) - 1 - ge + 3*pi*tmpw./(2*tmpE).^(3/2)/4); 

curve_TCLC_11 = TCLC_11(2:end)*24.2; 
curve_TCLC_15 = TCLC_15(2:end)*24.2; 
curve_TCLC_18 = TCLC_18(2:end)*24.2; 
subcurve_TCLC_11 = delay + TCLC_11(2:end)*24.2; 
subcurve_TCLC_15 = delay + TCLC_15(2:end)*24.2; 
subcurve_TCLC_18 = delay + TCLC_18(2:end)*24.2; 

figure; hold on; 

subplot(1,2,1); hold on; 
plot(E(2:end), curve_CCP, 'r', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'CC P'); 
plot(E(2:end), curve_CCPA, 'r', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'CC P+A'); 
plot(E(2:end), curve_CCPAp, 'r', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', 'CC P+Ap'); 
plot(E(2:end), curve_TCLC_11, 'b', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'CLC Ivanov'); 
plot(E(2:end), curve_TCLC_15, 'b', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'CLC Pazourek'); 
plot(E(2:end), curve_TCLC_18, 'b', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', 'CLC Serov'); 
plot(E(2:end), delay, 'k', 'DisplayName', 'WignerDelay'); 
% xlim([min(measured_energy), max(measured_energy)]); 
xlim([1.5 3])
ylim([curve_TCLC_18(1), 1000])
xlabel('photoelectron energy (eV)'); ylabel('time delay (as)'); 
legend; 
title('Delays')
hold off; 

subplot(1,2,2); hold on; 
plot(E(2:end), subcurve_CCP, 'r', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'CC P'); 
plot(E(2:end), subcurve_CCPA, 'r', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'CC P+A'); 
plot(E(2:end), subcurve_CCPAp, 'r', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', 'CC P+Ap'); 
plot(E(2:end), subcurve_TCLC_11, 'b', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'CLC Ivanov'); 
plot(E(2:end), subcurve_TCLC_15, 'b', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'CLC Pazourek'); 
plot(E(2:end), subcurve_TCLC_18, 'b', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', 'CLC Serov'); 
% plot(E(2:end), delay, 'k', 'DisplayName', 'WignerDelay'); 
% xlim([min(measured_energy), max(measured_energy)]); 
xlim([1.5 3])
% ylim([0, 10^4])
xlabel('photoelectron energy (eV)'); ylabel('time delay (as)'); 
legend; 
title('Wigner delay + CC delay')
hold off; 

hold off; 
%% 
measured_delays = measured_phases.*wavelength/2/0.3/(2*pi); 
errors_delay   = errors_phase.*wavelength/2/0.3/(2*pi); 

tmp_data  = cat(3, ...
            [measured_energy(1:6); measured_delays(1:6)], ...
            [measured_energy(7:12); measured_delays(7:12)], ...
            [measured_energy(13:18); measured_delays(13:18)]); 
tmp_error = cat(3, ...
            [measured_energy(1:6); errors_delay(1:6)-100], ...
            [measured_energy(7:12); errors_delay(7:12)], ...
            [measured_energy(13:18); errors_delay(13:18)]); 

theory_list = [subcurve_CCP; subcurve_CCPA; subcurve_CCPAp; ...
               subcurve_TCLC_11; subcurve_TCLC_15; subcurve_TCLC_18]; 
theory_name = ["CC P"; "CC P+A"; "CC P+Ap"; "CLC Ivanov"; "CLC Pazourek"; "CLC Serov"]; 

newfolder = strcat('/Users/annaliw/Documents/lab/plots/', date, '/');  
if ~exist(newfolder, 'dir')
   mkdir(newfolder)
end

% theory_list = theory_list(1); 

% for ii=1:1:size(theory_list,1)
%     fig = plotfun_compareToTheory(tmp_data, tmp_error, 3, [E(2:end); theory_list(ii,:)], theory_name(ii));         
%     saveas(fig, newfolder + regexprep(theory_name(ii),' ','_') + '.png'); 
% end

fig = plotfun_compareToTheory(tmp_data(:,2:end,:), tmp_error(:,2:end,:), 3, E(2:end), theory_list, theory_name); 


%%

figure; hold on; 

subplot(2, 3, 1); hold on; 
% model
plot(E(2:end), subcurve_CCP, 'LineWidth', 2, 'DisplayName', 'CCP'); 
%data
plot(measured_energy(1:6), measured_delays(1:6)-1, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(7:11), measured_delays(7:11)-5.5, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(12:16), measured_delays(12:16)-6, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
legend; 
xlim([min(measured_energy), max(measured_energy)]); 
xlabel('photoelectron energy (a.u.)'); ylabel('delay (fs)'); 
hold off; 

subplot(2, 3, 2); hold on; 
% model
plot(E(2:end), subcurve_CCPA, 'LineWidth', 2, 'DisplayName', 'CCPA'); 
%data
plot(measured_energy(1:6), measured_delays(1:6)-6.8, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(7:11), measured_delays(7:11)-6.65, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(12:16), measured_delays(12:16)-6.45, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
legend; 
xlim([min(measured_energy), max(measured_energy)]); 
xlabel('photoelectron energy (a.u.)'); ylabel('delay (fs)'); 
hold off; 

subplot(2, 3, 3); hold on; 
% model
plot(E(2:end), subcurve_CCPAp, 'LineWidth', 2, 'DisplayName', 'CCPAp');  
%data
plot(measured_energy(1:6), measured_delays(1:6)-3.2, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(7:11), measured_delays(7:11)-1.4, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(12:16), measured_delays(12:16)-0.9, 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
legend; 
xlim([min(measured_energy), max(measured_energy)]); 
xlabel('photoelectron energy (a.u.)'); ylabel('delay (fs)'); 
hold off; 

subplot(2, 3, 4); hold on; 
% model
plot(E(2:end), subcurve_TCLC_11, 'LineWidth', 2, 'DisplayName', 'TCLC 11');  
%data
plot(measured_energy(1:6), measured_delays(1:6), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(7:11), measured_delays(7:11), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(12:16), measured_delays(12:16), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
legend; 
xlim([min(measured_energy), max(measured_energy)]); 
xlabel('photoelectron energy (a.u.)'); ylabel('delay (fs)'); 
hold off; 

subplot(2, 3, 5); hold on; 
% model
plot(E(2:end), subcurve_TCLC_15, 'LineWidth', 2, 'DisplayName', 'TCLC 15');  
%data
plot(measured_energy(1:6), measured_delays(1:6), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(7:11), measured_delays(7:11), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(12:16), measured_delays(12:16), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
legend; 
xlim([min(measured_energy), max(measured_energy)]); 
xlabel('photoelectron energy (a.u.)'); ylabel('delay (fs)'); 
hold off; 

subplot(2, 3, 6); hold on; 
% model
plot(E(2:end), subcurve_TCLC_18, 'LineWidth', 2, 'DisplayName', 'TCLC 18');  
%data
plot(measured_energy(1:6), measured_delays(1:6), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(7:11), measured_delays(7:11), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
plot(measured_energy(12:16), measured_delays(12:16), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'data');
legend; 
xlim([min(measured_energy), max(measured_energy)]); 
xlabel('photoelectron energy (a.u.)'); ylabel('delay (fs)'); 
hold off; 

%% CO2 scratch
% load('CO2_20190414.mat')

measured_energy = [mean(X12_param(:,2,:),3), ...
                   mean(X14_param(:,2,:),3), ...
                   mean(X16_param(:,2,:),3), ...
                   mean(X18_param(:,2,:),3), ...
                   mean(A14_param(:,2,:),3), ...
                   mean(A16_param(:,2,:),3), ...
                   mean(A18_param(:,2,:),3), ...
                   mean(B14_param(:,2,:),3), ...
                   mean(B16_param(:,2,:),3), ...
                   mean(B18_param(:,2,:),3)]; 
measured_phase  = [X12_phase, X14_phase, X16_phase, X18_phase, ...
                   A14_phase, A16_phase, A18_phase, ...
                   B14_phase, B16_phase, B18_phase]; 
measured_error  = [X12_phase_std, X14_phase_std, X16_phase_std, X18_phase_std, ...
                   A14_phase_std, A16_phase_std, A18_phase_std, ...
                   B14_phase_std, B16_phase_std, B18_phase_std]; 

measured_phase = unwrap(measured_phase); 
measured_delay = measured_phase.*wavelength/2/0.3/(2*pi); 

figure; hold on; 
errorbar(measured_energy(1:4), measured_phase(1:4), measured_error(1:4), '-o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k', ...
    'MarkerSize', 4, 'LineWidth', 2, 'DisplayName', 'X'); 
errorbar(measured_energy(5:7), measured_phase(5:7), measured_error(5:7), '-o', ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Color', 'r', ...
    'MarkerSize', 4, 'LineWidth', 2, 'DisplayName', 'A'); 
errorbar(measured_energy(8:10), measured_phase(8:10), measured_error(8:10), '-o', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'Color', 'b', ...
    'MarkerSize', 4, 'LineWidth', 2, 'DisplayName', 'B'); 
legend; 
xlabel('photoelectron energy (eV)'); ylabel('phase'); 
hold off; 

