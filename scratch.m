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
% first run Jame's C_C_DelayPlot.m (with the correct l, Z values)

l=0; Z=1; % make sure these match the CCdelay calculations
[sigma, delay] = coulombScatteringPhase(l, Z, E); 

curve_CCP = fliplr(CCP(2:end)).*(T_L*1000/2/(2*pi)); 
curve_CCPA = fliplr(CCPA(2:end)).*(T_L*1000/2/(2*pi)); 
curve_CCPAp = fliplr(unwrap(CCPAp(2:end))).*(T_L*1000/2/(2*pi)); 
subcurve_CCP = delay + fliplr(CCP(2:end)).*(T_L*1000/2/(2*pi)); 
subcurve_CCPA = delay + fliplr(CCPA(2:end)).*(T_L*1000/2/(2*pi)); 
subcurve_CCPAp = delay + fliplr(unwrap(CCPAp(2:end))).*(T_L*1000/2/(2*pi)); 


% CLC delay formulas
tmpE = E/E_AU; 
tmpl = wavelength/0.0529; 
tmpw = 2*pi*137/tmpl; 
ge = 0.5772; 
alpha = 0; 
% tmpE = E; 
% w = 2*pi*0.29979*24.2/(wavelength*0.001)/1000; 

TCLC_11 = Ivanov_curve(l, Z, E); 
TCLC_15 = -(1./(2*tmpE).^(3/2)).*(log(0.37*2*pi*tmpE/tmpw) - 1); 
TCLC_18 = Serov_curve(Z, E); 

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
ylim([-2000, 1000])
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
%% proccess data
ref_n = 3; 
nstates = 6; 
phase_data = cat(3, phase_SB12_slope, phase_SB14_slope, phase_SB16_noslope); 
mean_data = cat(3, mean(paramout_SB12_slope(:,1:4,:),3), mean(paramout_SB14_slope(:,1:4,:),3), mean(paramout_SB16_noslope,3)); 

measured_phases = [unwrap(phase_data(1:nstates,1,1))', ...
                   unwrap(phase_data(1:nstates,1,2))', ...
                   unwrap(phase_data(1:nstates,1,3))']; 
errors_phase =    [phase_data(1:nstates,2,1)', ...
                   phase_data(1:nstates,2,2)', ...
                   phase_data(1:nstates,2,3)'];
               
measured_delays = measured_phases.*(T_L*1000/2/(2*pi)); 
errors_delay = errors_phase.*(T_L*1000/2/(2*pi)); 

measured_energy = [mean_data(:,2,1)', mean_data(:,2,2)', mean_data(:,2,3)']; 


tmp_data  = cat(3, ...
            [measured_energy(1:6); measured_delays(1:6)], ...
            [measured_energy(7:12); measured_delays(7:12)], ...
            [measured_energy(13:18); measured_delays(13:18)]); 
tmp_error = cat(3, ...
            [measured_energy(1:6); errors_delay(1:6)], ...
            [measured_energy(7:12); errors_delay(7:12)], ...
            [measured_energy(13:18); errors_delay(13:18)]); 
        
%%
% theory_list = [subcurve_CCP; subcurve_CCPA; subcurve_CCPAp; ...
%                subcurve_TCLC_11; subcurve_TCLC_15; subcurve_TCLC_18]; 
% theory_name = ["CC P"; "CC P+A"; "CC P+Ap"; "CLC Ivanov"; "CLC Pazourek"; "CLC Serov"]; 
% 
% newfolder = strcat('/Users/annaliw/Documents/lab/plots/', date, '/');  
% if ~exist(newfolder, 'dir')
%    mkdir(newfolder)
% end
% 
% % theory_list = theory_list(1); 
% 
% % for ii=1:1:size(theory_list,1)
% %     fig = plotfun_compareToTheory(tmp_data, tmp_error, 3, [E(2:end); theory_list(ii,:)], theory_name(ii));         
% %     saveas(fig, newfolder + regexprep(theory_name(ii),' ','_') + '.png'); 
% % end
% 
% fig = plotfun_compareToTheory(tmp_data(:,2:end,:), tmp_error(:,2:end,:), 3, E(2:end), theory_list, theory_name); 

%%
[sigma0, delay0] = coulombScatteringPhase(0, 1, E); 
[sigma1, delay1] = coulombScatteringPhase(1, 1, E); 
% plotfun_compareToTheory(tmp_data(:,2:end,:), tmp_error(:,2:end,:), 3, E(2:end), [Serov_curve(1,E(2:end))*24.2+delay; Serov_curve(2,E(2:end))*24.2+delay], ["Z=1", "Z=2"]); 
plotfun_compareToTheory(tmp_data(:,2:end,:), tmp_error(:,2:end,:), 3, E(2:end), ...
    [Ivanov_curve(0, 1,E(2:end))*24.2+delay0; Ivanov_curve(1,1,E(2:end))*24.2+delay1; ...
     Ivanov_curve(0, 2,E(2:end))*24.2+delay1; Ivanov_curve(1,1,E(2:end))*24.2+delay0], ...
     ["l=0 Wigner delay, l=0 tIR", "l=1 Wigner delay, l=1 tIR", ...
      "l=1 Wigner Delay, l=0 tIR", "l=0 Wigner delay, l=1 tIR"]); 

