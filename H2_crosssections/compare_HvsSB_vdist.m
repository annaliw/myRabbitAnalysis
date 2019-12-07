% region = [0.05 1.6]; % harmonic 11
region = [3.2 4.8]; % harmonic 13
tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 
[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, ...
    E(start:stop), abs(XUV_only(start:stop)), XUV_only(start:stop), 1, 0); 
% save as labeled variables
paramout_gauss_H13 = paramout_gauss; 

text_size = 14; 
line_weight = 2; 
f1 = figure; hold on; 
axH = gca; 
axH_pos = axH.Position; 
        axH.XColor = 'k'; 
        axH.XLabel.String = 'Photoelectron Energy (eV)'; 
        axH.YLabel.String = 'Amplitude (arbitrary units)';
        axH.FontSize = text_size; 
        axH.YColor = 'k'; 
        axH.LineWidth = line_weight*0.5;
line(E(start:stop)-13*1240/810+IP(1), XUV_only(start:stop)./sum(XUV_only(start:stop)), ...
    'Parent', axH, 'Color', 'k', 'DisplayName', 'H13'); 


region = [1.55 3.05]; % sideband 12
tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 
[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, ...
    E(start:stop), abs(twoOmega_signal(start:stop)), twoOmega_signal(start:stop), 1, 0); 
% save as labeled variables
paramout_gauss_SB12 = paramout_gauss; 

axSB = axes(f1, 'Position', axH_pos, ...
            'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
            'Color', 'none', 'LineWidth', line_weight*0.5); 
line(E(start:stop)-12*1240/810+IP(1), abs(twoOmega_signal(start:stop))./sum(abs(twoOmega_signal(start:stop))), ...
    'Parent', axSB, 'Color', 'b', 'DisplayName', 'SB12'); 
axis off; 

%% plot cross sections

v = 0:1:5; 

figure; hold on; 
plot(v, paramout_gauss_H13(:,1).*paramout_gauss_H13(:,3), 'o-', ...
    'DisplayName', 'H13'); 
plot(v, paramout_gauss_SB12(:,1).*paramout_gauss_SB12(:,3), 'o-', ...
    'DisplayName', 'SB12'); 
xlabel('v state'); 
ylabel('cross section (arb)'); 
legend; 
goodplot(); 

figure; hold on; 
plot(v, paramout_gauss_H13(:,1), 'o-', ...
    'DisplayName', 'H13'); 
plot(v, paramout_gauss_SB12(:,1), 'o-', ...
    'DisplayName', 'SB12'); 
xlabel('v state'); 
ylabel('amplitue (arb)'); 
legend; 
goodplot(); 



