% plot for FERMI proposal demonstrating state resolution in CO2

subplot_width = 1.5; 
subplot_height = 0.2; 
hoffset = 0.05; 

% figure('Units', 'inches', 'Position', [1 1 1.2*subplot_width 6*subplot_height]); 
f = figure; 
hold on; 

pos1 = [0.1 0.3+hoffset 0.8 0.3];
subplot('Position',pos1); hold on; 
% yyaxis right; 
plot(E, XUV_only./sum(XUV_only), '-', 'Color', [0.4 0.4 0.4], ...
    'DisplayName', 'XUV only'); 
set(gca,'YTickLabel',[]);
% yyaxis left; 
plot(E, sum(sum(E_SpectraArray,2),3)./sum(E_SpectraArray(:)), 'k-', ...
    'DisplayName', 'IR + XUV (delay averaged)'); 
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
legend('Color', [1 1 1], 'EdgeColor', 'none'); 
goodplot(); 
ax1 = gca; 
AddHarmonicAxis(ax1, IP, IP_label, 810, 1); 

% pos2 = [0.1 0.3+hoffset 0.8 0.2];
% subplot('Position',pos2); hold on; 
% plot(E, sum(sum(E_SpectraArray,2),3), 'k-', ...
%     'DisplayName', '\tau averaged RABBITT'); 
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% % legend; 
% goodplot(); 
% ax2 = gca; 
% AddHarmonicAxis(ax2, IP, IP_label, 810, 0); 

pos3 = [0.1 0.1+hoffset 0.8 0.2];
subplot('Position',pos3); hold on; 
plot(E, abs(twoOmega_signal), 'b-', ...
    'DisplayName', '2\omega RABBITT amplitude'); 
set(gca,'YTickLabel',[]);
xlabel('electron kinetic energy (eV)'); 
legend('Color', [1 1 1], 'EdgeColor', 'none'); 
goodplot(); 
ax3 = gca; 
AddHarmonicAxis(ax3, IP, IP_label, 810, 0); 



