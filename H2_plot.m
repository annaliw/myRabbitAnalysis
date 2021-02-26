% load H2_plot_20200119.mat
%%
ax_font_size = 12; 
axv_font_size = 10; 
legend_font_size = 12; 
label_font_size = 10; 
%%
xdata = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
phase_data = squeeze(reshape(SB_delay_data(2,:,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, :,:), [1 15]));

% % use the mean of the measured H2 values
% xdata = squeeze(reshape(mean(SB_delay_data(1,2:end,:),2), [1 3])); 
% phase_data = squeeze(reshape(mean(SB_delay_data(2,2:end,:),2), [1 3])); 
% phase_error = squeeze(reshape(mean(SB_delay_error(2,2:end,:),2), [1 3])); 

% subtract out XUV contribution by referencing to Ar data and adding in
% provided Ar TDSE values
plot_data = phase_data - reshape(repmat(XUV_delay(1:3)-XUV_delay(1), [1 5])', [1 15]) + 140; 
% plot_data = phase_data - reshape(repmat(XUV_delay(1:3)-XUV_delay(1), [1 5])', [1 15])./2 - interp1(E_CC, CCPAp, xdata)*(T_L*1000/2/(2*pi)); 
plot_error = phase_error; 


figure; hold on; 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.96]);

ax1 = subplot(2, 3, 1); hold on; 
errorbar(xdata(1:5), plot_data(1:5), plot_error(1:5), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(1), -H2TDSE_810_140.t(1), H2TDSE_810_140.err(1), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(1), -H2TDSE_810_145.t(1), H2TDSE_810_145.err(1), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(1), -H2TDSE_810_150.t(1), H2TDSE_810_150.err(1), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_145, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_150, 'r--', 'HandleVisibility', 'off'); 
ylabel('delay (as)'); 
xlim([1.4 3])
ylim([100 180]); 
window = 60; 

text(1.4+1, 100+10, 'Sideband 12', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 

box 'on'; 
ax1.FontSize = ax_font_size; 
ax1.FontWeight = 'bold'; 
ax1.LineWidth = 1; 
ax1.Color = 'none'; 
ax1.XMinorTick = 'on'; 
ax1.XAxisLocation = 'bottom'; 
ax1.YMinorTick = 'on'; 

% add v state axis
POS = ax1.Position; 
ax1v = axes('Position',POS);
ax1v.XAxisLocation = 'top'; 
ax1v.Color = 'none'; 
ax1v.XTick = xdata(1:5); 
ax1v.XTickLabel = split(int2str(fliplr(1:5))); 
% xlabel('\nu'); 
ax1v.XGrid = 'on'; 
ax1v.GridAlpha = 0.7; 
ax1v.YTick = []; 
ax1v.FontSize = axv_font_size; 
% ax1v.FontWeight = 'bold'; 

linkaxes([ax1,ax1v],'x'); 
xlim([1.4 3])
ylim([100 180]); 
uistack(ax1, 'top'); 


ax2 = subplot(2, 3, 2); hold on; 
errorbar(xdata(6:10), plot_data(6:10), plot_error(6:10), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(2), -H2TDSE_810_140.t(2), H2TDSE_810_140.err(2), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(2), -H2TDSE_810_145.t(2), H2TDSE_810_145.err(2), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(2), -H2TDSE_810_150.t(2), H2TDSE_810_150.err(2), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_145, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_150, 'r--', 'HandleVisibility', 'off'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([4.4 6.2]); 
ylim([120-window 120+window]); 

text(4.4+1, 120-40+10, 'Sideband 14', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 
% annotation('textbox',[0.45 0.15 0.5 0.5],'String','Sideband 14','FitBoxToText','on');

box 'on'; 
ax2.FontSize = ax_font_size; 
ax2.FontWeight = 'bold'; 
ax2.LineWidth = 1; 
ax2.Color = 'none'; 
ax2.XMinorTick = 'on'; 
ax2.XAxisLocation = 'bottom'; 
ax2.YMinorTick = 'on'; 

% add v state axis
POS = ax2.Position; 
ax2v = axes('Position',POS);
ax2v.XAxisLocation = 'top'; 
ax2v.Color = 'none'; 
ax2v.XTick = xdata(6:10); 
ax2v.XTickLabel = split(int2str(fliplr(1:5))); 
% xlabel('\nu'); 
ax2v.XGrid = 'on'; 
ax2v.GridAlpha = 0.7; 
ax2v.YTick = []; 
ax2v.FontSize = axv_font_size; 
% ax2.FontWeight = 'bold'; 

linkaxes([ax2,ax2v],'x'); 
xlim([4.4 6.2]); 
ylim([120-window 120+window]); 
uistack(ax2, 'top'); 


ax3 = subplot(2, 3, 3); hold on; 
errorbar(xdata(11:15), plot_data(11:15), plot_error(11:15), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(3), -H2TDSE_810_140.t(3), H2TDSE_810_140.err(3), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(3), -H2TDSE_810_145.t(3), H2TDSE_810_145.err(3), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(3), -H2TDSE_810_150.t(3), H2TDSE_810_150.err(3), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_145, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_150, 'r--', 'HandleVisibility', 'off'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([7.5 9.2]); 
ylim([55-window 55+window]); 

text(7.5+1, 55-40+10, 'Sideband 16', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 

box 'on'; 
ax3.FontSize = ax_font_size; 
ax3.FontWeight = 'bold'; 
ax3.LineWidth = 1; 
ax3.Color = 'none'; 
ax3.XMinorTick = 'on'; 
ax3.XAxisLocation = 'bottom'; 
ax3.YMinorTick = 'on'; 

% add v state axis
POS = ax3.Position; 
ax3v = axes('Position',POS);
ax3v.XAxisLocation = 'top'; 
ax3v.Color = 'none'; 
ax3v.XTick = xdata(11:15); 
ax3v.XTickLabel = split(int2str(fliplr(1:5))); 
xlabel('\nu'); 
ax3v.XGrid = 'on'; 
ax3v.GridAlpha = 0.7; 
ax3v.YTick = []; 
ax3v.FontSize = axv_font_size; 
% ax3.FontWeight = 'bold'; 

linkaxes([ax3,ax3v],'x'); 
xlim([7.5 9.2]); 
ylim([55-window 55+window]); 
uistack(ax3, 'top'); 


ax4 = subplot(2, 3, [4 5 6]); hold on; 
errorbar(xdata, plot_data, plot_error, 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee, -H2TDSE_810_140.t, H2TDSE_810_140.err, ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee, -H2TDSE_810_145.t, H2TDSE_810_145.err, ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee, -H2TDSE_810_150.t, H2TDSE_810_150.err, ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_145, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_150, 'r--', 'HandleVisibility', 'off'); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
l = legend('Location', 'southeast'); 
xlim([1 10]); 
ylim([0 250]); 

% text(3.5, -200, 'Full View', 'FontSize', 14, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center'); 

box 'on'; 
ax4.FontSize = ax_font_size; 
ax4.FontWeight = 'bold'; 
ax4.LineWidth = 1; 
ax4.Color = 'none'; 

ax4.XTick = 2:1:10;
ax4.XTickLabel = split(int2str(2:1:10)); 
ax4.XMinorTick = 'on'; 
ax4.XAxisLocation = 'top'; 

ax4.YTick = 0:50:250; 
ax4.YTickLabel = split(int2str(0:50:250)); 
ax4.YMinorTick = 'on'; 

% add v state axis
POS = ax4.Position; 
ax4v = axes('Position',POS);
ax4v.XAxisLocation = 'bottom'; 
ax4v.Color = 'none'; 
ax4v.XTick = xdata; 
ax4v.XTickLabel = split(int2str([fliplr(1:5) fliplr(1:5) fliplr(1:5)])); 
xlabel('\nu'); 
ax4v.XGrid = 'on'; 
ax4v.GridAlpha = 0.7; 
ax4v.YTick = []; 
ax4v.FontSize = axv_font_size; 
% ax4.FontWeight = 'bold'; 

linkaxes([ax4,ax4v],'x'); 
xlim([1 10]); 
ylim([-0 250]); 
uistack(ax4, 'top'); set(l, 'color', 'white', 'FontSize', legend_font_size); 

% make large markers and lines
set(findall(gcf, 'Type', 'Line'), 'MarkerSize', 10, 'LineWidth', 2); 
set(findall(gcf, 'Type', 'ErrorBar'), 'MarkerSize', 10, 'LineWidth', 2); 
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [6 8]);
set(gcf,'PaperPosition',[0.5 0.5 7 7]);
set(gcf,'PaperPositionMode','Manual');
