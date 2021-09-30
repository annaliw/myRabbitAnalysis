% load H2_plot_20200119.mat
%%
ax_font_size = 12; 
axv_font_size = 10; 
legend_font_size = 12; 
label_font_size = 10; 
%%
% vx = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
xdata = reshape(SB_delay_data(1,:,:), [1 15]); xdata(11:15) = xdata(11:15)+1.5+0.2286; 
vx = [fliplr(xdata(1:5)) fliplr(xdata(6:10)) xdata(11) fliplr(xdata(12:15))]; 
phase_data = squeeze(reshape(SB_delay_data(2,:,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2,:,:), [1 15]));

% % use the mean of the measured H2 values
% xdata = squeeze(reshape(mean(SB_delay_data(1,2:end,:),2), [1 3])); 
% phase_data = squeeze(reshape(mean(SB_delay_data(2,2:end,:),2), [1 3])); 
% phase_error = squeeze(reshape(mean(SB_delay_error(2,2:end,:),2), [1 3])); 

% subtract out XUV contribution by referencing to Ar data and adding in
% provided Ar TDSE values
% offset12 = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
%     - ((mean(unwrap(H2_0V_SB18_phase(2:6,1)))+2*pi)*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 
% offset14 = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
%     - (mean(unwrap(H2_3V_SB18_phase(:,1)))*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 
% offset16 = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
%     - (mean(unwrap(H2_5V_SB18_phase(:,1)))*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 
offset = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
    - ((mean(unwrap(H2_0V_SB18_phase(2:5)))+mean(unwrap(H2_3V_SB18_phase(2:5)))+mean(unwrap(H2_5V_SB18_phase(2:5)))+2*pi)/3*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 
% offset = - (2*pi*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 

plot_data = phase_data - reshape(repmat(XUV_delay(1:3), [1 5])', [1 15])+offset; 
% plot_data = phase_data - reshape(repmat(XUV_delay(1:3)-XUV_delay(1), [1 5])', [1 15])./2 - interp1(E_CC, CCPAp, xdata)*(T_L*1000/2/(2*pi)); 
plot_error = phase_error; 

SB_shift_err = [0.0354 0.0223 0.0224]*(T_L*1000/4/pi); 
SB_shade_err = [interp1(xdata(1:5), plot_data(1:5), E_CC, 'spline'); 
                interp1(xdata(6:10), plot_data(6:10), E_CC, 'spline'); 
                interp1(xdata(11:15), plot_data(11:15), E_CC, 'spline')]; 

figure; hold on; 
plotpos = [0, 0, 20, 14]/2.75; 
set(gcf, 'units', 'inch', 'position', plotpos);

ax1 = subplot(2, 3, 1); hold on; 
set(ax1, 'position', [ax1.Position(1) ax1.Position(2)+ax1.Position(4)*0.15 ax1.Position(3) ax1.Position(4)*0.85]); 
% patch('XData', [E_CC(332:454) fliplr(E_CC(332:454))], ...
%     'YData', [SB_shade_err(1,332:454)+SB_shift_err(1)/2 fliplr(SB_shade_err(1,332:454))-SB_shift_err(1)/2], ...
%     'FaceColor', [1 1 1]*0.85, 'EdgeColor', 'none');
errorbar(xdata(1:5), plot_data(1:5), plot_error(1:5), 'ko', 'DisplayName', 'H2 measurement', ...
    'MarkerSize', 4, 'LineWidth', 2); 
errorbar(H2TDSE_810_140.x_Ee(1), -H2TDSE_810_140.t(1), H2TDSE_810_140.err(1), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_145.x_Ee(1), -H2TDSE_810_145.t(1), H2TDSE_810_145.err(1), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_150.x_Ee(1), -H2TDSE_810_150.t(1), H2TDSE_810_150.err(1), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50', ...
    'MarkerSize', 8, 'LineWidth', 2); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
plot(E_CC, -interp_145, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
plot(E_CC, -interp_150, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
ylabel('delay (as)'); 
xlim([1.4 3])
ylim([50 160]); 

text(1.4+1, 50+10, 'Sideband 12', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 

box 'on'; 
ax1.FontSize = ax_font_size; 
ax1.FontWeight = 'bold'; 
ax1.LineWidth = 1; 
ax1.Color = 'none'; 
ax1.XMinorTick = 'on'; 
ax1.XAxisLocation = 'bottom'; 
ax1.YMinorTick = 'on'; 
ax1.YTick = 50:20:160; 
ax1.YTickLabel = 50:20:160; 

% add v state axis
POS = ax1.Position; 
ax1v = axes('Position',POS);
ax1v.XAxisLocation = 'top'; 
ax1v.Color = 'none'; 
ax1v.XTick = vx(1:5); 
ax1v.XTickLabel = split(int2str(fliplr(1:5))); 
% xlabel('\nu'); 
ax1v.XGrid = 'on'; 
ax1v.GridAlpha = 0.7; 
ax1v.YTick = []; 
ax1v.FontSize = axv_font_size; 
% ax1v.FontWeight = 'bold'; 

linkaxes([ax1,ax1v],'x'); 
xlim([1.4 3])
ylim([50 160]); 
uistack(ax1, 'top'); 


ax2 = subplot(2, 3, 2); hold on; 
set(ax2, 'position', [ax2.Position(1) ax2.Position(2)+ax2.Position(4)*0.15 ax2.Position(3) ax2.Position(4)*0.85])
% patch('XData', [E_CC(582:656) fliplr(E_CC(582:656))], ...
%     'YData', [SB_shade_err(2,582:656)+SB_shift_err(2)/2 fliplr(SB_shade_err(2,582:656))-SB_shift_err(2)/2], ...
%     'FaceColor', [1 1 1]*0.85, 'EdgeColor', 'none');
errorbar(xdata(6:10), plot_data(6:10), plot_error(6:10), 'ko', 'DisplayName', 'H2 measurement', ...
    'MarkerSize', 4, 'LineWidth', 2); 
errorbar(H2TDSE_810_140.x_Ee(2), -H2TDSE_810_140.t(2), H2TDSE_810_140.err(2), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_145.x_Ee(2), -H2TDSE_810_145.t(2), H2TDSE_810_145.err(2), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_150.x_Ee(2), -H2TDSE_810_150.t(2), H2TDSE_810_150.err(2), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50', ...
    'MarkerSize', 8, 'LineWidth', 2); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r:', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_145, 'r:', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_150, 'r:', 'HandleVisibility', 'off'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([4.4 6.2]); 
ylim([90-window 90+window]); 

text(4.4+1, 90-window+10, 'Sideband 14', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 
% annotation('textbox',[0.45 0.15 0.5 0.5],'String','Sideband 14','FitBoxToText','on');

box 'on'; 
% ax2.YTickMode = 'Manual'; 
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
ax2v.XTick = vx(6:10); 
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
set(ax3, 'position', [ax3.Position(1) ax3.Position(2)+ax3.Position(4)*0.15 ax3.Position(3) ax3.Position(4)*0.85])
% patch('XData', [E_CC(750:809) fliplr(E_CC(750:809))], ...
%     'YData', [SB_shade_err(3,750:809)+SB_shift_err(3)/2 fliplr(SB_shade_err(3,750:809))-SB_shift_err(3)/2], ...
%     'FaceColor', [1 1 1]*0.85, 'EdgeColor', 'none');
errorbar(xdata(11:15), plot_data(11:15), plot_error(11:15), 'ko', 'DisplayName', 'H2 measurement', ...
    'MarkerSize', 4, 'LineWidth', 2); 
errorbar(H2TDSE_810_140.x_Ee(3), -H2TDSE_810_140.t(3), H2TDSE_810_140.err(3), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_145.x_Ee(3), -H2TDSE_810_145.t(3), H2TDSE_810_145.err(3), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_150.x_Ee(3), -H2TDSE_810_150.t(3), H2TDSE_810_150.err(3), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50', ...
    'MarkerSize', 8, 'LineWidth', 2); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
plot(E_CC, -interp_145, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
plot(E_CC, -interp_150, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([7.5 9.2]); 
ylim([55-window 55+window]); 

text(7.5+1, 55-window+10, 'Sideband 16', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 

box 'on'; 
% ax3.YTickMode = 'Manual'; 
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
ax3v.XTick = vx(11:15); 
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
% patch('XData', [E_CC(332:454) fliplr(E_CC(332:454))], ...
%     'YData', [SB_shade_err(1,332:454)+SB_shift_err(1)/2 fliplr(SB_shade_err(1,332:454))-SB_shift_err(1)/2], ...
%     'FaceColor', [1 1 1]*0.85, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% patch('XData', [E_CC(582:656) fliplr(E_CC(582:656))], ...
%     'YData', [SB_shade_err(2,582:656)+SB_shift_err(2)/2 fliplr(SB_shade_err(2,582:656))-SB_shift_err(2)/2], ...
%     'FaceColor', [1 1 1]*0.85, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% patch('XData', [E_CC(750:809) fliplr(E_CC(750:809))], ...
%     'YData', [SB_shade_err(3,750:809)+SB_shift_err(3)/2 fliplr(SB_shade_err(3,750:809))-SB_shift_err(3)/2], ...
%     'FaceColor', [1 1 1]*0.85, 'EdgeColor', 'none', 'HandleVisibility', 'off');
errorbar(xdata, plot_data, plot_error, 'ko', 'DisplayName', 'H2 measurement', ...
    'MarkerSize', 4, 'LineWidth', 2); 
errorbar(H2TDSE_810_140.x_Ee, -H2TDSE_810_140.t, H2TDSE_810_140.err, ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_145.x_Ee, -H2TDSE_810_145.t, H2TDSE_810_145.err, ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45', ...
    'MarkerSize', 8, 'LineWidth', 2); 
errorbar(H2TDSE_810_150.x_Ee, -H2TDSE_810_150.t, H2TDSE_810_150.err, ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50', ...
    'MarkerSize', 8, 'LineWidth', 2); 
% plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
plot(E_CC, -interp_140, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
plot(E_CC, -interp_145, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
plot(E_CC, -interp_150, 'r:', 'HandleVisibility', 'off', 'LineWidth', 1); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
% l = legend('Location', 'northeast'); 
xlim([1 10]); 
ylim([20 170]); 

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

% ax4.YTickMode = 'Manual'; 
ax4.YTick = 0:50:250; 
ax4.YTickLabel = split(int2str(0:50:250)); 
ax4.YMinorTick = 'on'; 

% add v state axis
POS = ax4.Position; 
ax4v = axes('Position',POS);
ax4v.XAxisLocation = 'bottom'; 
ax4v.Color = 'none'; 
ax4v.XTick = vx; 
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
uistack(ax4, 'top'); 
% set(l, 'color', 'white', 'FontSize', legend_font_size); 

% make large markers and lines
% set(findall(gcf, 'Type', 'Line'), 'MarkerSize', 4, 'LineWidth', 2); 
% set(findall(gcf, 'Type', 'ErrorBar'), 'MarkerSize', 4, 'LineWidth', 1); 
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [8.5 11]);
set(gcf,'PaperPosition',[0.5 0.5 5 7]);
set(gcf,'PaperPositionMode','Manual');

% save
set(gcf, 'units', 'inch', 'position', plotpos);
set(gcf,'PaperUnits','inches','PaperPosition',plotpos)
% saveas(gcf,'/Users/annawang/Box/writing/H2/paper/figures/H2_comparewtheory.eps', 'epsc')


%% raw spectra summary plot (carpet on bottom panel)
subplot_width = 1.5; 
subplot_height = 0.2; 
hoffset = 0.3; 

xr = [1 9.5]; 

tolerance = abs(Ebins(2)-Ebins(1)); 

region_12 = [1.55 2.8]; % sideband 12
start_12 = find(abs(Ebins-region_12(1))<tolerance, 1, 'last'); 
stop_12 = find(abs(Ebins-region_12(2))<tolerance, 1, 'first'); 
region_14 = [4.65 5.85]; % sideband 14
start_14 = find(abs(Ebins-region_14(1))<tolerance, 1, 'last'); 
stop_14 = find(abs(Ebins-region_14(2))<tolerance, 1, 'first'); 
region_16 = [5.9 7]+0.25; % SB16 5V; 
start_16 = find(abs(Ebins-region_16(1))<tolerance, 1, 'last'); 
stop_16 = find(abs(Ebins-region_16(2))<tolerance, 1, 'first'); 

% figure('Units', 'inches', 'Position', [1 1 1.2*subplot_width 6*subplot_height]); 
f = figure; hold on; 
plotpos = [0, 0, 12, 7]/1.3750;
% plotpos = [0, 0, 20, 14]/2.75; 
set(gcf, 'units', 'inch', 'position', plotpos);

pos1 = [0.1 0.2+hoffset 0.8 0.2];
ax1 = subplot('Position',pos1); hold on; 
plot(Ebins(start_12:stop_12), sum(sum(ESpectra_0V(start_12:stop_12,:,:),2),3)./sum(sum(sum(ESpectra_0V(start_12:stop_12,:,:),1),2),3), ...
    'k-.', 'DisplayName', 'time averaged yield', 'LineWidth', 2); 
plot(ones(20)*3.81, (-4:15)*0.02, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');  
plot(Ebins(start_14:stop_14), sum(sum(ESpectra_3V(start_14:stop_14,:,:),2),3)./sum(sum(sum(ESpectra_3V(start_14:stop_14,:,:),1),2),3), ...
    'k-.', 'HandleVisibility', 'off', 'LineWidth', 2); 
plot(ones(20)*6.65, (-4:15)*0.02, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');  
plot(Ebins(start_16:stop_16)+1.5+0.2286, sum(sum(ESpectra_5V(start_16:stop_16,:,:),2),3)./sum(sum(sum(ESpectra_5V(start_16:stop_16,:,:),1),2),3), ...
    'k-.', 'HandleVisibility', 'off', 'LineWidth', 2); 
ylabel({'Yield';'(Norm.)'}, 'FontSize', 9, 'FontWeight', 'bold'); 
xlim([xr(1) xr(2)]); 
% ylim([0.005 0.03])
ylim([0 max(abs(twoOmega_0V)./sum(abs(twoOmega_0V(start_12:stop_12))))*1.5]); 
ax1.YAxisLocation = 'left';
ax1.Color = 'none'; 
ax1.FontWeight = 'bold'; 
ax1.XTick = 2:1:10;
ax1.YTick = [0.005 max(abs(twoOmega_0V)./sum(abs(twoOmega_0V(start_12:stop_12))))*1]; 
ax1.YTickLabel = ["0", "1"]; 
ax1.XTickLabel = []; 
ax1.XAxisLocation = 'bottom'; 
ax1.XMinorTick = 'on'; 
box 'on'; 
% set(gca,'YTickLabel',[]);
% 2w plots
% pos1_2w = pos1;
% ax1_2w = subplot('Position',pos1); hold on; 
ax1_2w = axes('Position', pos1, ...
            'XAxisLocation', 'bottom', 'YAxisLocation', 'right', ...
            'Color', 'none'); hold on; 
plot(Ebins(start_12:stop_12), abs(twoOmega_0V(start_12:stop_12))./sum(abs(twoOmega_0V(start_12:stop_12))), ...
    'b-', 'DisplayName', '2\omega amplitude', 'LineWidth', 1.5); 
plot(ones(20)*3.81, (-4:15)*0.02, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off'); 
plot(Ebins(start_14:stop_14), abs(twoOmega_3V(start_14:stop_14))./sum(abs(twoOmega_3V(start_14:stop_14))), ...
    'b-', 'HandleVisibility', 'off', 'LineWidth', 1.5); 
plot(Ebins(start_16:stop_16)+1.5+0.2286, abs(twoOmega_5V(start_16:stop_16))./sum(abs(twoOmega_5V(start_16:stop_16))), ...
    'b-', 'HandleVisibility', 'off', 'LineWidth', 1.5); 
plot(ones(20)*6.65, (-4:15)*0.02, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
ylabel('2\omega Amp.', 'FontSize', 9, 'FontWeight', 'bold'); 
% goodplot(20); 
xlim([xr(1) xr(2)]); 
% ylim([0.005 0.03])
ylim([0 max(abs(twoOmega_0V)./sum(abs(twoOmega_0V(start_12:stop_12))))*1.5]); 
box 'on'; 
ax1_2w.Color = 'none'; 
ax1_2w.YColor = 'b'; 
ax1_2w.XColor = 'none'; 
ax1_2w.FontWeight = 'bold'; 
% % ax1.YTick = max(abs(twoOmega_0V)./sum(abs(twoOmega_0V(start_12:stop_12))))*0.5;
ax1_2w.YTick = [0.005 max(abs(twoOmega_0V)./sum(abs(twoOmega_0V(start_12:stop_12))))*1]; 
% ax1.YRuler.TickLabelGapOffset = -12; 
% ax1_2w.YRuler.TickLabelGapOffset = -15; 
ax1_2w.YTickLabel = ["0", "0.03"]; 
ax1_2w.XTick = 2:1:10; 
ax1_2w.XTickLabel = []; 
ax1_2w.XMinorTick = 'on'; 
ax1_2w.XAxisLocation = 'bottom'; 
ax1_2w.YAxisLocation = 'right'; 

% add v state axis
POS = ax1.Position; 
ax1v = axes('Position',POS);
ax1v.XAxisLocation = 'top'; 
ax1v.Color = 'none'; 
ax1v.XTick = [fliplr(12*1240/810-IP) fliplr(14*1240/810-IP) fliplr(16*1240/810-IP)]; 
ax1v.XTickLabel = split(int2str(fliplr(0:5))); 
xlabel('\nu'); 
ax1v.XGrid = 'on'; 
ax1v.GridAlpha = 0.7; 
ax1v.YTick = []; 
% finalize plot
linkaxes([ax1,ax1v],'x');  
uistack(ax1, 'top'); % set(l, 'color', 'white', 'FontSize', legend_font_size); 
xlim([xr(1) xr(2)]); 
POS = ax1_2w.Position; 
ax1v_2w = axes('Position',POS);
ax1v_2w.XAxisLocation = 'top'; 
ax1v_2w.Color = 'none'; 
ax1v_2w.XTick = [fliplr(12*1240/810-IP) fliplr(14*1240/810-IP) fliplr(16*1240/810-IP)]; 
ax1v_2w.XTickLabel = []; 
% xlabel('\nu'); 
ax1v_2w.XGrid = 'on'; 
ax1v_2w.GridAlpha = 0.7; 
ax1v_2w.YTick = []; 
% finalize plot
linkaxes([ax1_2w,ax1v_2w],'x');  
uistack(ax1_2w, 'top'); % set(l, 'color', 'white', 'FontSize', legend_font_size); 
xlim([xr(1) xr(2)]); 

% %%%%%% second panel %%%%%%
% make color map
numc = 50; 
cm = ones([numc 3]); 
cm(1:(numc/2),1:2) = repmat((1:1:(numc/2))' ./ (numc/2), 1, 2); 
cm((numc/2+1):end,2:3) = flipud(cm(1:(numc/2),1:2)); 

boxwidth = 0.8/3; 

pos2_0V = [0.1 0+hoffset boxwidth 0.2];
ax2_0V = subplot('Position',pos2_0V); hold on; 
ESpectra_norm = ESpectra_0V ./ repmat(sum(ESpectra_0V,1), size(ESpectra_0V,1), 1, 1);
XUV_norm = XUV_only_0V ./ sum(XUV_only_0V);
% tmp = sum(abs(ESpectra_norm - repmat(XUV_norm, 1, 223, 37)),3);
tmp = ESpectra_norm; 
imagesc(E, stageTimes*1e15, ((mean(sum(abs(tmp),3),2)-sum(abs(tmp),3)) ./ mean(sum(abs(tmp),3),2))');
colormap(cm); 
% cb = colorbar; 
% cb.Label.String = 'XUV/IR - XUV only'; 
% cb.Ticks = []; 
% cb.Position = 'eastoutside'; 
caxis([-0.25 0.25]);  
xlim([xr(1) xr(1)+(xr(2)-xr(1))/3]); 
ylim([-4 4]); 
box 'on'; 
ax2_0V.FontSize = ax_font_size; 
ax2_0V.FontWeight = 'bold'; 
ax2_0V.LineWidth = 1; 
ax2_0V.Color = 'none'; 
ax2_0V.XTick = 2:1:10;
ax2_0V.XTickLabel = split(int2str(2:1:10)); 
ax2_0V.XMinorTick = 'on'; 
ax2_0V.XAxisLocation = 'bottom'; 
ax2_0V.YTick = [-3 0 3]; 
ax2_0V.YTickLabel = split(int2str([-3 0 3])); 
ax2_0V.YMinorTick = 'on'; 
% xlabel('electron kinetic energy (eV)'); 
ylabel('XUV/IR delay (fs)', 'FontSize', 10); 
text(2.3, -1.5, '0V', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 
plot(ones(20)*3.81, (-10:9)*1e15, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');  

pos2_3V = [0.1+boxwidth 0+hoffset boxwidth 0.2];
ax2_3V = subplot('Position',pos2_3V); hold on; 
ESpectra_norm = ESpectra_3V ./ repmat(sum(ESpectra_3V,1), size(ESpectra_3V,1), 1, 1);
XUV_norm = XUV_only_3V ./ sum(XUV_only_3V);
% tmp = sum(abs(ESpectra_norm - repmat(XUV_norm, 1, 223, 61)),3);
tmp = ESpectra_norm; 
imagesc(E, stageTimes*1e15, ((mean(sum(abs(tmp),3),2)-sum(abs(tmp),3)) ./ mean(sum(abs(tmp),3),2))');
colormap(cm); 
% cb = colorbar; 
% cb.Label.String = 'XUV/IR - XUV only'; 
% cb.Ticks = []; 
% cb.Position = 'eastoutside'; 
caxis([-0.25 0.25]);  
xlim([xr(1) xr(1)+(xr(2)-xr(1))/3]+(xr(2)-xr(1))/3); 
ylim([-4 4]); 
box 'on'; 
ax2_3V.FontSize = ax_font_size; 
ax2_3V.FontWeight = 'bold'; 
ax2_3V.LineWidth = 1; 
ax2_3V.Color = 'none'; 
ax2_3V.XTick = 2:1:10;
ax2_3V.YTick = []; 
ax2_3V.XTickLabel = split(int2str(2:1:10)); 
ax2_3V.XMinorTick = 'on'; 
ax2_3V.XAxisLocation = 'bottom'; 
% ax2.YTick = 0:50:250; 
% ax2.YTickLabel = split(int2str(0:50:250)); 
% ax2_3V.YMinorTick = 'on'; 
xlabel('electron kinetic energy (eV)'); 
% ylabel('XUV/IR delay (fs)'); 
text(5.3, -1.5, '3V', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 
plot(ones(20)*6.65, (-10:9)*1e15, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');  

pos2_5V = [0.1+2*boxwidth 0+hoffset boxwidth 0.2];
ax2_5V = subplot('Position',pos2_5V); hold on; 
ESpectra_norm = ESpectra_5V ./ repmat(sum(ESpectra_5V,1), size(ESpectra_5V,1), 1, 1);
XUV_norm = XUV_only_5V ./ sum(XUV_only_5V);
% tmp = sum(abs(ESpectra_norm - repmat(XUV_norm, 1, 223, 37)),3);
tmp = ESpectra_norm; 
imagesc(E+1.5+0.2286, stageTimes*1e15, ((mean(sum(abs(tmp),3),2)-sum(abs(tmp),3)) ./ mean(sum(abs(tmp),3),2))');
% imagesc(E+1.5+0.2286, stageTimes*1e15, tmp);
colormap(cm); 
caxis([-0.25 0.25]);  
xlim([xr(1) xr(1)+(xr(2)-xr(1))/3]+2*(xr(2)-xr(1))/3); 
ylim([-4 4]); 
box 'on'; 
ax2_5V.FontSize = ax_font_size; 
ax2_5V.FontWeight = 'bold'; 
ax2_5V.LineWidth = 1; 
ax2_5V.Color = 'none'; 
ax2_5V.XTick = 2:1:10;
ax2_5V.YTick = []; 
ax2_5V.XTickLabel = split(int2str(2:1:10)); 
ax2_5V.XMinorTick = 'on'; 
ax2_5V.XAxisLocation = 'bottom'; 
% ax2.YTick = 0:50:250; 
% ax2.YTickLabel = split(int2str(0:50:250)); 
% ax2_3V.YMinorTick = 'on'; 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('XUV/IR delay (fs)'); 
text(8.1, -1.5, '5.5V', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center'); 
cb = colorbar; 
% cb.Position = [0.065     0.58    0.02    0.12]; # upper left position
cb.Position = [0.92 0.31 0.02 0.15]; 
cb.Ticks = [-0.25 0.25]; 
cb.TickLabels = ["-25" "+25"];
title(cb, '%'); 
% goodplot(14); 
ax1v.FontSize = axv_font_size;  

% save
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [8.5 11]);
set(gcf,'PaperPosition',[0.5 0.5 5 7]);
set(gcf,'PaperPositionMode','Manual');

set(gcf, 'units', 'inch', 'position', plotpos);
set(gcf,'PaperUnits','inches','PaperPosition',plotpos)
saveas(gcf,'/Users/annawang/Box/writing/H2/paper/figures/H2_raw_summary.eps', 'epsc')

%% plot a 2w example
% text and color settings
text_size = 20; 
line_weight = 1.5; 
abs_color = [0 0 0.8]; 
phi_color = [0 0.7 1]; 
lgnd_pad = 1.4; 
lgnd_pos = 'best'; 
    
region = [1.55 2.8]; % sideband 12
tolerance = abs(Ebins(2)-Ebins(1)); 
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

figure; hold on; 
ax1 = gca; 
ax1_pos = ax1.Position; 
ax1.XColor = 'k'; 
ax1.XLabel.String = 'Photoelectron Energy (eV)'; 
ax1.YTickLabel = []; 
ax1.YLabel.String = 'Amplitude';
ax1.FontSize = text_size; 
ax1.YColor = abs_color; 
ax1.LineWidth = line_weight*0.5; 
hold(ax1, 'on'); 
s1 = plot(ax1, Ebins(start:stop), abs(twoOmega_0V(start:stop)), '-'); 
s1.MarkerEdgeColor = abs_color; 
s1.Color = abs_color; 
s1.LineWidth = line_weight; 
s1.DisplayName = 'amplitude data'; 

ax2_0V = axes('Position', ax1_pos, ...
        'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
        'Color', 'none', 'LineWidth', line_weight*0.5); 
ax2_0V.XColor = 'none'; 
ax2_0V.YColor = phi_color; 
ax2_0V.YLabel.String = 'phase (radians)'; 
ax2_0V.FontSize = text_size; 
hold(ax2_0V, 'on'); 
s2 = plot(Ebins(start:stop), unwrap(angle(twoOmega_0V(start:stop)))+pi, '-'); 
s2.MarkerEdgeColor = phi_color; 
s2.Color = phi_color; 
s2.LineWidth = line_weight; 
s2.DisplayName = 'phase data'; 
ax2_0V.YLabel.String = 'Phase (radians)';
ax2_0V.YLim = [0 2*pi]; 
ax2_0V.YTick = (0:(pi/4):2*pi); 
ax2_0V.YTickLabel = ["-\pi", "-3\pi/4", "-\pi/2", "-\pi/4", ...
                    "0", "\pi/4", "\pi/2", "3\pi/4", "\pi"];  

linkaxes([ax1,ax2_0V],'x'); 
ax2_0V.YLim = [0 2*pi]; 
% ax1.XLim = [xout(1), xout(end)]; 
%     ax1.YLim = [0, max(yout_abs)*lgnd_pad]; 
%     lgnd = legend([s1, l1, s2, l2]) %, 'Location', lgnd_pos); 
%     
% legend; 

% goodplot(22); 

%% referee requested plots
%% full XUV only, XUV+IR, and 2w spectrum for each data set

subplot_width = 1.5; 
subplot_height = 0.2; 
hoffset = 0.05; 
font_size = 12; 

% f = figure('Units', 'inches', 'Position', [1 1 1.2*subplot_width 6*subplot_height]); 
f = figure; 
hold on; 
set(gcf, 'units', 'inch', 'position', plotpos);

pos1 = [0.1 0.3+hoffset 0.8 0.3];
subplot('Position',pos1); hold on; 
plot(Ebins, movmean(XUV_only_0V,5)./sum(XUV_only_0V), '-', 'Color', [0.4 0.4 0.4]*1.2, ...
    'DisplayName', 'XUV only', ...
    'LineWidth', 2); 
plot(Ebins, sum(sum(ESpectra_0V,2),3)./sum(ESpectra_0V(:)), 'k-', ...
    'DisplayName', 'IR + XUV (delay averaged)', ...
    'LineWidth', 2); 
set(gca,'XTickLabel',[], 'XTick', 0:1:20);
set(gca,'YTickLabel',[], 'YTick', []);
legend('Color', [1 1 1], 'EdgeColor', 'none'); 
% goodplot(font_size); 
ax1 = gca; 
ntrue = 2:2:20;
%Harmonic Labels
for ii=1:numel(ntrue)
    nlbl{ii} = sprintf('%2.0f',ntrue(ii));
end
% axh = axes('Position', [ax1.Position(1) ax1.Position(2) ax1.Position(3)*1.5 ax1.Position(4)*1.5]); 
% linkaxes([ax1, axh])
% AddHarmonicAxis(ax1, IP, IP_label, 810, 1); 
for ii=1:numel(IP)   
    for jj=12:2:18
        x = jj*1240/810 - IP(ii); 
        plot([x, x], [0 0.01-ii*0.0005], 'k', 'HandleVisibility', 'off', 'LineWidth', 0.2); 
        text(x-0.1, 0.01-(ii-1)*0.0005, IP_label(ii))
    end
end
xlim([0 15]); 
ylim([0 0.011]); 
ylabel('Amplitude (arb.)')
box on; 

% pos2 = [0.1 0.3+hoffset 0.8 0.2];
% subplot('Position',pos2); hold on; 
% plot(E, sum(sum(E_SpectraArray,2),3), 'k-', ...
%     'DisplayName', '\tau averaged RABBITT'); 
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% % legend; 
% goodplot(font_size); 
% ax2 = gca; 
% AddHarmonicAxis(ax2, IP, IP_label, 810, 0); 

pos3 = [0.1 0.1+hoffset 0.8 0.2];
subplot('Position',pos3); hold on; 
plot(Ebins, abs(twoOmega_0V), 'b-', ...
    'DisplayName', '2\omega RABBITT amplitude', 'LineWidth', 2); 
set(gca, 'XTick', 0:1:20)
set(gca,'YTickLabel',[], 'YTick', []);
xlabel('electron kinetic energy (eV)'); 
legend('Color', [1 1 1], 'EdgeColor', 'none'); 
% goodplot(font_size); 
ax3 = gca; 
% AddHarmonicAxis(ax3, IP, IP_label, 810, 0); 
for ii=1:numel(IP)   
    for jj=12:2:18
        x = jj*1240/810 - IP(ii); 
        plot([x, x], [0 2e4], 'k', 'HandleVisibility', 'off', 'LineWidth', 0.2); 
    end
    
%     POS(4) = POS(4)-del;  
end
xlim([0 15]); 
ylim([0 1.5e4]); 
box on; 

set(gcf,'color','w');

%% save
% make large markers and lines
% set(findall(gcf, 'Type', 'Line'), 'MarkerSize', 4, 'LineWidth', 2); 
% set(findall(gcf, 'Type', 'ErrorBar'), 'MarkerSize', 4, 'LineWidth', 1); 
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [8.5 11]);
set(gcf,'PaperPosition',[0.5 0.5 7 1]);
set(gcf,'PaperPositionMode','Manual');

% save
set(gcf, 'units', 'inch', 'position', plotpos);
set(gcf,'PaperUnits','inches','PaperPosition',plotpos)
saveas(gcf,'/Users/annawang/Box/writing/H2/paper/figures/ref_0Vtimeaverage.eps', 'epsc')
% saveas(gcf,'/Users/annawang/Box/writing/H2/paper/figures/ref_fullspectrum_0V.fig')

%%

figure; hold on; 
plotpos = [0, 0, 10, 2];
% plotpos = [0, 0, 20, 14]/2.75; 
set(gcf, 'units', 'inch', 'position', plotpos);

plot(-(Ebins+IP(1)), fliplr(sum(sum(ESpectra_0V,2),3)), 'k-'); 
xlabel('Photon Energy (eV)'); 
% ylabel('Photoelectron Yield'); 
set(gca, 'Ytick', []); 
xlim([-35 -15]); 
set(gca, 'XTick', -35:5:-15, 'XTickLabel', split(int2str(fliplr(15:5:35)))); 
goodplot(20); 
