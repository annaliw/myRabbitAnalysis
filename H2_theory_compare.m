%% Import Nishi

tmp = csvread('/Users/annawang/Documents/data/H2plus_Nishi/H2plus_vphase_H11_Nishi.csv'); 
Nishi_H11_dvv = tmp(:,2); 

tmp = csvread('/Users/annawang/Documents/data/H2plus_Nishi/H2plus_vphase_H13_Nishi.csv'); 
Nishi_H13_dvv = tmp(:,2); 

tmp = csvread('/Users/annawang/Documents/data/H2plus_Nishi/H2plus_vphase_H15_Nishi.csv'); 
Nishi_H15_dvv = tmp(:,2); 

%% Import Anatoli's TDSE values for Argon, H2

% import from CSV and convert to table to match H2 calculations
tmp = csvread('/Users/annawang/Documents/plots/Calculations/H2_TDSE/ArTDSE_810.csv',0,0); 
x_Ee = tmp(:,1)-15.736; 
% t = 2*flipud(tmp(:,2)); 
t = tmp(:,2); 
n = 12:2:24; 
ArTDSE_810 = table(x_Ee, n', t); 

% import H2 as tables
H2TDSE_810_145 = readtable('/Users/annawang/Documents/plots/Calculations/H2_TDSE/H2-R145.dat');
H2TDSE_810_150 = readtable('/Users/annawang/Documents/plots/Calculations/H2_TDSE/H2-R150.dat');
H2TDSE_810_140 = readtable('/Users/annawang/Documents/plots/Calculations/H2_TDSE/H2-omega056.dat');
H2TDSE_785_140 = readtable('/Users/annawang/Documents/plots/Calculations/H2_TDSE/H2-omega058.dat');
H2TDSE_760_140 = readtable('/Users/annawang/Documents/plots/Calculations/H2_TDSE/H2-omega06.dat');

%% fit lines to different H2 r0 sets separately

fun = @(x,xdata) 1./xdata .* (x(1)*log(xdata) + x(2)); 
p0 = [100, -500]; 
[p_r140,resnorm_r140] = lsqcurvefit(fun, p0, H2TDSE_810_140.x_Ee, H2TDSE_810_140.t); 
[p_r145,resnorm_r145] = lsqcurvefit(fun, p0, H2TDSE_810_145.x_Ee, H2TDSE_810_145.t); 
[p_r150,resnorm_r150] = lsqcurvefit(fun, p0, H2TDSE_810_150.x_Ee, H2TDSE_810_150.t); 

figure; hold on; 

plot(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t, 'ro', 'DisplayName', 'r0=140'); 
plot(E, fun(p_r140,E), 'r', 'DisplayName', ...
    strcat('p=[',num2str(p_r140),'], resnorm=',num2str(resnorm_r140))); 

plot(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t, 'go', 'DisplayName', 'r0=145'); 
plot(E, fun(p_r145,E), 'g', 'DisplayName', ...
    strcat('p=[',num2str(p_r145),'], resnorm=',num2str(resnorm_r145))); 

plot(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t, 'bo', 'DisplayName', 'r0=150'); 
plot(E, fun(p_r150,E), 'b', 'DisplayName', ...
    strcat('p=[',num2str(p_r150),'], resnorm=',num2str(resnorm_r150))); 

ylim([-300 10]); 
xlabel('photoelectron energy (eV)'); 
ylabel('delay (as)'); 
legend; 
goodplot(22)

fit_140 = fun(p_r140,E); 
fit_145 = fun(p_r145,E); 
fit_150 = fun(p_r150,E);

%% interpolate between Anatoli's points
interp_145 = interp1(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t, H2TDSE_810_150.x_Ee); 
interp_140 = interp1(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t, H2TDSE_810_145.x_Ee); 

% figure; hold on; 
% plot(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t, 'bo-', 'DisplayName', '145 orignal'); 
% plot(H2TDSE_810_140.x_Ee, interp_145, 'bo--', 'DisplayName', '145 interp to 140'); 
% plot(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t, 'ro-', 'DisplayName', '150 orignal'); 
% plot(H2TDSE_810_145.x_Ee, interp_150, 'ro--', 'DisplayName', '150 interp to 145'); 

%% get delta v values from Anatoli's data and fits

% en_ind_140 = zeros([1 numel(H2TDSE_810_145.x_Ee)]); 
% en_ind_145 = zeros([1 numel(H2TDSE_810_140.x_Ee)]); 
% en_ind_150 = zeros([1 numel(H2TDSE_810_150.x_Ee)]); 
% for ii=1:numel(H2TDSE_810_145.x_Ee)
%     [~,ind] = min(abs(E - H2TDSE_810_140.x_Ee(ii))); 
%     en_ind_140(ii) = ind; 
%     
%     [~,ind] = min(abs(E - H2TDSE_810_145.x_Ee(ii))); 
%     en_ind_145(ii) = ind; 
%     
%     [~,ind] = min(abs(E - H2TDSE_810_150.x_Ee(ii))); 
%     en_ind_150(ii) = ind; 
% end
% 
% dvv_anatoli = [fit_150(en_ind_150) - fit_145(en_ind_150), fit_145(en_ind_150) - fit_140(en_ind_150); 
%                fit_150(en_ind_145) - fit_145(en_ind_145), fit_145(en_ind_145) - fit_140(en_ind_145); 
%                fit_150(en_ind_140) - fit_145(en_ind_140), fit_145(en_ind_140) - fit_140(en_ind_140)];
dvv_anatoli = [H2TDSE_810_150.t - interp_145, H2TDSE_810_145.t - interp_140]; 
dvv_anatoli = dvv_anatoli ./ (T_L*1000/2/(2*pi)); 

dvv_Nishi = -[Nishi_H13_dvv-Nishi_H11_dvv, Nishi_H15_dvv-Nishi_H13_dvv] * pi; 

figure; hold on; 

subplot(2,2,1); hold on; 
plot(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t ./ (T_L*1000/2/(2*pi)), 'ro', 'DisplayName', 'r0=140'); 
plot(E, interp1(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t, E, 'spline') ./ (T_L*1000/2/(2*pi)), 'r', 'HandleVisibility', 'off'); 
plot(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t ./ (T_L*1000/2/(2*pi)), 'go', 'DisplayName', 'r0=145'); 
plot(E, interp1(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t, E, 'spline') ./ (T_L*1000/2/(2*pi)), 'g', 'HandleVisibility', 'off'); 
plot(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t ./ (T_L*1000/2/(2*pi)), 'bo', 'DisplayName', 'r0=150'); 
plot(E, interp1(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t, E, 'spline') ./ (T_L*1000/2/(2*pi)), 'b', 'HandleVisibility', 'off'); 
xlabel('electron kinetic energy (eV)'); 
ylabel('two photon \Delta\phi (radians)');  
xlim([2 3]); 
title('SB 12'); 
legend; 
goodplot(18)

subplot(2,2,3); hold on; 
plot(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t ./ (T_L*1000/2/(2*pi)), 'ro', 'DisplayName', 'r0=140'); 
plot(E, interp1(H2TDSE_810_140.x_Ee, H2TDSE_810_140.t, E, 'spline') ./ (T_L*1000/2/(2*pi)), 'r', 'HandleVisibility', 'off'); 
plot(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t ./ (T_L*1000/2/(2*pi)), 'go', 'DisplayName', 'r0=145'); 
plot(E, interp1(H2TDSE_810_145.x_Ee, H2TDSE_810_145.t, E, 'spline') ./ (T_L*1000/2/(2*pi)), 'g', 'HandleVisibility', 'off'); 
plot(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t ./ (T_L*1000/2/(2*pi)), 'bo', 'DisplayName', 'r0=150'); 
plot(E, interp1(H2TDSE_810_150.x_Ee, H2TDSE_810_150.t, E, 'spline') ./ (T_L*1000/2/(2*pi)), 'b', 'HandleVisibility', 'off'); 
xlabel('electron kinetic energy (eV)'); 
ylabel('two photon \Delta\phi (radians)'); 
xlim([5 6.5]); 
title('SB 14'); 
goodplot(18)

subplot(2,2,2); hold on; 
plot([0 1], [dvv_anatoli(1,1) dvv_anatoli(1,2)], 'rd', 'DisplayName', 'TDSE')
plot((1:numel(Nishi_H13_dvv))-1, dvv_Nishi(:,1), 'ko', 'DisplayName', 'Nishi'); 
xlabel('v state'); 
ylabel('\Delta_{\nu,\nu+1} (radians)'); 
xlim([0 5]); 
title('SB 12'); 
legend; 
goodplot(18)

subplot(2,2,4); hold on;
plot([0 1], [dvv_anatoli(2,1) dvv_anatoli(2,2)], 'rd', 'DisplayName', 'TDSE'); 
plot((1:numel(Nishi_H13_dvv))-1, dvv_Nishi(:,2), 'ko', 'DisplayName', 'Nishi');
xlabel('v state'); 
ylabel('\Delta_{\nu,\nu+1} (radians)'); 
xlim([0 5]); 
title('SB 14'); 
goodplot(18)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRY VECTOR QUANTIFICATION OF THEORY DISCREPANCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xdata = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
phase_data = squeeze(reshape(SB_delay_data(2,:,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, :,:), [1 15]));
offset = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
    - ((mean(unwrap(H2_0V_SB18_phase(2:5)))+mean(unwrap(H2_3V_SB18_phase(2:5)))+mean(unwrap(H2_5V_SB18_phase(2:5)))+2*pi)/3*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 
plot_data = phase_data - reshape(repmat(XUV_delay(1:3), [1 5])', [1 15]) + offset; 
plot_error = phase_error;  
plot_data = plot_data ./ (T_L*1000/2/(2*pi)); 
plot_error = phase_error ./ (T_L*1000/2/(2*pi)); 

% SB12
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(1) H2TDSE_810_145.x_Ee(1) H2TDSE_810_150.x_Ee(1)]; 
yt = -[H2TDSE_810_140.t(1)    H2TDSE_810_145.t(1)    H2TDSE_810_150.t(1)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(1:5)) - exp(1j*yt_int(1:5)); 

figure; hold on; 

subplot(3,3,1); hold on; 
plot(xdata(1:5), abs(diff), 'bo'); 
xlabel('photoelectron energy (eV)'); 
ylabel('amplitude'); 
goodplot(18); 
subplot(3,3,2); hold on; 
plot(xdata(1:5), angle(diff) , 'kd'); 
xlabel('photoelectron energy (eV)'); 
ylabel('phase'); 
goodplot(18); 

% SB14
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(2) H2TDSE_810_145.x_Ee(2) H2TDSE_810_150.x_Ee(2)]; 
yt = -[H2TDSE_810_140.t(2)    H2TDSE_810_145.t(2)    H2TDSE_810_150.t(2)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(6:10)) - exp(1j*yt_int(6:10)); 

subplot(3,3,4); hold on; 
plot(xdata(6:10), [abs(diff(1)) -abs(diff(2:end))], 'bo'); 
xlabel('photoelectron energy (eV)'); 
ylabel('amplitude'); 
goodplot(18); 
subplot(3,3,5); hold on; 
plot(xdata(6:10), [angle(diff(1)) angle(diff(2:end))-pi], 'kd'); 
xlabel('photoelectron energy (eV)'); 
ylabel('phase'); 
goodplot(18); 

% SB16
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(3) H2TDSE_810_145.x_Ee(3) H2TDSE_810_150.x_Ee(3)]; 
yt = -[H2TDSE_810_140.t(3)    H2TDSE_810_145.t(3)    H2TDSE_810_150.t(3)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(11:15)) - exp(1j*yt_int(11:15)); 

subplot(3,3,7); hold on; 
plot(xdata(11:15), abs(diff), 'bo'); 
xlabel('photoelectron energy (eV)'); 
ylabel('amplitude'); 
goodplot(18); 
subplot(3,3,8); hold on; 
plot(xdata(11:15), angle(diff), 'kd'); 
xlabel('photoelectron energy (eV)'); 
ylabel('phase'); 
goodplot(18); 

% add the old plots for sanity check
ax1 = subplot(3,3,3); hold on; 
errorbar(xdata(1:5), plot_data(1:5)*(T_L*1000/2/(2*pi)), plot_error(1:5)*(T_L*1000/2/(2*pi)), 'ko', 'DisplayName', 'H2 measurement'); 
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
ylim([100 270]); 
window = 170; 

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
ylim([100 270]); 
uistack(ax1, 'top'); 

ax2 = subplot(3, 3, 6); hold on; 
errorbar(xdata(6:10), plot_data(6:10)*(T_L*1000/2/(2*pi)), plot_error(6:10)*(T_L*1000/2/(2*pi)), 'ko', 'DisplayName', 'H2 measurement'); 
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


ax3 = subplot(3, 3, 9); hold on; 
errorbar(xdata(11:15), plot_data(11:15)*(T_L*1000/2/(2*pi)), plot_error(11:15)*(T_L*1000/2/(2*pi)), 'ko', 'DisplayName', 'H2 measurement'); 
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

% figure; hold on; 
% for ii=1:numel(diff)
%     plot([0 abs(diff(ii))*cos(angle(diff(ii)))], [0 abs(diff(ii))*sin(angle(diff(ii)))]); 
% end

%% M dot delta

% SB12
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(1) H2TDSE_810_145.x_Ee(1) H2TDSE_810_150.x_Ee(1)]; 
yt = -[H2TDSE_810_140.t(1)    H2TDSE_810_145.t(1)    H2TDSE_810_150.t(1)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(1:5)) - exp(1j*yt_int(1:5)); 
figure; hold on; 
plot(xdata(1:5), abs(diff) .* cos(angle(diff) - yt_int(1:5) + pi/2), 'bo'); 

% SB14
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(2) H2TDSE_810_145.x_Ee(2) H2TDSE_810_150.x_Ee(2)]; 
yt = -[H2TDSE_810_140.t(2)    H2TDSE_810_145.t(2)    H2TDSE_810_150.t(2)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(6:10)) - exp(1j*yt_int(6:10)); 
plot(xdata(6:10), abs(diff) .* cos([angle(diff(1)) angle(diff(2:5))-pi] - yt_int(6:10) + pi/2), 'bo'); 

% SB16
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(3) H2TDSE_810_145.x_Ee(3) H2TDSE_810_150.x_Ee(3)]; 
yt = -[H2TDSE_810_140.t(3)    H2TDSE_810_145.t(3)    H2TDSE_810_150.t(3)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(11:15)) - exp(1j*yt_int(11:15)); 
plot(xdata(11:15), abs(diff) .* cos(angle(diff) - yt_int(11:15)), 'bo'); 

xlabel('photoelectron energy (eV)'); 
ylabel('Dot Product'); 
goodplot(18); 

%% plain subraction of exp and theory
%% subtract exp and theory plot

vx = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
xdata = reshape(mean_data(:,2,:), [1 15]); xdata(11:15) = xdata(11:15)+1.5; %+0.2286; 
phase_data = squeeze(reshape(SB_delay_data(2,:,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, :,:), [1 15]));

offset = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
    - ((mean(unwrap(H2_0V_SB18_phase(2:5)))+mean(unwrap(H2_3V_SB18_phase(2:5)))+mean(unwrap(H2_5V_SB18_phase(2:5)))+2*pi)/3*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 

plot_data = phase_data - reshape(repmat(XUV_delay(1:3), [1 5])', [1 15])+offset; plot_data = plot_data/(T_L*1000/4/pi); 
plot_error = phase_error; plot_error = plot_error/(T_L*1000/4/pi); 

diff = []; 

figure; hold on; 
plotpos = [0, 0, 6, 15]; 
set(gcf, 'units', 'inch', 'position', plotpos);
% set(gcf,'PaperUnits','inches','PaperPosition',plotpos/5); 

% SB12
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(1) H2TDSE_810_145.x_Ee(1) H2TDSE_810_150.x_Ee(1)]; 
yt = -[H2TDSE_810_140.t(1)    H2TDSE_810_145.t(1)    H2TDSE_810_150.t(1)] / (T_L*1000/4/pi);
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = [diff plot_data(1:5)-yt_int(1:5)]; 

ax1 = subplot(3, 1, 1); hold on; 
plot(xdata(1:5), diff, 'ko'); 
xlim([1.4 3])
ylim([-pi/8 0]); 
window = pi/16; 

text(1.4+1, -pi/16+window-pi/32, 'Sideband 12', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
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
ax1v = axes('Position',get(ax1,'Position'));
ax1v.XAxisLocation = 'top'; 
ax1v.Color = 'none'; 
ax1v.XTick = vx(1:5); 
ax1v.XTickLabel = split(int2str(fliplr(1:5))); 
xlabel('\nu'); 
ax1v.XGrid = 'on'; 
ax1v.GridAlpha = 0.7; 
ax1v.YTick = []; 
ax1v.FontSize = axv_font_size; 
% ax1v.FontWeight = 'bold'; 

linkaxes([ax1,ax1v],'x'); 
xlim([1.4 3])
uistack(ax1, 'top'); 


% SB14
ax2 = subplot(3, 1, 2); hold on; 
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(2) H2TDSE_810_145.x_Ee(2) H2TDSE_810_150.x_Ee(2)]; 
yt = -[H2TDSE_810_140.t(2)    H2TDSE_810_145.t(2)    H2TDSE_810_150.t(2)] / (T_L*1000/4/pi); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = [diff plot_data(6:10)-yt_int(6:10)]; 
plot(xdata(6:10), diff(6:10), 'ko'); 
ylabel('\theta^{\prime}-\theta (radians)');
% ylabel('delay (as)'); 
% legend; 
xlim([4.4 6.2]); 
ylim([0-window 0+window]); 

text(4.4+1, 0+window-pi/32, 'Sideband 14', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
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
uistack(ax2, 'top'); 

% SB16
ax3 = subplot(3, 1, 3); hold on; 
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(3) H2TDSE_810_145.x_Ee(3) H2TDSE_810_150.x_Ee(3)]; 
yt = -[H2TDSE_810_140.t(3)    H2TDSE_810_145.t(3)    H2TDSE_810_150.t(3)] / (T_L*1000/4/pi); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = [diff plot_data(11:15)-yt_int(11:15)]; 
plot(xdata(11:15), diff(11:15), 'ko'); 
xlabel('photoelectron energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([7.7 9.2]); 
ylim([0-window 0+window]); 

text(7.7+1, 0+window-pi/32, 'Sideband 16', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
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
ax3v.XTick = vx(11:15); 
ax3v.XTickLabel = split(int2str(fliplr(1:5))); 
ax3v.XGrid = 'on'; 
ax3v.GridAlpha = 0.7; 
ax3v.YTick = []; 
ax3v.FontSize = axv_font_size; 
% ax3.FontWeight = 'bold'; 

linkaxes([ax3,ax3v],'x'); 
xlim([7.5 9.2]); 
ylim([-15-window 15+window]); 
uistack(ax3, 'top'); 

% make large markers and lines
set(findall(gcf, 'Type', 'Line'), 'MarkerSize', 10, 'LineWidth', 2); 
set(findall(gcf, 'Type', 'ErrorBar'), 'MarkerSize', 10, 'LineWidth', 2); 
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [6 8]);
set(gcf,'PaperPosition',[0.5 0.5 7 7]);
set(gcf,'PaperPositionMode','Manual');

% save
set(gcf, 'units', 'inch', 'position', plotpos);
set(gcf,'PaperUnits','inches','PaperPosition',plotpos)
saveas(gcf,'/Users/annawang/Box/writing/H2/paper/figures/thetasub.png')

%% M dot d = M^2 sin^2(theta/2)

vx = reshape(repmat(12:2:16, [5 1])*1240/810 - repmat(fliplr(IP(2:end))', [1 3]), [1 15]); 
xdata = reshape(mean_data(:,2,:), [1 15]); xdata(11:15) = xdata(11:15)+1.5; %+0.2286; 
phase_data = squeeze(reshape(SB_delay_data(2,:,:), [1 15])); 
phase_error = squeeze(reshape(SB_delay_error(2, :,:), [1 15]));

offset = -mean([H2TDSE_810_140.t(4) H2TDSE_810_145.t(4) H2TDSE_810_150.t(4)]) ...
    - ((mean(unwrap(H2_0V_SB18_phase(2:5)))+mean(unwrap(H2_3V_SB18_phase(2:5)))+mean(unwrap(H2_5V_SB18_phase(2:5)))+2*pi)/3*(T_L*1000/2/(2*pi)) - XUV_delay(4)); 

plot_data = phase_data - reshape(repmat(XUV_delay(1:3), [1 5])', [1 15])+offset; plot_data = plot_data/(T_L*1000/4/pi); 
plot_error = phase_error; plot_error = plot_error/(T_L*1000/4/pi); 

diff = []; 

figure; hold on; 
% plotpos = [0, 0, 20, 6]/1.4; 
plotpos = [0, 0, 6, 15]; 
set(gcf, 'units', 'inch', 'position', plotpos);

% SB12
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(1) H2TDSE_810_145.x_Ee(1) H2TDSE_810_150.x_Ee(1)]; 
yt = -[H2TDSE_810_140.t(1)    H2TDSE_810_145.t(1)    H2TDSE_810_150.t(1)] / (T_L*1000/4/pi);
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = [diff plot_data(1:5)-yt_int(1:5)]; 

ax1 = subplot(3, 1, 1); hold on; 
plot(xdata(1:5), sin(abs(diff(1:5)/2)).^2, 'ko'); 
xlim([1.4 3])
ylim([0.005 0.02])
window = 50; 

text(1.4+1, 1.2, 'Sideband 12', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
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
ax1v.XTick = vx(1:5); 
ax1v.XTickLabel = split(int2str(fliplr(1:5))); 
xlabel('\nu'); 
ax1v.XGrid = 'on'; 
ax1v.GridAlpha = 0.7; 
ax1v.YTick = []; 
ax1v.FontSize = axv_font_size; 
% ax1v.FontWeight = 'bold'; 

linkaxes([ax1,ax1v],'x'); 
xlim([1.4 3])
% ylim([70 180]); 
uistack(ax1, 'top'); 


% SB14
ax2 = subplot(3, 1, 2); hold on; 
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(2) H2TDSE_810_145.x_Ee(2) H2TDSE_810_150.x_Ee(2)]; 
yt = -[H2TDSE_810_140.t(2)    H2TDSE_810_145.t(2)    H2TDSE_810_150.t(2)] / (T_L*1000/4/pi); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = [diff plot_data(6:10)-yt_int(6:10)]; 
plot(xdata(6:10), sin(abs(diff(6:10)/2)).^2, 'ko'); 
xlabel('photoelectron energy (eV)'); 
ylabel('M\cdot\delta');
% legend; 
xlim([4.4 6.2]); 
ylim([-0.001 0.002]) 

text(4.4+1, 1.2, 'Sideband 14', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
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
% ylim([120-window 120+window]); 
uistack(ax2, 'top'); 

% SB16
ax3 = subplot(3, 1, 3); hold on; 
% fit a line to theory points
xt = [H2TDSE_810_140.x_Ee(3) H2TDSE_810_145.x_Ee(3) H2TDSE_810_150.x_Ee(3)]; 
yt = -[H2TDSE_810_140.t(3)    H2TDSE_810_145.t(3)    H2TDSE_810_150.t(3)] / (T_L*1000/4/pi); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = [diff plot_data(11:15)-yt_int(11:15)]; 
plot(xdata(11:15), sin(abs(diff(11:15)/2)).^2, 'ko'); 
xlabel('photoelectron energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([7.5 9.2]); 
ylim([-0.001 0.004])

text(7.5+1, 1.2, 'Sideband 16', 'FontSize', label_font_size, 'FontWeight', 'bold', ...
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
ax3v.XTick = vx(11:15); 
ax3v.XTickLabel = split(int2str(fliplr(1:5))); 
% xlabel('\nu'); 
ax3v.XGrid = 'on'; 
ax3v.GridAlpha = 0.7; 
ax3v.YTick = []; 
ax3v.FontSize = axv_font_size; 
% ax3.FontWeight = 'bold'; 

linkaxes([ax3,ax3v],'x'); 
xlim([7.5 9.2]); 
% ylim([55-window 55+window]); 
uistack(ax3, 'top'); 

% make large markers and lines
set(findall(gcf, 'Type', 'Line'), 'MarkerSize', 10, 'LineWidth', 2); 
set(findall(gcf, 'Type', 'ErrorBar'), 'MarkerSize', 10, 'LineWidth', 2); 
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [6 8]);
set(gcf,'PaperPosition',[0.5 0.5 7 7]);
set(gcf,'PaperPositionMode','Manual');

% save
set(gcf, 'units', 'inch', 'position', plotpos);
set(gcf,'PaperUnits','inches','PaperPosition',plotpos)
saveas(gcf,'/Users/annawang/Box/writing/H2/paper/figures/Mdotdelta.png')

%%

xt = [H2TDSE_810_140.x_Ee(1) H2TDSE_810_145.x_Ee(1) H2TDSE_810_150.x_Ee(1)]; 
yt = -[H2TDSE_810_140.t(1)    H2TDSE_810_145.t(1)    H2TDSE_810_150.t(1)] ./ (T_L*1000/2/(2*pi)); 
yt_int = interp1(xt, yt, xdata, 'spline'); 
diff = exp(1j*plot_data(1:5)) - exp(1j*yt_int(1:5)); 

ii=5; 
figure; hold on; 
plot([0 cos(plot_data(ii))], [0 sin(plot_data(ii))], 'k-'); 
plot([0 cos(yt_int(ii))], [0 sin(yt_int(ii))], 'r-'); 



    
    