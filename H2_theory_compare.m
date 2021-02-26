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
    
    