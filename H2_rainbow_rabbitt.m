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

%% combine H2 TDSE into single data set

x_Ee = [H2TDSE_810_140.x_Ee; H2TDSE_810_145.x_Ee; H2TDSE_810_150.x_Ee]; 
n = [H2TDSE_810_140.n; H2TDSE_810_145.n; H2TDSE_810_150.n]; 
t = [H2TDSE_810_140.t; H2TDSE_810_145.t; H2TDSE_810_150.t]; 
err = [H2TDSE_810_140.err; H2TDSE_810_145.err; H2TDSE_810_150.err]; 

H2TDSE_810_all = table(x_Ee, n, t, err); 

%% fit line to combined H2 data set

% function given by Anatoli with free parameters A, B, N and independent
% variable E
fun = @(x,xdata) 1./xdata.^1 .* (x(1)*log(xdata) + x(2)); 
p0 = [100, -500]; 
[p,resnorm] = lsqcurvefit(fun, p0, H2TDSE_810_all.x_Ee, H2TDSE_810_all.t); 

figure; hold on; 
plot(H2TDSE_810_all.x_Ee, H2TDSE_810_all.t, 'bo'); 
plot(E, fun(p,E), 'r'); 
ylim([-300 10]); 
xlabel('photoelectron energy (eV)'); 
ylabel('delay (as)'); 
legend('data', strcat('p=[',num2str(p),'], resnorm=',num2str(resnorm))); 
goodplot()

% integrate
E_AU = 27.2114; 
T_AU = 1/24.189; 
xtmp = E ./ E_AU; 
ytmp = fun(p,E) .* T_AU; 
true_phase = fliplr(cumtrapz(fliplr(xtmp), fliplr(ytmp))); 

figure; hold on; 
plot(E, true_phase); 
hw = 1240/810; 

RABBITT_phase = fliplr(cumtrapz(fliplr((E+hw)/E_AU), fliplr(fun(p,E+hw)*T_AU))) ...
    - fliplr(cumtrapz(fliplr((E-hw)/E_AU), fliplr(fun(p,E-hw)*T_AU))); 
plot(E, RABBITT_phase); 

xlabel('photoelectron energy (eV)'); 
ylabel('phase'); 
legend('"true" phase', 'RABBITT phase'); 
goodplot()

%% compare above to data

signal = squeeze(sum(twoOmega_signal,2)); 

region = [1.6 2.8]; % sideband 12
% region = [4.7 6.15]; % sideband 14
% region = [7.7548 9.18]; % sideband 16
% region = [7.72 9.18]; % sideband 16
% region = [10.8258 12.2]; % sideband 18

% look inside specific vibrational states
v_energy = 14*hw - IP(3); 
% region = v_energy + [-0.2 0.2]; 

tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

% [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), ...
%     abs(signal(start:stop)), signal(start:stop), 1, 1); 
% % save as labeled variables
% paramout_original = paramout; 
% fval_original = fval 

xdata = E(start:stop); 
ydata = twoOmega_signal(start:stop); 
ytheory = RABBITT_phase(start:stop); 

% ar_background = Ar_SB14_phase(1) + Ar_SB14_slope(1)*xdata; 
ar_background = 0; 

figure; hold on; 
yyaxis left; ylabel('amplitude'); 
plot(xdata, abs(ydata)); 
yyaxis right; ylabel('phase'); 
tmp = angle(ydata)' - ar_background; 
plot(xdata, tmp - mean(tmp), 'x-'); 
plot(xdata, ytheory - mean(ytheory)); 
xlabel('photoelectron energy (eV)'); 
title('SB 12'); 
goodplot()

%% fit lines to different H2 r0 sets separately

fun = @(x,xdata) 1./xdata.^(x(3)) .* (x(1)*log(xdata) + x(2)); 
p0 = [100, -500, 1]; 
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
goodplot()


% integrate
E_AU = 27.2114; 
T_AU = 1/24.189; 

figure; hold on; 
hw = 1240/810; 

RABBITT_phase_r140 = fliplr(cumtrapz(fliplr((E+hw)/E_AU), fliplr(fun(p_r140,E+hw)*T_AU))) ...
    - fliplr(cumtrapz(fliplr((E-hw)/E_AU), fliplr(fun(p_r140,E-hw)*T_AU))); 
plot(E, RABBITT_phase_r140, 'r', 'DisplayName', 'r=140'); 
RABBITT_phase_r145 = fliplr(cumtrapz(fliplr((E+hw)/E_AU), fliplr(fun(p_r145,E+hw)*T_AU))) ...
    - fliplr(cumtrapz(fliplr((E-hw)/E_AU), fliplr(fun(p_r145,E-hw)*T_AU))); 
plot(E, RABBITT_phase_r145, 'g', 'DisplayName', 'r=145'); 
RABBITT_phase_r150 = fliplr(cumtrapz(fliplr((E+hw)/E_AU), fliplr(fun(p_r150,E+hw)*T_AU))) ...
    - fliplr(cumtrapz(fliplr((E-hw)/E_AU), fliplr(fun(p_r150,E-hw)*T_AU))); 
plot(E, RABBITT_phase_r150, 'b', 'DisplayName', 'r=150'); 

xlabel('photoelectron energy (eV)'); 
ylabel('phase'); 
legend; 
goodplot()

%% compare above to data

signal = squeeze(sum(twoOmega_signal,2)); 

% region = [1.6 2.8]; % sideband 12
region = [4.7 6.15]; % sideband 14
% region = [7.72 9.18]; % sideband 16
% region = [10.8258 12.2]; % sideband 18

% look inside specific vibrational states
v_energy = 14*hw - IP(3); 
% region = v_energy + [-0.2 0.2]; 

tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

% [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), ...
%     abs(signal(start:stop)), signal(start:stop), 1, 1); 
% % save as labeled variables
% paramout_original = paramout; 
% fval_original = fval 

xdata = E(start:stop); 
ydata = twoOmega_signal(start:stop); 
ytheory_r140 = RABBITT_phase_r140(start:stop); 
ytheory_r145 = RABBITT_phase_r145(start:stop); 
ytheory_r150 = RABBITT_phase_r150(start:stop); 

% ar_background = Ar_SB14_phase(1) + Ar_SB14_slope(1)*xdata; 
ar_background = 0; 

figure; hold on; 

yyaxis left; ylabel('amplitude'); 
plot(xdata, abs(ydata)); 

yyaxis right; ylabel('phase'); 
tmp = angle(ydata)' - ar_background; 
plot(xdata, tmp - mean(tmp), 'x-'); 
plot(xdata, ytheory_r140 - mean(ytheory_r140), 'r'); 
plot(xdata, ytheory_r145 - mean(ytheory_r145), 'g'); 
plot(xdata, ytheory_r150 - mean(ytheory_r150), 'b'); 
xlabel('photoelectron energy (eV)'); 
title('SB 14'); 
goodplot()











