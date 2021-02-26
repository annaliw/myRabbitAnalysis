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
goodplot(24)

% integrate
E_AU = 27.2114; 
T_AU = 1/24.189; 
w_tmp = E ./ E_AU; 
t_interp = fun(p,E) .* T_AU; 
true_phase = fliplr(cumtrapz(fliplr(w_tmp), fliplr(t_interp))); 

figure; hold on; 
plot(E, true_phase); 
hw = 1240/810; 

RABBITT_phase = fliplr(cumtrapz(fliplr((E+hw)/E_AU), fliplr(fun(p,E+hw)*T_AU))) ...
    - fliplr(cumtrapz(fliplr((E-hw)/E_AU), fliplr(fun(p,E-hw)*T_AU))); 
plot(E, RABBITT_phase); 

xlabel('photoelectron energy (eV)'); 
ylabel('phase'); 
legend('"true" phase', 'RABBITT phase'); 
goodplot(24)

%% compare above to data

signal = squeeze(sum(twoOmega_signal,2)); 
SB=12; 
rV = 0;

if rV == 0
    Ar_SB12_slope = Ar_0V_SB12_slope; 
    Ar_SB14_slope = Ar_0V_SB14_slope; 
    Ar_SB16_slope = Ar_0V_SB16_slope; 
    Ar_SB18_slope = Ar_0V_SB18_slope; 
    
    Ar_SB12_phase = Ar_0V_SB12_phase; 
    Ar_SB14_phase = Ar_0V_SB14_phase; 
    Ar_SB16_phase = Ar_0V_SB16_phase; 
    Ar_SB18_phase = Ar_0V_SB18_phase; 
elseif rV == 3
    Ar_SB12_slope = Ar_3V_SB12_slope; 
    Ar_SB14_slope = Ar_3V_SB14_slope; 
    Ar_SB16_slope = Ar_3V_SB16_slope; 
    Ar_SB18_slope = Ar_3V_SB18_slope; 
    
    Ar_SB12_phase = Ar_3V_SB12_phase; 
    Ar_SB14_phase = Ar_3V_SB14_phase; 
    Ar_SB16_phase = Ar_3V_SB16_phase; 
    Ar_SB18_phase = Ar_3V_SB18_phase; 
end
    

if SB==12
    region = [1.6 2.8]; % sideband 12
    ar_slope = Ar_SB12_slope(1); 
    ar_phase = Ar_SB12_phase(1); 
    H2_phase = H2_0V_SB12_phase(1); 
elseif SB==14
    region = [4.7 6.15]; % sideband 14
    ar_slope = Ar_SB14_slope(1); 
elseif SB==16
    region = [7.7548 9.18]; % sideband 16
    % region = [7.72 9.18]; % sideband 16
    ar_slope = Ar_SB16_slope(1); 
elseif SB==18
    region = [10.8258 12.2]; % sideband 18
    ar_slope = Ar_SB18_slope(1); 
else
    region = [1.6 2.8]; 
    SB = 12; 
end

tolerance = 0.05; 
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 
xdata = E(start:stop); 
ydata = signal(start:stop); 

figure; hold on; 
yyaxis left; ylabel('amplitude'); 
p = plot(xdata, abs(ydata), 'k'); 
p.Color(4) = 0.2; 

% look inside specific vibrational states
for ii=2:numel(IP)
    v_energy = SB*1240/810 - IP(ii); 
    v_region = v_energy + [-0.14 0.14]; 
    v_start = find(abs(E-v_region(1))<tolerance, 1, 'last'); 
    v_stop = find(abs(E-v_region(2))<tolerance, 1, 'first'); 
    xtheory = E(v_start:v_stop); 
%     ytheory = RABBITT_phase(v_start:v_stop); 
    ydata = signal(v_start:v_stop); 

    ar_background = ar_slope * (xtheory - mean(xtheory)); 
%     ar_background = 0; 
    
    yyaxis right; ylabel('phase'); 
    tmp = angle(ydata)' - ar_background; 
    plot(xtheory, angle(ydata), 'c*-'); 
    plot(xtheory, tmp, 'bx-'); 
%     plot(xtheory, ytheory - mean(ytheory) + mean(tmp), 'k-'); 

end

xlabel('photoelectron energy (eV)'); 
title(strcat('SB ',num2str(SB))); 
goodplot(24)

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
goodplot(22)


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
SB=12; 

if SB==12
    region = [1.6 2.8]; % sideband 12
    ar_slope = Ar_SB12_slope(1); 
elseif SB==14
    region = [4.7 6.15]; % sideband 14
    ar_slope = Ar_SB14_slope(1); 
elseif SB==16
    region = [7.7548 9.18]; % sideband 16
    % region = [7.72 9.18]; % sideband 16
    ar_slope = Ar_SB16_slope(1); 
elseif SB==18
    region = [10.8258 12.2]; % sideband 18
    ar_slope = Ar_SB18_slope(1); 
else
    region = [1.6 2.8]; 
    SB = 12; 
end

tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first');  
xdata = E(start:stop); 
ydata = signal(start:stop); 

figure; hold on; 
yyaxis left; ylabel('amplitude'); 
plot(xdata, abs(ydata), 'k'); 

% look inside specific vibrational states
for ii=2:numel(IP)
    v_energy = SB*hw - IP(ii); 
    v_region = v_energy + [-0.14 0.14]; 
    v_start = find(abs(E-v_region(1))<tolerance, 1, 'last'); 
    v_stop = find(abs(E-v_region(2))<tolerance, 1, 'first'); 
    xtheory = E(v_start:v_stop); 
    ytheory_r140 = RABBITT_phase_r140(v_start:v_stop); 
    ytheory_r145 = RABBITT_phase_r145(v_start:v_stop); 
    ytheory_r150 = RABBITT_phase_r150(v_start:v_stop); 
    ydata = signal(v_start:v_stop); 

    ar_background = ar_slope*xtheory; 
%     ar_background = 0; 
    
    yyaxis right; ylabel('phase'); 
    tmp = angle(ydata)' - ar_background; 
    plot(xtheory, tmp, 'cx-'); 
    plot(xtheory, ytheory_r140 - mean(ytheory_r140) + mean(tmp), 'r.-'); 
    plot(xtheory, ytheory_r145 - mean(ytheory_r145) + mean(tmp), 'go-'); 
    plot(xtheory, ytheory_r150 - mean(ytheory_r150) + mean(tmp), 'b--'); 

end

xlabel('photoelectron energy (eV)'); 
title(strcat('SB ',num2str(SB))); 
goodplot()

%% instead of integrating for rainbow-RABBITT, plot different wavelengths: 

signal = squeeze(sum(twoOmega_signal,2)); 
SB=12; 

if SB==12
    region = [1.6 2.8]; % sideband 12
    ar_slope = Ar_SB12_slope(1); 
    ind = 1; 
elseif SB==14
    region = [4.7 6.15]; % sideband 14
    ar_slope = Ar_SB14_slope(1); 
    ind = 2; 
elseif SB==16
    region = [7.7548 9.18]; % sideband 16
    % region = [7.72 9.18]; % sideband 16
    ar_slope = Ar_SB16_slope(1); 
    ind = 3; 
elseif SB==18
    region = [10.8258 12.2]; % sideband 18
    ar_slope = Ar_SB18_slope(1); 
    ind = 4; 
else
    region = [1.6 2.8]; 
    SB = 12; 
    ind = 1; 
end

w = 1239.8 / 810 / E_AU; %Photon Energy
T_L = 2*pi/w / T_AU; % laser period in fs
% look inside v=1
v_energy = SB*hw - IP(2); 
v_region = v_energy + [-0.14 0.14]; 
tolerance = 0.05; 
v_start = find(abs(E-v_region(1))<tolerance, 1, 'last'); 
v_stop = find(abs(E-v_region(2))<tolerance, 1, 'first'); 
xtheory = [H2TDSE_810_140.x_Ee(ind), H2TDSE_785_140.Ee(ind), H2TDSE_760_140.Ee(ind)]; 
ytheory = [H2TDSE_810_140.t(ind), H2TDSE_785_140.t(ind), H2TDSE_760_140.t(ind)]./(T_L*1000/2/(2*pi)); 
ydata = signal(v_start:v_stop); 
xdata = E(v_start:v_stop); 

ar_background = ar_slope*xdata; 
%     ar_background = 0; 

figure; hold on; 
yyaxis left; ylabel('amplitude'); 
plot(xdata, abs(ydata), 'b'); 
yyaxis right; ylabel('phase'); 
tmp = angle(ydata)' - ar_background; 
plot(xdata, tmp, 'cx-'); 
plot(SB*1240./[810 785 760] - IP(3), ytheory - mean(ytheory) + mean(tmp), 'ko-'); 

xlabel('photoelectron energy (eV)'); 
title(strcat('SB ',num2str(SB))); 
goodplot()


%% R-RABBITT data analysis

% run fiterror_bootstrap for desired fit region. 
% new region, full sideband data: 
region = [1.2 3.5]; SB=12; 
% region = [4.4 6.5]; SB=14; 
if SB==12
    Ar_slope = Ar_SB12_slope(1); 
elseif SB==14
    Ar_slope = Ar_SB14_slope(1); 
else
    Ar_slope = 0; 
end


[~,start] = min(abs(E-region(1))); 
[~, stop] = min(abs(E-region(2))); 
xdata = E(start:stop); 
ydata = twoOmega_signal(start:stop).'; 

% loop over each peak
amp_array = zeros([numel(xdata), size(paramout_gauss,1)]); 
amp_array_shift = amp_array; 
ydata_shift = amp_array; 
FT_array_nophase = amp_array; 
FT_array_phase = FT_array_nophase; 
p_nophase = zeros([1 size(paramout_gauss,1)]); 
p_phase = p_nophase; 

for ii=1:size(paramout_gauss,1)
    p = paramout_gauss(ii,:); 
%     amp_array(:,ii) = 1/p(3)/2/pi .* exp(-(xdata-p(2)).^2/(2*p(3).^2)); 
    amp_array(:,ii) = p(1) .* exp(-(xdata-p(2)).^2/(2*p(3).^2)); 
    [~, c_ind] = min(abs(E-p(2))); 
    xshift = circshift(xdata, -(c_ind-ceil(numel(xdata)/2))); 
    xcenter = xdata-xdata(floor(numel(xdata)/2)); 
    if SB==12
        amp_array_shift(:,ii) = fftshift(circshift(amp_array(:,ii), -(c_ind-ceil(numel(xdata)/2)))); 
        ydata_shift(:,ii) = fftshift(circshift(ydata, -(c_ind-ceil(numel(xdata)/2)))); 
    elseif SB==14
        amp_array_shift(:,ii) = circshift(amp_array(:,ii), -(c_ind-ceil(numel(xdata)/2))); 
        ydata_shift(:,ii) = circshift(ydata, -(c_ind-ceil(numel(xdata)/2))); 
    end
    FT_array_nophase(:,ii) = fftshift(fft(ifftshift(amp_array_shift(:,ii) .* exp(1j*Ar_slope*xcenter).'))); 
    FT_array_phase(:,ii) = fftshift(fft(ifftshift(amp_array_shift(:,ii) .* exp(1j*angle(ydata_shift(:,ii)))))); 
    
    % measure change in peak width due to phase variation
    [~,tmp1] = min(abs(amp_array_shift(:,ii) - 0.8*max(amp_array_shift(:,ii)))); 
    [~, tmp2]  = min(abs(flipud(amp_array_shift(:,ii)) - 0.8*max(amp_array_shift(:,ii)))); 
    if tmp1>tmp2
        tmp_start = tmp2; 
        tmp_stop = tmp1; 
    elseif tmp1<tmp2
        tmp_start = tmp1; 
        tmp_stop = tmp2; 
    end
    
    fun = @(p, xdata) 1./p/2/pi .* exp(-(xdata-xdata(floor(numel(xdata)/2))).^2/(2*p.^2)); 
    p0 = paramout_gauss(ii,3); 
    tmp = fun(p0, xdata); 
    p_nophase(ii) = lsqcurvefit(fun, p0, xdata, abs(FT_array_nophase(:,ii))');
    p_phase(ii) = lsqcurvefit(fun, p0, xdata, abs(FT_array_phase(:,ii))');
 
end

wdata = xdata ./ 4.135; % h in eV/fs
dw = wdata(2)-wdata(1); 
dt = 1/(wdata(end)-wdata(1)); 
ft_xdata = (-0.5/dw):dt:(0.5/dw); 
if numel(ft_xdata)>numel(wdata)
    t_tmp = (-0.5/dw):dt:(0.5/dw - dt);
elseif numel(ft_xdata)<numel(wdata)
    t_tmp = (-0.5/dw):dt:(0.5/dw + dt);
end

figure; hold on; plot(xdata, amp_array); 
yyaxis right; plot(xdata, angle(ydata)); 
    xlabel('electron kinetic energy (eV)'); 
    goodplot(); legend('0', '1', '2', '3', '4', '5', '6');
    title(strcat('SB', num2str(SB))); 
figure; hold on; plot(ft_xdata, abs(FT_array_nophase)); 
    xlabel('time (fs)'); 
    goodplot(); legend('0', '1', '2', '3', '4', '5', '6');
    title(strcat('FT SB', num2str(SB), ', without phase')); 
figure; hold on; plot(ft_xdata, abs(FT_array_phase)); 
    xlabel('time (fs)'); 
    goodplot(); legend('0', '1', '2', '3', '4', '5', '6');
    title(strcat('FT SB', num2str(SB), ', with phase')); 
% figure; hold on; plot(xshift, amp_array_shift);
%     xlabel('electron kinetic energy (eV)'); 
%     goodplot(); legend('0', '1', '2', '3', '4', '5', '6');
%     title(strcat('SB', num2str(SB), ', shifted')); 
figure; hold on; plot((1:1:6)-1, p_nophase, 'o-'); 
                 plot((1:1:6)-1, p_phase, 's-'); 
    xlabel('v-state'); ylabel('peak width (eV)'); 
    goodplot(); legend('no phase FT', 'phase FT'); 
    title(strcat('SB', num2str(SB), ', peak broadening')); 


%% plotting in a grid

figure; hold on; 
subplot(3, 2, 1); hold on; 
plot(-100, 100, 'k-'); 
plot(-101, 101, 'k--'); 
plot(-102, 102, 'bo'); 
plot(-103, 103, 'ro'); 
legend('FT with phase', 'FT without phase', 'SB12', 'SB14'); 
xlim([-40 40]); 
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
goodplot(); 
for ii=2:size(paramout_gauss,1)
    subplot(3, 2, ii); hold on; 
    plot(ft_xdata_SB12, abs(FT_array_phase_SB12(:,ii)), 'b-'); 
    plot(ft_xdata_SB12, abs(FT_array_nophase_SB12(:,ii)), 'b--'); 
    plot(ft_xdata_SB14, abs(FT_array_phase_SB14(:,ii)), 'r-'); 
    plot(ft_xdata_SB14, abs(FT_array_nophase_SB14(:,ii)), 'r--'); 
    xlabel('time (fs)'); ylabel('amplitude'); 
    title(strcat('v=', num2str(ii-1))); 
    goodplot()
    xlim([-40 40]); 
end

figure; hold on; 
subplot(3, 2, 1); hold on; 
plot(-100, 100, 'k-'); 
plot(-102, 102, 'bo'); 
plot(-103, 103, 'ro'); 
legend('subtract FT amp with no phase from FT amp with phase', 'SB12', 'SB14'); 
xlim([-40 40]); 
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
goodplot(); 
for ii=2:size(paramout_gauss,1)
    subplot(3, 2, ii); hold on; 
    plot(ft_xdata_SB12, (abs(FT_array_phase_SB12(:,ii))-abs(FT_array_nophase_SB12(:,ii)))./abs(FT_array_phase_SB12(:,ii)), 'b-'); 
    plot(ft_xdata_SB14, (abs(FT_array_phase_SB14(:,ii))-abs(FT_array_nophase_SB14(:,ii)))./abs(FT_array_phase_SB14(:,ii)), 'r-'); 
    xlabel('time (fs)'); ylabel('difference'); 
    title(strcat('v=', num2str(ii-1))); 
    goodplot()
    xlim([-40 40]); 
end


%% Short Time FT

t_range = 50; % fs
t = 0:0.1:t_range; 
hbar = 4.135/2/pi; 
Delta_t = 15; % STFT window (gaussian envelope)

S_STFT_phase = zeros([numel(xdata), numel(t), size(paramout_gauss,1)]); 
S_STFT_nophase = S_STFT_phase; 
xcenter = xdata-xdata(floor(numel(xdata)/2)); 
for ii=1:size(paramout_gauss,1)
    for jj=1:numel(t)
        t_tmp = ft_xdata - t(jj); 
        alpha = exp(-t_tmp.^6/(2*Delta_t.^6)); 
        S_STFT_phase(:,jj,ii) = sum(FT_array_phase(:,ii).' .* alpha .* exp(1j* xcenter.' .*ft_xdata./hbar),2); 
        S_STFT_nophase(:,jj,ii) = sum(FT_array_nophase(:,ii).' .* alpha .* exp(1j* xcenter.' .*ft_xdata./hbar),2); 
%         S_LFT(:,jj,ii) = 
    end
%     figure; hold on; 
%     imagesc(xdata, t, abs(S_STFT_phase(:,:,ii).').^2 - abs(S_STFT_nophase(:,:,ii).').^2); colorbar; 
%     title(strcat('SB', num2str(SB), ' v=', num2str(ii))); 
%     ylim([-5 t_range]); 
end


%%
ii=6; 

figure; hold on; 
imagesc(xcenter, t, abs(S_STFT_phase(:,:,ii).').^2); 
% imagesc(xcenter, t, abs(S_STFT_phase(:,:,ii).').^2 - abs(S_STFT_nophase(:,:,ii).').^2); 
% colormap(flipud(hot)); colorbar; 
colormap(bone); colorbar; 
title(strcat('SB', num2str(SB), ' v=', num2str(ii-1))); 
ylim([-5 t_range]); 
xlim([-1 1])
xlabel('relative electron energy (eV)'); 
ylabel('time (fs)'); 
goodplot()

path = '/Users/annawang/Documents/presentations/DAMOP_2020/'; 
f = gcf;
export_fig(f, strcat(path, 'R-RABBITT_STFT_SB', num2str(SB), '_v', num2str(ii-1)), '-pdf', '-eps');

% figure; hold on; 
% imagesc(xdata, t, abs(S_STFT_phase(:,:,ii).' - S_STFT_nophase(:,:,ii).').^2); 
% colormap bone; colorbar; 
% title(strcat('SB', num2str(SB), ' v=', num2str(ii))); 
% ylim([-5 t_range]); 

%% 
ii=5; 

tpoints = 1:20:numel(t); 
cmap = zeros([numel(tpoints) 3]); 
cmap(:,3) = 1; 
cmap(:,2) = tpoints ./ max(tpoints); 

figure; hold on; 
for jj=tpoints
    plot(xcenter, abs(S_STFT_phase(:,jj,ii).').^2 - abs(S_STFT_nophase(:,jj,ii).').^2)
end
title(strcat('SB', num2str(SB), ' v=', num2str(ii-1))); 
xlim([-1 1])
xlabel('relative electron energy (eV)'); 
ylabel('time (fs)'); 
goodplot()
colororder(cmap); 
