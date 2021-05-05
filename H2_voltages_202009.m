%% load data

folderName_Ar = '/Users/annawang/Documents/data/2020_09_01-11Scan/'; % Argon 0V
% folderName_Ar = '/Users/annawang/Documents/data/2020_09_01-12Scan/'; % Argon 3V
% folderName_Ar = '/Users/annawang/Documents/data/2020_09_01-14Scan/'; % Argon 8.5V

wavelength = 810; 
global IP; IP = [15.7596];
global IP_label; IP_label = ["Ar 2P 3/2"];

[HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderName_Ar, 1, wavelength);  
figure; plot(sum(sum(HistTot_array,2),3)); 

%% ARGON ENERGY CALIBRATION

t0 = 23; 
n = 16:1:21; 
calibEnergy = 1240/810 * n - IP; 
tof_peak = [2035 1356 1081 930 824 752]; 
calibType = 'Ar'; 

config.calibEnergy = calibEnergy; 
config.tofPeaks = tof_peak;   % redo with peak finding
config.IPcal = IP; 
config.Plot = 1; 

A = ECalibrate(t0, n, wavelength, calibType, config); % TO DO: hard set t0 into ECalibrate
filename = strcat(folderName_Ar, 'calibration/', calibType, '_calibration'); 
save(filename, 'A', 'calibEnergy', 'tof_peak', 'wavelength'); 

%% CHECK ARGON CALIBRATION

tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); 
E_vec = [0 20 900]; 

[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be used in loop
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
Ebins = E; % save for later in case E gets overwritten
norm = sum(E_SpectraArray,1); 

plotfun_rabbitspectrum(11:1:21, 810, E, sum(sum(E_SpectraArray,2),3), 'average');

%% H2 analysis

global IP; IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206] + 0.05;
global IP_label; IP_label = ["0", "1", "2", "3", "4", "5"]; % start 1

% load H2 file
% folderName_H2 = '/Users/annawang/Documents/data/2020_08_27-15Scan/'; v=0; % H2 0V
% folderName_H2 = '/Users/annawang/Documents/data/2020_08_28-15Scan/'; v=3; % H2 3V
folderName_H2 = '/Users/annawang/Documents/data/2020_08_30-15Scan/'; v=5; % H2 5V
% folderName_H2 = '/Users/annawang/Documents/data/2020_08_31-18Scan/'; v=8.5; % H2 8.5V

% load correct energy calibration
load(strcat(folderName_Ar, 'calibration/Ar_calibration.mat')); 

[HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderName_H2, 1, wavelength);  
% figure; plot(sum(sum(HistTot_array,2),3)); 

tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); 
E_vec = [0 20 900]; 

[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be used in loop
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
Ebins = E; % save for later in case E gets overwritten
norm = sum(E_SpectraArray,1); 

plotfun_rabbitspectrum(11:1:21, 810, E, sum(sum(E_SpectraArray,2),3), 'average');

if v==0
    ESpectra_0V = E_SpectraArray; 
elseif v==3
    ESpectra_3V = E_SpectraArray; 
elseif v==5
    ESpectra_5V = E_SpectraArray; 
elseif v==8.5
    ESpectra_8V = E_SpectraArray; 
end

%% BEGIN 2w ANALYSIS
% choose your data
data = ESpectra_5V; 
% data = E_SpectraArray; 

% filter 2w
tmp = fftshift(fft(ifftshift(data - repmat(mean(data,2),[1 numel(stageTimes)]),2),[],2),2); 
twoOmega_ind = 122; 
twoOmega_signal = squeeze(sum(tmp(:,twoOmega_ind,:),2)); % + conj(squeeze(sum(tmp(:,92:94,:),2)));
twoOmega_nosum = twoOmega_signal; 

% drift compensation
sideband_list = [17]*1240/wavelength-IP(3); % select 16th harmonic of lowest ionization state 
peak_phase = 0; 
for ii=1:1:length(sideband_list)
    [~, index] = min(abs(E - sideband_list(ii)));
    window_center = index; 
    window = 25; 
    histogram_windows = twoOmega_signal((window_center-window):(window_center+window),:); 
    % integrate over window
    peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); %./squeeze(sum(sum(E_SpectraArray,1),2))'; 
end
peak_amp = abs(peak_phase); 
peak_phase = angle(peak_phase);
phase_shift = peak_phase; 

for ii=1:1:length(phase_shift)
   twoOmega_signal(:,ii) = twoOmega_signal(:,ii).*exp(-1j*phase_shift(ii)); 
end
% sum phase matched values
twoOmega_signal = squeeze(sum(twoOmega_signal, 2)); 
plotfun_rabbitspectrum(11:1:21, 810, E, twoOmega_signal, 'twoOmega');

region = [4.65 5.85]; % sideband 14
tolerance = abs(E(2)-E(1)); 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first');
int2w = [int2w sum(abs(twoOmega_signal))]; 
max14 = [max14 max(abs(twoOmega_signal(start:stop)))]; 
int14 = [int14 sum(abs(twoOmega_signal(start:stop)))/sum(abs(twoOmega_signal))]; 

%% PHASE EXTRACTION BY FITTING
% signal = squeeze(sum(twoOmega_0V,2)); 
% signal = twoOmega_0V; 
signal = twoOmega_signal; 

% region = [1.55 2.8]; % sideband 12
% region = [4.65 5.85]; % sideband 14
% region = [7.6 9.15]; % sideband 16
% region = [10.65 12.2]; % sideband 18
% region = [5.95 7.25]; % SB16 5V; 
region = [9.0 10.3]; % SB18 5V; 

% argon
% region = [2 3]; % SB12
% region = [5.2 6.2]; % SB14
% region = [8 9.4]; % SB16
% region = [10.8 12.4]; % SB18

tolerance = abs(E(2)-E(1)); 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), abs(signal(start:stop)), signal(start:stop), 1, 1); 
% save as labeled variables
paramout_original = paramout;  
fval_original = fval

%% CASE RESAMPLING (BOOTSTRAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array
trials = 100; 
% numfiles = size(E_SpectraArray, 3); 
numfiles = numel(phase_shift); 
% array to save fit parameters
paramout_array = cat(2, zeros([size(paramout_gauss), trials]), zeros([size(paramout_gauss,1), 2, trials])); 
fval_array = zeros([1, trials]); 
sample_2w_array = zeros([size(twoOmega_signal,1), trials]); 
sample_ES_array = zeros([size(twoOmega_signal,1), trials]); 
% yfit_array = zeros([numel(xdata) trials]); 

for nn=1:trials
%     resample files with replacement
    sample_list = randsample(numfiles, numfiles, true); 
    sample_ESpectra = data(:,:,sample_list); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(sum(tmp(:,twoOmega_ind,:),2));
    
    phase_offset = angle(twoOmega_signal(start)) - angle(sample_twoOmega(start)); 
    
    for mm=1:1:numfiles
        sample_twoOmega(:,mm) = sample_twoOmega(:,mm).*exp(-1j*phase_shift(sample_list(mm))); 
    end
    sample_twoOmega = sum(sample_twoOmega, 2); 
    sample_2w_array(:,nn) = squeeze(sample_twoOmega); 
    sample_ES_array(:,nn) = squeeze(sum(sum(sample_ESpectra,2),3)); 

    % fit data
    [paramout, fval] = complexfit_section_bootstrap(wavelength, E(start:stop), sample_twoOmega(start:stop), paramout_gauss, paramout_original, 0); 
    paramout_array(:,:,nn) = paramout; 
    fval_array(nn) = fval; 

end 
check_if_done = 'done!'

%% check fit convergence on resampled sets

xin = E(start:stop); 
tmp = zeros([1 100]); 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    p = plot(xin, abs(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))), 'Color', 'c', 'LineWidth', 2); 
%     p = plot(xin, abs(yfit_array(:,ii)), 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plot(xin, abs(Spectrum(xin, paramout_gauss, paramout_original)), 'Color', 'b', 'LineWidth', 2); 
% plot(xin, mean(abs(yfit_array),2), 'Color', 'b', 'LineWidth', 2); 
% scatter(xin, abs(twoOmega_signal(start:stop))/sum(abs(twoOmega_signal(start:stop))), 'bo'); 
title('Fit Amplitude', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Amplitude'); 
hold off; 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    plotphase = unwrap(mod(angle(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))),2*pi)); 
%     plotphase = unwrap(mod(angle(yfit_array(:,ii)),2*pi)); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
    s = plot(xin, unwrap(mod(angle(sample_2w_array(start:stop,ii)),2*pi)), 'cs');
    s.Color(4) = 0.1; 
    
    tmp1 = unwrap(mod(angle(sample_2w_array(start:stop,ii)),2*pi)); 
    tmp(ii) = tmp1(37); 
end

plotphase = unwrap(mod(angle(Spectrum(xin, paramout_gauss, paramout_original)),2*pi)); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
scatter(xin, unwrap(mod(angle(twoOmega_signal(start:stop)),2*pi)), 'bo'); 

plotphase = unwrap(mod(angle(Spectrum(xin, mean(paramout_array(:,1:3,:),3), mean(paramout_array(:,4:end,:),3))), 2*pi)); 
plot(xin, plotphase, 'Color', 'k', 'LineWidth', 2); 
scatter(xin, unwrap(mod(angle(mean(sample_2w_array(start:stop,:),2)), 2*pi)), 'ks'); 

title('Fit Phase', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Phase (referenced to v=6 phase)'); 
hold off; 

%% save

% bootstrap_phase = mean(paramout_array(:,4,:),3);
% bootstrap_slope = mean(paramout_array(:,5,:),3); 
% 
% bootstrap_phase_std = sqrt(trials/(trials-1))*std(paramout_array(:,4,:),0,3);
% bootstrap_slope_std = sqrt(trials/(trials-1))*std(paramout_array(:,5,:),0,3);

% H2_5V_SB16_paramout = paramout_array; 
% H2_5V_SB16_phase = [bootstrap_phase, bootstrap_phase_std]; 
% H2_5V_SB16_slope = [bootstrap_slope, bootstrap_slope_std]; 

H2_3V_SB18_err = [mean(phase_array,2), sqrt(trials/(trials-1))*std(phase_array,0,2)]; 
% H2_5V_SB16_peaks = peak_centers; 

%% CASE RESAMPLING (BOOTSTRAP) for SI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SB = 15; 
% peak_centers = paramout_gauss(:,2); 
peak_centers = mean(E(start:stop)); 
trials = 100; 
numfiles = numel(phase_shift); 
fval_array = zeros([1, trials]); 
sample_2w_array = zeros([size(twoOmega_signal,1), trials]); 
sample_ES_array = zeros([size(twoOmega_signal,1), trials]); 
phase_array = zeros([size(paramout_gauss,1), trials]); 
for nn=1:trials
%     resample files with replacement
    sample_list = randsample(numfiles, numfiles, true); 
    sample_ESpectra = data(:,:,sample_list); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(sum(tmp(:,twoOmega_ind,:),2));
    
    phase_offset = angle(twoOmega_signal(start)) - angle(sample_twoOmega(start)); 
    
    for mm=1:1:numfiles
        sample_twoOmega(:,mm) = sample_twoOmega(:,mm).*exp(-1j*phase_shift(sample_list(mm))); 
    end
    sample_twoOmega = sum(sample_twoOmega, 2); 
    sample_2w_array(:,nn) = squeeze(sample_twoOmega); 
    sample_ES_array(:,nn) = squeeze(sum(sum(sample_ESpectra,2),3)); 

    % fit data
    phase_array(:,nn) = spectral_integration(E(start:stop), abs(sample_twoOmega(start:stop)), angle(sample_twoOmega(start:stop)), peak_centers, 10); 

end 
check_if_done = 'done!'
%% look inside specific vibrational states
% run first three cells of H2_rainbow_rabbitt.m

H2_params = H2_5V_SB16_paramout; 
Ar_phase = Ar_3V_SB16_phase; 
Ar_slope = Ar_3V_SB16_slope; 

xdata = E(start:stop); 
ydata = twoOmega_signal(start:stop); 

figure; hold on; 
yyaxis left; ylabel('amplitude'); 
p = plot(xdata, abs(ydata), 'k', 'DisplayName', '2w osc. amplitude'); 
p.Color(4) = 0.2; 
for ii=1:size(H2_params,1)
    v_energy = mean(H2_params(ii,2,:),3); 
    v_region = v_energy + [-0.14 0.14]; 
    v_start = find(abs(E-v_region(1))<tolerance, 1, 'last'); 
    v_stop = find(abs(E-v_region(2))<tolerance, 1, 'first'); 
    xtheory = E(v_start:v_stop); 
    ydata = signal(v_start:v_stop); 
    ytheory = RABBITT_phase(v_start:v_stop); 
    
    ar_background = ar_slope * (xtheory - mean(xtheory)); 
    
    yyaxis right; ylabel('phase'); 
    tmp = angle(ydata)' - ar_background; 
%     plot(xtheory, angle(ydata), 'c*-', 'HandleVisibility', 'off'); 
    plot(xtheory, tmp, 'bx-', 'HandleVisibility', 'off'); 
    plot(xtheory, ytheory - mean(ytheory) + mean(tmp), 'k--', 'HandleVisibility', 'off'); 

end
xlim([E(start) E(stop)]); 
% plot(100, 100, 'c*-', 'DisplayName', '2w phase no sub'); 
plot(101, 101, 'bx-', 'DisplayName', '2w phase Ar sub'); 
plot(102, 102, 'k--', 'DisplayName', 'theory'); 

xlabel('photoelectron energy (eV)'); 
% title(strcat('SB ',num2str(SB))); 
legend
goodplot(24)


%% directly compare phases to theory
%% prepare k, E, etc. 
% run C_C_DelayPlot.m
% this will also generate perturbative CC models

% trim NaN's (if starting at E=0)
E = E(2:end); 
k = k(2:end); 
CCP = CCP(2:end); 
CCPA = CCPA(2:end); 
CCPAp = CCPAp(2:end); 
%%
nstates = 5; 
glob_phase = 0; 
% phase_data = flipud(cat(3, H2_0V_SB12_phase, ...
%                            H2_3V_SB14_phase, ...
%                            H2_5V_SB16_phase)); 
% phase_data(:,:,3) = flipud(H2_5V_SB16_SI); 
% phase_data(:,:,2) = flipud(H2_3V_SB14_SI); 
% phase_data(:,:,1) = flipud(H2_0V_SB12_SI); 
phase_data = flipud(cat(3, H2_0V_SB12_SI, ...
                           H2_3V_SB14_SI, ...
                           H2_5V_SB16_SI)); 
% phase_data(:,1,3) = phase_data(:,1,3) + (mean(H2_5V_SB16_phase(:,1))+2*pi - mean(H2_5V_SB16_SI(:,1))); 

Ar_phase = [Ar_0V_SB12_phase; ...
            Ar_3V_SB14_phase; ...
            Ar_3V_SB16_phase; ...
            Ar_3V_SB18_phase]; 
Ar_slope = [Ar_0V_SB12_slope; ...
            Ar_3V_SB14_slope; ...
            Ar_3V_SB16_slope; ...
            Ar_3V_SB18_slope]; 

% mean_data = flipud(cat(3, mean(H2_0V_SB12_paramout(1:nstates,:,:),3), ...
%                           mean(H2_3V_SB14_paramout(1:nstates,:,:),3), ...
%                           mean(H2_5V_SB16_paramout(1:nstates,:,:),3))); 
% SB_phase_data = cat(3, ...
%                [mean_data(1:nstates,2,1)'; unwrap(phase_data(1:nstates,1,1))'], ...
%                [mean_data(1:nstates,2,2)'; unwrap(phase_data(1:nstates,1,2))'], ...
%                [mean_data(1:nstates,2,3)'; unwrap(phase_data(1:nstates,1,3))']); 
% SB_phase_error = cat(3, ...
%                [mean_data(1:nstates,2,1)'; phase_data(1:nstates,2,1)'], ...
%                [mean_data(1:nstates,2,2)'; phase_data(1:nstates,2,2)'], ...
%                [mean_data(1:nstates,2,3)'; phase_data(1:nstates,2,3)']); 
% SB_delay_data = SB_phase_data; SB_delay_data(2,:) = SB_delay_data(2,:).*(T_L*1000/2/(2*pi)); %20191009 data processing accidentally added a 120as phase offset between Ar and H2
% SB_delay_error = cat(3, ...
%                [mean_data(1:nstates,2,1)'; phase_data(1:nstates,2,1)'.*(T_L*1000/2/(2*pi))], ...
%                [mean_data(1:nstates,2,2)'; phase_data(1:nstates,2,2)'.*(T_L*1000/2/(2*pi))], ...
%                [mean_data(1:nstates,2,3)'; phase_data(1:nstates,2,3)'.*(T_L*1000/2/(2*pi))]); 

SB_phase_data = cat(3, ...
               [H2_0V_SB12_peaks'; unwrap(H2_0V_SB12_SI(:,1))'], ...
               [H2_3V_SB14_peaks'; unwrap(H2_3V_SB14_SI(:,1))'], ...
               [H2_5V_SB16_peaks'; unwrap(H2_5V_SB16_SI(:,1))']); 
SB_phase_error = cat(3, ...
               [H2_0V_SB12_peaks'; H2_0V_SB12_SI(:,2)'], ...
               [H2_3V_SB14_peaks'; H2_3V_SB14_SI(:,2)'], ...
               [H2_5V_SB16_peaks'; H2_5V_SB16_SI(:,2)']); 
SB_delay_data = SB_phase_data; SB_delay_data(2,:) = SB_delay_data(2,:).*(T_L*1000/2/(2*pi)); %20191009 data processing accidentally added a 120as phase offset between Ar and H2
SB_delay_error = SB_phase_error; SB_delay_error(2,:) = SB_delay_error(2,:).*(T_L*1000/2/(2*pi)); 
           
% Ar_phase(:,1) = Ar_phase(:,1) - pi; 
% Ar_phase(3:4,1) = Ar_phase(3:4,1)+2*pi; 
Ar_phase(:,1) = unwrap(Ar_phase(:,1)); 
Ar_delay = Ar_phase.*(T_L*1000/2/(2*pi)); 
Ar_energy = [12 14 16 18]*1240/810 - 15.763; 

%% get and plot XUV phase

% Ar_theory = 2*cat(2, [Ar_energy(1); -48], Ar_Dahlstrom(:,1:3)); 
Ar_theory = ArTDSE_810.t(1:4); 

figure; errorbar(Ar_energy, Ar_delay(:,1)-pi, Ar_delay(:,2), 'ko', 'DisplayName', 'Argon raw measurement');
hold on; plot(ArTDSE_810.x_Ee(1:4), Ar_theory, 'rv', 'DisplayName', 'Argon TDSE'); 
xlabel('electron kinetic energy'); 
ylabel('delay (as)'); 
xlim([2 12.5]); 
legend; 
goodplot(24); 
hold off; 

XUV_delay = Ar_delay(:,1) - Ar_theory; 
figure; errorbar(Ar_energy, XUV_delay, Ar_delay(:,2), 'ko', 'DisplayName', 'XUV delay'); 
xlabel('electron kinetic energy'); 
ylabel('delay (as)'); 
xlim([2 12.5]); 
legend; 
goodplot(24); 
hold off; 
%% MAKE THE PLOTS
%% compare phase to theory
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

subplot(2, 3, [1 2 3]); hold on; 
errorbar(xdata, plot_data, plot_error, 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee, -H2TDSE_810_140.t, H2TDSE_810_140.err, ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee, -H2TDSE_810_145.t, H2TDSE_810_145.err, ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee, -H2TDSE_810_150.t, H2TDSE_810_150.err, ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
xlabel('electron kinetic energy (eV)'); 
ylabel('delay (as)'); 
legend; 
xlim([1 10]); 
ylim([0 500]); 
goodplot(24); 

subplot(2, 3, 4); hold on; 
errorbar(xdata(1:5), plot_data(1:5), plot_error(1:5), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(1), -H2TDSE_810_140.t(1), H2TDSE_810_140.err(1), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(1), -H2TDSE_810_145.t(1), H2TDSE_810_145.err(1), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(1), -H2TDSE_810_150.t(1), H2TDSE_810_150.err(1), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([1.4 3])
% ylim([-200 -110]); 
goodplot(24); 

subplot(2, 3, 5); hold on; 
errorbar(xdata(6:10), plot_data(6:10), plot_error(6:10), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(2), -H2TDSE_810_140.t(2), H2TDSE_810_140.err(2), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(2), -H2TDSE_810_145.t(2), H2TDSE_810_145.err(2), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(2), -H2TDSE_810_150.t(2), H2TDSE_810_150.err(2), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([4.4 6.2]); 
% ylim([-80-45 -80+45]); 
goodplot(24); 

subplot(2, 3, 6); hold on; 
errorbar(xdata(11:15), plot_data(11:15), plot_error(11:15), 'ko', 'DisplayName', 'H2 measurement'); 
errorbar(H2TDSE_810_140.x_Ee(3), -H2TDSE_810_140.t(3), H2TDSE_810_140.err(3), ...
    'rv', 'DisplayName', 'H2 TDSE 810nm, r=1.40'); 
errorbar(H2TDSE_810_145.x_Ee(3), -H2TDSE_810_145.t(3), H2TDSE_810_145.err(3), ...
    'r*', 'DisplayName', 'H2 TDSE 810nm, r=1.45'); 
errorbar(H2TDSE_810_150.x_Ee(3), -H2TDSE_810_150.t(3), H2TDSE_810_150.err(3), ...
    'rs', 'DisplayName', 'H2 TDSE 810nm, r=1.50'); 
plot(E_CC, -t_interp ./ T_AU, 'r--', 'HandleVisibility', 'off'); 
% xlabel('electron kinetic energy (eV)'); 
% ylabel('delay (as)'); 
% legend; 
xlim([7.5 9.2]); 
% ylim([-58-45 -58+45]); 
goodplot(24); 

%% look at RR
figure; hold on; 
% subplot(1,3,1); hold on; 
%%
% run correct voltage twoOmega analysis
signal = squeeze(sum(twoOmega_signal,2)); 
SB=16; 
rV = 5;
% subplot(1,3,3); hold on; 

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
%     Ar_SB12_slope = Ar_3V_SB12_slope; 
    Ar_SB14_slope = Ar_3V_SB14_slope; 
    Ar_SB16_slope = Ar_3V_SB16_slope; 
    Ar_SB18_slope = Ar_3V_SB18_slope; 
    
%     Ar_SB12_phase = Ar_3V_SB12_phase; 
    Ar_SB14_phase = Ar_3V_SB14_phase; 
    Ar_SB16_phase = Ar_3V_SB16_phase; 
    Ar_SB18_phase = Ar_3V_SB18_phase; 
end
    

if SB==12
    region = [1.6 2.8]; % sideband 12
    ar_slope = Ar_SB12_slope(1); 
    ar_phase = Ar_SB12_phase(1); 
    H2_phase = H2_0V_SB12_phase(1); 
    H2_param = H2_0V_SB12_paramout; 
elseif SB==14
    region = [4.7 6.15]; % sideband 14
    ar_slope = Ar_SB14_slope(1); 
    H2_param = H2_3V_SB14_paramout; 
elseif SB==16
%     region = [7.7548 9.18]; % sideband 16
    region = [5.9 7.3];
    ar_slope = Ar_SB16_slope(1); 
    H2_param = H2_5V_SB16_paramout; 
%     E = E_5V; 
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


yyaxis left; ylabel('amplitude'); 
p = plot(xdata, abs(ydata), 'k'); 
p.Color(4) = 0.2; 

% look inside specific vibrational states
for ii=1:size(H2_param,1)
    v_energy = mean(H2_param(ii,2,:),3); 
    v_region = v_energy + [-0.1 0.1]; 
    v_start = find(abs(E-v_region(1))<tolerance, 1, 'last'); 
    v_stop = find(abs(E-v_region(2))<tolerance, 1, 'first'); 
    xtheory = E(v_start:v_stop); 
%     ytheory = RABBITT_phase(v_start:v_stop); 
    ydata = signal(v_start:v_stop); 

    ar_background = ar_slope * (xtheory - mean(xtheory)); 
%     ar_background = 0; 
    
    yyaxis right; ylabel('phase'); 
    tmp = angle(ydata)' - ar_background; 
%     plot(xtheory, angle(ydata), 'c*-'); 
    plot(xtheory, tmp, 'bx-'); 
%     plot(xtheory, ytheory - mean(ytheory) + mean(tmp), 'k-'); 

end

xlabel('photoelectron energy (eV)'); 
title(strcat('SB ',num2str(SB))); 
goodplot(24)
%% spectra plot
figure; hold on; 
plot(Ebins, sum(sum(ESpectra_0V,2),3)./sum(ESpectra_0V(:)), 'k-', 'DisplayName', '0V'); 
plot(Ebins, sum(sum(ESpectra_3V,2),3)./sum(ESpectra_3V(:)), 'b-', 'DisplayName', '3V'); %'color', [0.25 0.25 0.25]); 
% plot(E, sum(sum(ESpectra_5V,2),3)./sum(ESpectra_5V(:)), 'g-'); %'color', [0.5 0.5 0.5]); 
plot(Ebins+1240/810+0.224, sum(sum(ESpectra_5V,2),3)./sum(ESpectra_5V(:)), 'g-', 'DisplayName', '5V'); %'color', [0.5 0.5 0.5]); 
legend; 
xlabel('photoelectron energy (eV)'); 
% ylabel('normalized counts'); 
xlim([0 18])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
goodplot(24); 

ax1 = gca; 
axl = AddHarmonicAxis(ax1, IP, IP_label, wavelength, 1);

%% carpet plot
h3 = figure; hold on; 
tmp = fft(ESpectra_0V,[],2) .* exp(1j*permute(repmat(peak_phase, [900 1 size(ESpectra_0V,2)]), [1 3 2])); 
tmp = ifft(fftshift(tmp,2),[],2);
carpet = abs(sum(tmp, 3));% - median(median(abs(sum(tmp,3)),1),2); 
% carpet = carpet - repmat(mean(carpet,1),[900 1]);
% carpet = carpet - repmat(XUV_only./sum(XUV_only), [1 223]); 
% carpet = carpet - repmat(mean(carpet,2),[1 223]); 
imagesc(Ebins, stageTimes, carpet'); 
% colormap(map); 
c = colorbar; 
c.Label.String = 'photoelectron yield';
xlim([0 9.5]); 
ylim([-6 -1]*10^(-15))
xlabel('photoelectron energy (eV)'); 
ylabel('IR/XUV time delay (s)'); 
goodplot(24); 

ax1 = gca; 
axl = AddHarmonicAxis(ax1, IP, IP_label, wavelength, 1);

%% cartoon

xdata = -10:0.01:20; 

De = 20; 
a = 0.1; 
re=2; 
morse = De * (exp(-2*a*(xdata-re)) - 2*exp(-a*(xdata-re))); 
figure; hold on; 
plot(xdata, morse, 'k-'); 
ylim([-30 20]); 

figure; hold on; 
for ii=0:50
    E = 5 * (ii+0.5) - (5 * (ii+0.5)).^2/4/De; 
    plot(xdata, E*ones(size(xdata)), 'k-'); 
end

%% functions
function Yout = Spectrum(E, gaussian, p)
    Yout = 0; 
    Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
    if size(p,2) == 1
        Phase = @(x,b,mu) exp(1j .* b); 
        % sum the 2w signal
        for n = 1:size(p,1)
            Amp = gaussian(n,1); 
%             Amp = 1; 
            E0 = gaussian(n,2); 
            wid = gaussian(n,3); 

            b = mod(p(n,1),2*pi); 

            Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,b,E0);
        end
    elseif size(p,2) == 2
        Phase = @(x,b,c,mu) exp(1j .* (b + c.*(x-mu)) ); 
        % sum the 2w signal
%             clist = [0.17, 0.64, 0.64]; 
        for n = 1:size(p,1)
            Amp = gaussian(n,1); 
%             Amp = 1; 
            E0 = gaussian(n,2); 
            wid = gaussian(n,3); 

            b = mod(p(n,1),2*pi);
            c = p(n,2); 
%                 c = clist(n); 


            Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,b,c,E0);
        end
    else
        error('invalid guess input'); 
    end

    yout_mat = [abs(Yout); angle(Yout)]; 

end
