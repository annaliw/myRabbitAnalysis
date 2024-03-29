%%%% bootstrap method for getting errorbars on fit
%%%% case resampling 

% folderName = '/Users/annawang/Documents/data/2018_07_31-16Scan/';
% folderName = '/Users/annawang/Documents/data/2018_04_18-18Scan/'; 
% folderName = '/Users/annawang/Documents/data/2018_07_27-17Scan/'; % NO
% folderName = '/Users/annawang/Documents/data/2018_07_23-19Scan/'; % probably CO2, use Kr for calibration (start 1)
%     folderName = '/Users/annaliw/code/KrCO2_scan/'; 
% folderName = '/Users/annawang/Documents/data/2019_12_13-23Scan/'; % Argon
% folderName = '/Users/annawang/Documents/data/2019_12_14-16Scan/'; % H2 long long long scan
% folderName = '/Users/annawang/Documents/data/2020_03_07-19Scan/'; % H2 
% folderName = '/Users/annawang/Documents/data/NOscan'; 

% % James' scans: 
% folderName = '/Users/annawang/Documents/data/2020_03_15-11Scan_Ar/'; % Ar 
% % folderName = '/Users/annawang/Documents/data/2020_03_15-13Scan_H2/'; % H2 
% % folderName = '/Users/annawang/Documents/data/2020_03_14-21Scan_H2/'; % H2 

folderName = '/Users/annawang/Documents/data/2020_08_25-16Scan/'; % Argon long
% folderName = '/Users/annawang/Documents/data/2020_08_27-13Scan/'; % H2 short
% folderName = '/Users/annawang/Documents/data/2020_08_27-15Scan/'; % H2 0V
% folderName = '/Users/annawang/Documents/data/2020_08_28-15Scan/'; % H2 3V
% folderName = '/Users/annawang/Documents/data/2020_08_30-15Scan/'; % H2 5V
% folderName = '/Users/annawang/Documents/data/2020_08_31-18Scan/'; % H2 8.5V
% folderName = '/Users/annawang/Documents/data/2020_09_01-11Scan/'; % Argon 0V
% folderName = '/Users/annawang/Documents/data/2020_09_01-12Scan/'; % Argon 3V
% folderName = '/Users/annawang/Documents/data/2020_09_01-14Scan/'; % Argon 8.5V


alternate = [1 1]; 
wavelength=810; 
% global IP; IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206];
% global IP_label; IP_label = ["0", "1", "2", "3", "4", "5"]; % start 1
% global IP; IP = [13.776   17.7   18.0770   19.3760]; % CO2
% global IP_label; IP_label = ["X", "A", "B", "C"]; 
% global IP; IP = fliplr([14 14.665]); % Krypton
% global IP_label; IP_label = fliplr(["14", "14.665"]); 
global IP; IP = [15.7596];
global IP_label; IP_label = ["Ar 2P 3/2"]; % start with 2
% global IP; IP = [15.763, 15.763+0.17749];
% global IP_label; IP_label = ["Ar 2P 3/2", "Ar 2P 1/2"]; % start with 2
% IP = [9.553, 16.56, 18.318, 21.722]; % NO
% global IP; IP = [9.55, 15.66305, 16.55995, 16.87514, 17.59988, 17.81864, 18.07566, 18.32547, 21.72947, 22.73139]; 
% global IP_label; IP_label = ["X", "a", "b", "w", "bp", "Ap", "W", "A", "Bc", "BpB"]; 
% global IP; IP = [9.55, 16.55995, 16.87514, 17.59988, 18.32547, 21.72947]; 
% global IP_label; IP_label = ["X", "b", "w", "bp", "A", "Bc"]; 

[HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderName, 1, wavelength);  
HistTot_array = HistTot_array(:,:,alternate(1):alternate(2):end); % might need to alternate files
% HistTot_array = HistTot_array./sum(HistTot_array(:)); % normalize!
XUV_only = sum(XUV_only(:,alternate(1):alternate(2):end),2); 
XUV_only_raw = XUV_only; 
%     % do data padding for better fitting
%     pad_data = zeros([size(HistTot_array,1)*10-9, size(HistTot_array,2), size(HistTot_array,3)]); % not sure how the sizing works here
%     for ii=1:1:size(HistTot_array, 2)
%         for jj=1:1:size(HistTot_array, 3)
%             pad_data(:,ii,jj) = interp1(0:1:(size(HistTot_array,1)-1), HistTot_array(:,ii,jj), 0:0.1:(size(HistTot_array,1)-1)); 
%             pad_data(:,ii,jj) = poissrnd(pad_data(:,ii,jj)); 
%         end
%     end
%     HistTot_array = pad_data; 
figure; plot(sum(sum(HistTot_array,2),3)); 
%% load calibration
% load(strcat(folderName, 'calibration/Ar_calibration.mat')); 
% load('/Users/annawang/Documents/data/2018_07_27-17Scan/calibration/Kr_calibration.mat'); 
load('/Users/annawang/Documents/data/2018_07_23-19Scan/calibration/Kr_calibration.mat'); 
% load('/Users/annaliw/code/2020_01_22-16Scan/calibration/Ar_calibration.mat'); 
% load('/Users/annaliw/code/KrCO2_scan/calibration/Kr_calibration.mat'); 
%% OR redo it

t0=25; 
wavelength = 810; 
n = 10:1:21; 
calibEnergy = repmat(n, [size(IP,2) size(IP,1)])*1240/wavelength - repmat(IP, [size(n,2), size(n,1)]).'; 
% calibEnergy = calibEnergy(:).'; 
% calibEnergy = n*1240/wavelength - IP(3); 
% tof_peak = [2107 1354 1077 929 821 750 695 653 614 581 554]; 
tof_peak = [2614 1429 1118 947 839 762 705 657 619 585 559 534; ...
            1828 1261 1034 894 803 736 684 639 604 573 548 525]; 
% tof_peak = tof_peak(:).'; 
calibType = 'Kr'; 

% tof_peak(5)=[]; tof_peak(7)=[]; 
% n(5)=[]; n(7)=[]; 
% calibEnergy(5)=[]; calibEnergy(7)=[]; 

% config.calibEnergy = n*1240/wavelength - IP; 
config.calibEnergy = calibEnergy; 
config.tofPeaks = tof_peak;   % redo with peak finding
config.IPcal = IP; 
% config.IPcal = IP(1); 
config.Plot = 1; 
% calibrate to Kr peaks
A = ECalibrate(t0, n, wavelength, calibType, config); % TO DO: hard set t0 into ECalibrate
filename = strcat(folderName, 'calibration/', calibType, '_calibration'); 
save(filename, 'A', 'calibEnergy', 'tof_peak', 'wavelength'); 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to do full analysis on the original file. this will also give
% structures to reference for size in the loop. 

% do energy conversion
tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); 
E_vec = [0 25 900]; 

%convert the ToF spectrum to energy (linear energy scale)
[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be used in loop
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
Ebins = E; % save for later in case E gets overwritten

% Difference signal
%     E_SpectraArray = E_SpectraArray - mean(E_SpectraArray,2); 
% fft data

% normalize on each delay step (dim 2 in array)
norm = sum(E_SpectraArray,1); 
% E_SpectraArray = E_SpectraArray./repmat(norm, size(E_SpectraArray, 1), 1, 1); 

%%
% % subtract out slope
% tmp = mean(mean(E_SpectraArray(217:270,:,:),1),3); 
% p = polyfit(stageTimes, tmp, 1); 
% 
% tmp = E_SpectraArray; 
% pmat = zeros([numel(E), numel(stageTimes)]); 
% for ii=1:numel(E)
%     p = polyfit(stageTimes, squeeze(sum(E_SpectraArray(ii,:,:),3)), 1); 
%     pmat(ii,:) = polyval(p,stageTimes); 
%     tmp(ii,:,:) = E_SpectraArray(ii,:,:) - repmat(polyval(p,stageTimes), [1 1 size(E_SpectraArray,3)])/size(E_SpectraArray,3); 
% end
% tmp = tmp * (stageTimes(2)-stageTimes(1)); 
% tmp = fftshift(fft(ifftshift(tmp,2),[],2),2); 

% fft and filter
% tmp = fftshift(fft(E_SpectraArray(:,1:201,:), [], 2), 2); 
tmp = fftshift(fft(ifftshift(E_SpectraArray - repmat(mean(E_SpectraArray,2),[1 numel(stageTimes)]),2),[],2),2); 

% tmp = fftshift(fft(HistTot_array,[],2),2); 
twoOmega_ind = 130; % MAKE THIS AUTOMATICALLY DETECTED
% twoOmega_ind = 157; 
twoOmega_signal = squeeze(sum(tmp(:,twoOmega_ind,:),2)); % + conj(squeeze(sum(tmp(:,92:94,:),2)));

twoOmega_nosum = twoOmega_signal; 

%% drift compensation
twoOmega_signal = twoOmega_nosum; 
sideband_list = [18]*1240/wavelength-IP(1); % select 16th harmonic of lowest ionization state 
% phase_shift = zeros([1 size(twoOmega_signal,2)]); 

% figure; hold on; 
peak_phase = 0; 
for ii=1:1:length(sideband_list)
    [~, index] = min(abs(E - sideband_list(ii)));
    window_center = index; 
%     window_center = 1301; 
    window = 2; 
    histogram_windows = twoOmega_signal((window_center-window):(window_center+window),:); 
    % integrate over window
    peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); %./squeeze(sum(sum(E_SpectraArray,1),2))'; 
end
peak_amp = abs(peak_phase); 
peak_phase = angle(peak_phase);

% more drift comp
% ptmp = polyfit(1:numel(peak_phase), peak_phase, 2); 
% phase_shift = polyval(ptmp, 1:numel(peak_phase));
% plot(unwrap(peak_phase), 'o-'); 
% plot(1:0.1:numel(peak_phase), polyval(ptmp, 1:0.1:numel(peak_phase))); 
phase_shift = peak_phase; 

for ii=1:1:length(phase_shift)
% for ii=1:1:60
   twoOmega_signal(:,ii) = twoOmega_signal(:,ii).*exp(-1j*phase_shift(ii)); 
end
% sum phase matched values
twoOmega_signal = squeeze(sum(twoOmega_signal, 2)); 

%% w space stage drift comp

phase_shift = zeros([2 size(tmp,3)]); 
sideband_list = [12 14 16 18]; 

% for jj=1:size(tmp,3)
for jj=1:1
    jj=10; 
    [~, index] = min(abs(E - (16*1240/wavelength-IP(1))));
    window_center = index; 
    window = 20; 
    wtmp = squeeze(sum(tmp((window_center-window):(window_center+window),:,jj),1)); 
    
    % fit line to w space
    wcenter = find(abs(freqAxis+0.34) < 0.007,1); 
    window = 2; 
    p = polyfit(freqAxis((wcenter-window):(wcenter+window)), unwrap(angle(wtmp((wcenter-window):(wcenter+window)))), 1); 
    phase_shift(1, jj) = p(1); 
    phase_shift(2, jj) = p(2); 
%     phase_shift(3, jj) = p(3); 
    
    figure; hold on; 
    yyaxis left; 
    plot(freqAxis((wcenter-4):(wcenter+4)), abs(wtmp((wcenter-4):(wcenter+4)))); 
    yyaxis right; 
    plot(freqAxis((wcenter-4):(wcenter+4)), unwrap(angle(wtmp((wcenter-4):(wcenter+4)))), 'o-'); 
    plot(freqAxis((wcenter-4):(wcenter+4)), polyval(p, freqAxis((wcenter-4):(wcenter+4)))); 
end

% drift comp
w_signal = zeros([numel(sideband_list) size(tmp,2)]); 
for ii=numel(sideband_list)
    [~, index] = min(abs(E - (sideband_list(ii)*1240/wavelength-IP(1))));
    window_center = index; 
    window = 20; 
    wtmp = unwrap(angle(squeeze(sum(tmp((window_center-window):(window_center+window),:,:),1)))); 
    
    for jj=1:size(tmp,3)
        wtmp(:,jj) = wtmp(:,jj) .* exp(-1j * polyval([phase_shift(1,jj), phase_shift(2,jj)], freqAxis')); 
    end
    w_signal(ii,:) = squeeze(sum(wtmp,2)); 
end

%% E integrated drift comp check
% FT_SB = squeeze(sum(tmp(76:128,:,:),1)); % SB12
% FT_SB = squeeze(sum(tmp(215:269,:,:),1)); % SB14
FT_SB = squeeze(sum(tmp(354:417,:,:),1)); % SB16
peak_phase = angle(FT_SB(92,:)); 
% figure; hold on; 
% plot(unwrap(peak_phase), 'o-')

for ii=1:numel(peak_phase)
    FT_SB(:,ii) = FT_SB(:,ii) .* exp(-1j * peak_phase(ii)); 
end
FT_SB = squeeze(sum(FT_SB,2)); 
IFT_SB = fftshift(ifft(ifft(FT_SB))); 
figure; hold on; 
plot(stageTimes, real(IFT_SB)); 
%% convert XUV only data
tmp = XUV_only_raw; 

%Convert ToF to Energy using previously calculated Overlap Matrix (OM)
if (numel(tof) == size(tmp,1))
    Counts = zeros( size(tmp,2), numel(Ebins) );
    for ind = 1:size(tmp,2)
        %Counts(ind, :) = OM * [tcounts(2:end,ind);0];
        Counts(ind, :) = OM * circshift(tmp(:,ind),-1);
    end
elseif (numel(tof) == size(tmp,2))
    Counts = zeros( size(tmp,1), numel(Ebins) );
    for ind = 1:size(tmp,1)
        %Counts(ind, :) = OM * [tcounts(ind,2:end)';0];
        Counts(ind, :) = OM * circshift(tmp(ind,:)',-1);
    end
else
    Counts = 0;
end
XUV_only = Counts'; 


%% FIRST FIT of original data set (single sideband)

% signal = squeeze(sum(twoOmega_signal,2)); 
signal = twoOmega_signal; 
% signal = XUV_only; 

% Ar
% region = [2.44 2.8]; % sideband 12
% region = [5.4 5.85]; % sideband 14
% region = [8.2 9.2]; % sideband 16
% region = [11.3 12.4]; % sideband 18

% region = [2.1 3.1]; % sideband 12
% region = [5.1 6.5]; % sideband 14
% region = [8.2 9.4]; % sideband 16
% region = [10.8 12.5]; % sideband 18

% H2
% region = [1.5 2.75]; % sideband 12
% region = [4.5 5.85]; % sideband 14
% region = [7.65 9]; % sideband 16
% region = [10.7 11.9]; % sideband 18
% region = [1.77 3.15]; % SB14 3V
% region = [4.8 6.3]; % SB16 3V
% region = [3 4.5] % SB16 5V
% region = [6.2 7.7]; % SB18 5V
% region = [3 4.48]; % SB18 8.5V

% CO2
% region = [2.6 4.1419]; 
% region = [4.3742 5.7677]; 
% region = [5.8452 7.2387]; 
% region = [7.4 8.7]; 
% region = [8.9 10.08]; 
% region = [10.2 11.7]; 
% region = [11.8 15.5]; 

% Kr
% region = [3.4 5]; % SB12
% region = [6.2 7.9]; % SB14
% region = [9.3 10.8]; % SB16
region = [12.5 14]; %18
% region = [15.65 17.2]; % 20

% region = []

% NO
% region = [2.6 3.5]; IP = [8*1240/wavelength-3.15 13*1240/wavelength-3.4 13*1240/wavelength-2.85 18*1240/wavelength-2.5]; IP_label = ["other", "other", "X SB8", "other"]; 
% region = [5.7 6.6]; IP = [10*1240/wavelength-6 15*1240/wavelength-6.4]; IP_label = ["X SB10", "other"]; 
% region = [8.2 9.5]; IP = [12*1240/wavelength-8.75 17*1240/wavelength-8.5 17*1240/wavelength-9.4]; IP_label = ["X SB12", "other", "other"]; % SB12
% region = [10.6 12]; IP = [14*1240/wavelength-11.7 19*1240/wavelength-10.85]; IP_label = ["X SB14", "other"];% SB14
% region = [14 15.5]; IP = [16*1240/wavelength-14.9]; IP_label = ["X SB16"];% SB16
% region = [17.4 18.6]; IP = [18*1240/wavelength-18.1];  IP_label = ["X SB18"]; % SB18
% region = [10.75 12.2]; 
% region = [20.5 22]; IP = [20*1240/wavelength-21]; IP_label=["X SB20"]; 

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
% run poly2phase_20200528.m to make sure this works and initialize poly
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
    sample_ESpectra = E_SpectraArray(:,:,sample_list); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(sum(tmp(:,twoOmega_ind,:),2));
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

%     xdata = E(start:stop); 
%     ydata = tmp_signal(start:stop)'; 
%     y_amp = zeros([size(paramout_gauss,1) numel(xdata)]); 
%     polytmp = zeros([size(paramout_gauss,1) pdeg+1]); % pdeg is defined in poly2phase_20200528.m
%     for ii=1:size(paramout_gauss,1)
%     % for ii=1:2
%         y_amp(ii,:) = paramout_gauss(ii,1) .* exp(-(xdata-paramout_gauss(ii,2)).^2/(2*paramout_gauss(ii,3).^2)); 
% 
%         window = 4; 
%         [cval,center] = min(abs(xdata - paramout_gauss(ii,2))); 
%         if center-window < 1
%             ind1 = 1; 
%         else
%             ind1 = center-window;
%         end
%         if center+window > numel(xdata)
%             ind2 = numel(xdata); 
%         else
%             ind2 = center+window; 
%         end
%         polytmp(ii,:) = polyfit(xdata(ind1:ind2)-xdata(center), angle(ydata(ind1:ind2)), pdeg); % pdeg is defined in poly2phase_20200529.m
%         paramout_array(ii,1:3,nn) = paramout_gauss(ii,:); 
%         paramout_array(ii,4:end,nn) = polytmp(ii,:); 
% 
%         yfit_array(:,nn) = yfit_array(:,nn) + y_amp(ii,:)' .* exp(1j * polyval(polytmp(ii,:), xdata-xdata(center))'); 
%     end
    
%     paramout_array(:,:,nn) = paramout; 

end 
check_if_done = 'done!'

%% save
data = paramout_array; 

bootstrap_phase = mean(mod(data(:,4,:),2*pi),3);
bootstrap_slope = mean(data(:,5,:),3); 

bootstrap_phase_std = sqrt(trials/(trials-1))*std(mod(data(:,4,:),2*pi),0,3);
bootstrap_slope_std = sqrt(trials/(trials-1))*std(data(:,5,:),0,3);

% Ar_SB18_paramout = data; 
% Ar_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
% Ar_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
% H2_SB18_paramout = data; 
% H2_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
% H2_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
Kr_SB18_paramout = data; 
Kr_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
Kr_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
% tmp = [bootstrap_phase, bootstrap_phase_std]; 

% bootstrap_d1 = mean(data(:,end,:),3);
% bootstrap_d2 = mean(data(:,end-1,:),3); 
% % bootstrap_d3 = mean(data(:,end-2,:),3); 
% 
% bootstrap_d1_std = sqrt(trials/(trials-1))*std(data(:,end,:),0,3);
% bootstrap_d2_std = sqrt(trials/(trials-1))*std(data(:,end-1,:),0,3);
% % bootstrap_d3_std = sqrt(trials/(trials-1))*std(data(:,end-2,:),0,3);
% 
% Ar_SB12_paramout = data; 
% Ar_SB12_phase = cat(3, [bootstrap_d1, bootstrap_d1_std], ...
%                        [bootstrap_d2, bootstrap_d2_std]); 
% %                        [bootstrap_d3, bootstrap_d3_std]); 
% % H2_SB14_paramout = data; 
% % H2_SB14_phase = cat(3, [bootstrap_d1, bootstrap_d1_std], ...
% %                        [bootstrap_d2, bootstrap_d2_std], ...
% %                        [bootstrap_d3, bootstrap_d3_std]); 


%% check fit convergence on resampled sets

xin = E(start:stop); 

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
    plotphase = mod(angle(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))),2*pi); 
%     plotphase = unwrap(mod(angle(yfit_array(:,ii)),2*pi)); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
    s = plot(xin, mod(angle(sample_2w_array(start:stop,ii)),2*pi), 'cs');
    s.Color(4) = 0.1; 
end
plotphase = unwrap(mod(angle(Spectrum(xin, paramout_gauss, paramout_original)),2*pi)); 
% plotphase = unwrap(mod(mean(angle(yfit_array),2),2*pi)); 
plotphase(1)
% plotphase = plotphase - plotphase(1); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
plot(xin, unwrap(mod(angle(Spectrum(xin, mean(paramout_array(:,1:3,:),3), mean(paramout_array(:,4:end,:),3))),2*pi)), 'b--', 'Linewidth', 2);
scatter(xin, unwrap(mod(angle(twoOmega_signal(start:stop)),2*pi)), 'bo'); 
title('Fit Phase', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Phase (referenced to v=6 phase)'); 
hold off; 

%% JACKKNIFE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numfiles = size(E_SpectraArray, 3); 
% array to save fit parameters
paramout_array = zeros([size(paramout_original), numfiles]); 
fval_array = zeros([1, trials]); 
sample_2w_array = zeros([size(twoOmega_signal), numfiles]); 
sample_ES_array = zeros([size(twoOmega_signal), numfiles]); 

for nn=1:1:numfiles
    
    sample_ESpectra = E_SpectraArray; 
    sample_ESpectra(:,:,nn) = []; % cut out file nn
    sample_norm = sum(abs(sample_ESpectra),1); 
    sample_ESpectra = sample_ESpectra./repmat(sample_norm, size(sample_ESpectra, 1), 1, 1); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(tmp(:,twoOmega_ind,:));
 
    % restack 2w component
    sideband_list = [12, 14, 16, 18]*1240/wavelength-IP(3); % select 16th harmonic of lowest ionization state 
    peak_phase = 0; 
    for i=1:1:length(sideband_list)
        [~, index] = min(abs(E - sideband_list(i)));
        window_center = index; 
        window = 3; 
        histogram_windows = sample_twoOmega((window_center-window):(window_center+window),:); 
        % integrate over window
        peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); 
        % % fft wrt IR delay 
        % peak_fft = fftshift(fft(peak_vol, [], 1), 1); 
        % peak_phase = angle(peak_fft(twoOmega_ind, :)); 
    end
    peak_phase = angle(peak_phase); 
    for mm=1:1:(numfiles-1)
        sample_twoOmega(:,mm) = sample_twoOmega(:,mm).*exp(-1j*peak_phase(mm)); 
    end
    sample_twoOmega = sum(sample_twoOmega, 2); 
    sample_2w_array(:,nn) = sample_twoOmega; 
    sample_ES_array(:,nn) = sum(sum(sample_ESpectra,2),3); 

    % fit data
    [paramout, fval] = complexfit_section_bootstrap(wavelength, E(start:stop), sample_twoOmega(start:stop)', paramout_gauss, paramout_original, 0); 
    paramout_array(:,:,nn) = paramout; 
    fval_array(nn) = fval; 

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








