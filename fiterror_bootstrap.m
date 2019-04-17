%%%% bootstrap method for getting errorbars on fit
%%%% case resampling 

% folderName = '/Users/annaliw/code/2018_07_31-16Scan/';
    folderName = '/Users/annaliw/code/KrCO2_scan/'; 
% folderName = '/Users/annaliw/code/NOscan/'; 
alternate = [2 2]; % debug setting
t0=0; 
wavelength=810; 
% global IP; IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206]; 
% global IP_label; IP_label = ["0", "1", "2", "3", "4", "5"]; 
    IP = [13.9000   17.6000   18.0770   19.3760]; 
    IP_label = ["X", "A", "B", "C"]; 
%     IP = [14 14.665]; 
%     IP_label = ["14", "14.665"]; 
%     IP = [15.763]; 
%     IP_label = ['Argon']; 
% IP = [9.553, 16.56, 18.318, 21.722]; 
% IP_label = ['X', 'b', 'A', 'c']; 

[HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderName, 1, wavelength);  
HistTot_array = HistTot_array(:,:,alternate(1):alternate(2):end); % might need to alternate files
% HistTot_array = HistTot_array./sum(HistTot_array(:)); % normalize!
XUV_only = sum(XUV_only(:,alternate(1):alternate(2):end),2); 
%     % do data padding for better fitting
%     pad_data = zeros([size(HistTot_array,1)*10-9, size(HistTot_array,2), size(HistTot_array,3)]); % not sure how the sizing works here
%     for ii=1:1:size(HistTot_array, 2)
%         for jj=1:1:size(HistTot_array, 3)
%             pad_data(:,ii,jj) = interp1(0:1:(size(HistTot_array,1)-1), HistTot_array(:,ii,jj), 0:0.1:(size(HistTot_array,1)-1)); 
%             pad_data(:,ii,jj) = poissrnd(pad_data(:,ii,jj)); 
%         end
%     end
%     HistTot_array = pad_data; 

%% load calibration
% load(strcat(folderName, 'calibration/Ar_calibration.mat')); 
load('/Users/annaliw/code/KrCO2_scan/calibration/Kr_calibration.mat'); 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to do full analysis on the original file. this will also give
% structures to reference for size in the loop. 

% do energy conversion
tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); 
E_vec = [0 20 900]; 

%convert the ToF spectrum to energy (linear energy scale)
[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be used in loop
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 

% Difference signal
%     E_SpectraArray = E_SpectraArray - mean(E_SpectraArray,2); 
% fft data

% normalize on each delay step (dim 2 in array)
norm = sum(E_SpectraArray,1); 
E_SpectraArray = E_SpectraArray./repmat(norm, size(E_SpectraArray, 1), 1, 1); 

%%

% fft and filter
tmp = fftshift(fft(E_SpectraArray, [], 2), 2); 
twoOmega_location = 130; % MAKE THIS AUTOMATICALLY DETECTED
twoOmega_signal = squeeze(tmp(:,twoOmega_location,:));

twoOmega_nosum = twoOmega_signal; 

%% drift compensation
twoOmega_signal = twoOmega_nosum; 
%     phase_shift = interp1(1:length(peak_phase), peak_phase, 1.5:1:length(peak_phase));
sideband_list = [16]*1240/wavelength-IP(3); % select 16th harmonic of lowest ionization state 
peak_phase = 0; 
for i=1:1:length(sideband_list)
    [~, index] = min(abs(E - sideband_list(i)));
    window_center = index; 
    window = 3; 
    histogram_windows = twoOmega_signal((window_center-window):(window_center+window),:); 
    % integrate over window
    peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); 
    % % fft wrt IR delay 
    % peak_fft = fftshift(fft(peak_vol, [], 1), 1); 
    % peak_phase = angle(peak_fft(twoOmega_location, :)); 
end
peak_phase = angle(peak_phase); 

%% more drift comp
phase_shift = peak_phase; 
%     phase_shift = interp1(1:length(peak_phase), peak_phase, 1.5:1:length(peak_phase));
for ii=1:1:length(phase_shift)
   twoOmega_signal(:,ii) = twoOmega_signal(:,ii).*exp(-1j*phase_shift(ii)); 
end

% sum phase matched values
twoOmega_signal = sum(twoOmega_signal, 2); 

%% convert XUV only data
tmp = XUV_only; 
%Convert ToF to Energy using previously calculated Overlap Matrix (OM)
if (numel(tof) == size(tmp,1))
    Counts = zeros( size(tmp,2), numel(E) );
    for ind = 1:size(tmp,2)
        %Counts(ind, :) = OM * [tcounts(2:end,ind);0];
        Counts(ind, :) = OM * circshift(tmp(:,ind),-1);
    end
elseif (numel(tof) == size(tmp,2))
    Counts = zeros( size(tmp,1), numel(E) );
    for ind = 1:size(tmp,1)
        %Counts(ind, :) = OM * [tcounts(ind,2:end)';0];
        Counts(ind, :) = OM * circshift(tmp(ind,:)',-1);
    end
else
    Counts = 0;
end
XUV_only = Counts'; 


%% FIRST FIT of original data set (single sideband)

% H2
% region = [1.65 3.007]; % sideband 12
% region = [4.7 6.1]; % sideband 14
% region = [7.7548 9.18]; % sideband 16
% region = [10.8258 12.2]; % sideband 18

% CO2
% region = [2.8 4.1419]; 
% region = [4.3742 5.7677]; 
% region = [5.8452 7.2387]; 
% region = [7.4 8.7]; 
region = [8.9 10.08]; 
% region = [10.2 11.3]; 
% region = [11.8 15.5]; 

tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), abs(twoOmega_signal(start:stop)), twoOmega_signal(start:stop), 1); 
% save as labeled variables
paramout_original = paramout; 
fval_original = fval; 

%% CASE RESAMPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trials = 100; 
numfiles = size(E_SpectraArray, 3); 
% array to save fit parameters
paramout_array = zeros([size(paramout_original), trials]); 
paramout_array = cat(2, paramout_array, paramout_array); 
fval_array = zeros([1, trials]); 
sample_2w_array = zeros([size(twoOmega_signal,1), trials]); 
sample_ES_array = zeros([size(twoOmega_signal,1), trials]); 

for nn=1:1:trials
%     resample files with replacement
    sample_list = randsample(numfiles, numfiles, true); 
%     sample_list = ones([1 40]); 
    sample_ESpectra = E_SpectraArray(:,:,sample_list); 
    sample_norm = sum(abs(sample_ESpectra),1); 
%     sample_ESpectra = sample_ESpectra./repmat(sample_norm, size(sample_ESpectra, 1), 1, 1); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(sum(tmp(:,130,:),2));
 
%     % restack 2w component
%     sideband_list = [16]*1240/wavelength-IP(3); % select 16th harmonic of lowest ionization state 
%     peak_phase = 0; 
%     for i=1:1:length(sideband_list)
%         [~, index] = min(abs(E - sideband_list(i)));
%         window_center = index; 
%         window = 3; 
%         histogram_windows = sample_twoOmega((window_center-window):(window_center+window),:); 
%         % integrate over window
%         peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); 
%         % % fft wrt IR delay 
%         % peak_fft = fftshift(fft(peak_vol, [], 1), 1); 
%         % peak_phase = angle(peak_fft(twoOmega_location, :)); 
%     end
%     peak_shift = angle(peak_phase); 
    for mm=1:1:length(phase_shift)
        sample_twoOmega(:,mm) = sample_twoOmega(:,mm).*exp(-1j*phase_shift(sample_list(mm))); 
    end
    sample_twoOmega = sum(sample_twoOmega, 2); 
    sample_2w_array(:,nn) = squeeze(sample_twoOmega); 
    sample_ES_array(:,nn) = squeeze(sum(sum(sample_ESpectra,2),3)); 

    % fit data
    [paramout, fval] = complexfit_section_bootstrap(wavelength, E(start:stop), sample_twoOmega(start:stop).', paramout_gauss, paramout_original, 0); 
    paramout_array(:,:,nn) = paramout; 
    fval_array(nn) = fval; 

end 
check_if_done = 'done!'
% figure; plot(E, var(squeeze(sample_ES_array), 0, 2)); 
% figure; plot(E, std(squeeze(sample_ES_array), 0, 2)); 
% figure; plot(E, sum(sum(E_SpectraArray./repmat(sum(E_SpectraArray,1), [size(E_SpectraArray,1) 1]),2),3)); 


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
    sample_twoOmega = squeeze(tmp(:,twoOmega_location,:));
 
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
        % peak_phase = angle(peak_fft(twoOmega_location, :)); 
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
%% check fit convergence on resampled sets

% paramout_jackknife = paramout_array; 
% paramout_jackknife_referenced = paramout_jackknife; 
% paramout_jackknife_referenced(:,2,:) = paramout_jackknife(:,2,:) ...
%     - repmat(paramout_jackknife(end,2,:), [size(paramout_jackknife, 1), 1]);

% phaseloc = 5; 
% paramout_bootstrap = paramout_array; 
% paramout_bootstrap_referenced = paramout_bootstrap; 
% paramout_bootstrap_referenced(:,phaseloc,:) = paramout_bootstrap(:,phaseloc,:) ...
%     - repmat(paramout_bootstrap(end,phaseloc,:), [size(paramout_bootstrap, 1), 1]);
% 
% paramout_array = paramout_bootstrap; 

% tolerance = 0.02; 
% start = find(abs(E-region(1))<tolerance, 1, 'last'); 
% stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

xin = E(start:stop); 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    p = plot(xin, abs(Spectrum(xin, paramout_array(:,:,ii))), 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plot(xin, abs(Spectrum(xin, cat(2, paramout_gauss, paramout_original))), 'Color', 'b', 'LineWidth', 2); 
title('Fit Amplitude', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Amplitude'); 
hold off; 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    plotphase = unwrap(angle(Spectrum(xin, paramout_array(:,:,ii)))); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
paramout_original_referenced = paramout_original; 
paramout_original_referenced(:,2,:) = paramout_original(:,2,:) ...
    - repmat(paramout_original(end,2,:), [size(paramout_original, 1), 1]);
plotphase = angle(Spectrum(xin, cat(2, paramout_gauss, paramout_original))); 
plotphase(1)
% plotphase = plotphase - plotphase(1); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
title('Fit Phase', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Phase (referenced to v=6 phase)'); 
hold off; 

%%
figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    plotphase = abs(sample_2w_array(start:stop,ii)); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plotphase = abs(twoOmega_signal(start:stop)); 
plotphase(1)
% plotphase = plotphase - plotphase(1); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
title('Fit Phase', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('2w Amplitude'); 
hold off; 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    plotphase = unwrap(angle(sample_2w_array(start:stop,ii))); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plotphase = unwrap(angle(twoOmega_signal(start:stop))); 
plotphase(1)
% plotphase = plotphase - plotphase(1); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
title('Fit Phase', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Phase (referenced to v=6 phase)'); 
hold off; 

%% plot jackknife values with error bars
vstates = (1:1:size(paramout_jackknife_referenced, 1)) - 1; 

jackknife_phase = mean(paramout_jackknife_referenced(:,2,:),3);
jackknife_slope = mean(paramout_jackknife_referenced(:,3,:),3);

jackknife_phase_std = sqrt(numfiles-1)*std(paramout_jackknife_referenced(:,2,:),0,3);
jackknife_slope_std = sqrt(numfiles-1)*std(paramout_jackknife_referenced(:,3,:),0,3);

figure; hold on; grid on; 
errorbar(vstates, jackknife_phase, jackknife_phase_std, '-o', ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Color', 'r', ...
    'MarkerSize', 8, 'LineWidth', 2); 
xlim([-0.5 5.5]); 
xlabel('Vibrational State'); ylabel('Phase (referenced to v=6)'); 
title('Phase Extracted by Jackknife/Fitting'); 
hold off; 

figure; hold on; grid on; 
errorbar(vstates, jackknife_slope, jackknife_slope_std, '-o', ...
    'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'Color', 'm', ...
    'MarkerSize', 8, 'LineWidth', 2); 
xlim([-0.5 5.5]); 
xlabel('Vibrational State'); ylabel('Phase Slope'); 
title('Phase Slope Extracted by Jackknife/Fitting'); 
hold off; 

%% plot bootstrap values with error bars
data = paramout_array; 

vstates = (1:1:size(data, 1)) - 1; 

bootstrap_phase = mean(data(:,5,:),3);
bootstrap_slope = mean(data(:,6,:),3);

bootstrap_phase_std = sqrt(trials/(trials-1))*std(data(:,5,:),0,3);
bootstrap_slope_std = sqrt(trials/(trials-1))*std(data(:,6,:),0,3);

X18_param = paramout_array(2,:,:); 
X18_phase = mean(paramout_array(2,5,:),3); 
X18_slope = mean(paramout_array(2,6,:),3);
X18_phase_std = sqrt(trials/(trials-1))*std(data(2,5,:),0,3);
X18_slope_std = sqrt(trials/(trials-1))*std(data(2,6,:),0,3);

% A18_param = paramout_array(2,:,:); 
% A18_phase = mean(paramout_array(2,5,:),3); 
% A18_slope = mean(paramout_array(2,6,:),3);
% A18_phase_std = sqrt(trials/(trials-1))*std(data(2,5,:),0,3);
% A18_slope_std = sqrt(trials/(trials-1))*std(data(2,6,:),0,3);
% 
% B18_param = paramout_array(3,:,:); 
% B18_phase = mean(paramout_array(3,5,:),3); 
% B18_slope = mean(paramout_array(3,6,:),3);
% B18_phase_std = sqrt(trials/(trials-1))*std(data(3,5,:),0,3);
% B18_slope_std = sqrt(trials/(trials-1))*std(data(3,6,:),0,3);

figure; hold on; grid on; 
errorbar(vstates, bootstrap_phase, bootstrap_phase_std, '-o', ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Color', 'r', ...
    'MarkerSize', 4, 'LineWidth', 2); 
xlim([-0.5 5.5]); 
xlabel('Vibrational State'); ylabel('Phase'); 
title('Phase Extracted by Bootstrap/Fitting'); 
hold off; 

figure; hold on; grid on; 
errorbar(vstates, bootstrap_slope, bootstrap_slope_std, '-o', ...
    'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'Color', 'm', ...
    'MarkerSize', 4, 'LineWidth', 2); 
xlim([-0.5 5.5]); 
xlabel('Vibrational State'); ylabel('Phase Slope'); 
title('Phase Slope Extracted by Bootstrap/Fitting'); 
hold off; 




%% functions
function Yout = Spectrum(E, p)
    Yout = 0; 
    Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
    if size(p,2) == 5
        Phase = @(x,a,b,mu) a .* exp(1j .* b); 
        % sum the 2w signal
        for n = 1:size(p,1)
            Amp = p(n,1); 
%             Amp = 1; 
            E0 = p(n,2); 
            wid = p(n,3); 
 
            b = mod(p(n,5),2*pi); 

            Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,b,E0);
        end
    elseif size(p,2) == 6
        Phase = @(x,b,c,mu) exp(1j .* (b + c.*(x-mu)) ); 
        % sum the 2w signal
        for n = 1:size(p,1)
            Amp = p(n,1); 
%             Amp = 1; 
            E0 = p(n,2); 
            wid = p(n,3); 

            b = mod(p(n,5),2*pi);
            c = p(n,6); 

            Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,b,c,E0);
        end
    else
        error('invalid guess input'); 
    end

end








