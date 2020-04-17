%%%% bootstrap method for getting errorbars on fit
%%%% case resampling 

folderName = '/Users/annawang/Documents/data/2018_07_31-16Scan/';
% folderName = '/Users/annaliw/code/2018_04_18-18Scan/'; 
%     folderName = '/Users/annaliw/code/KrCO2_scan/'; 
% folderName = '/Users/annawang/Documents/data/2019_12_13-23Scan/'; % Argon
% folderName = '/Users/annawang/Documents/data/2019_12_14-16Scan/'; % H2 long long long scan
% folderName = '/Users/annawang/Documents/data/2020_03_07-19Scan/'; % H2 

alternate = [2 2]; 
t0=0; 
wavelength=810; 
% global IP; IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206];
% global IP_label; IP_label = ["0", "1", "2", "3", "4", "5"]; % start 1
% global IP; IP = [13.9000   17.6000   18.0770   19.3760]; 
% global IP_label; IP_label = ["X", "A", "B", "C"]; 
% global IP; IP = [14 14.665]; 
% global IP_label; IP_label = ["14", "14.665"]; 
global IP; IP = [15.763]; 
global IP_label; IP_label = ['Argon']; % start with 2
% IP = [9.553, 16.56, 18.318, 21.722]; 
% IP_label = ['X', 'b', 'A', 'c']; 

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
load(strcat(folderName, 'calibration/Ar_calibration.mat')); 
% load('/Users/annaliw/code/2020_01_22-16Scan/calibration/Ar_calibration.mat'); 
% load('/Users/annaliw/code/KrCO2_scan/calibration/Kr_calibration.mat'); 
%% OR redo it

t0=51; 
n = 11:1:21; 
% n=[13 15 17 19]; 
calibEnergy = n*1240/810 - 15.763; 
% tof_peak = [1874 1291 1044 905 805 739 683 642 602 573 545]; 
tof_peak = [1962 1338 1079 936 838 771 715 673 635 603 576]; 
% tof_peak = [1962 1079 838 715 635]; 
calibType = 'none'; 

tof_peak(5)=[]; tof_peak(7)=[]; 
n(5)=[]; n(7)=[]; 
calibEnergy(5)=[]; calibEnergy(7)=[]; 

% config.calibEnergy = n*1240/wavelength - IP; 
config.calibEnergy = calibEnergy; 
config.tofPeaks = tof_peak;   % redo with peak finding
config.IPcal = 15.736; 
% config.IPcal = IP(1); 
config.Plot = 1; 
% calibrate to Kr peaks
A = ECalibrate(t0, n, wavelength, calibType, config); % TO DO: hard set t0 into ECalibrate
calibType = 'Ar'; 
filename = strcat(folderName, 'calibration/', calibType, '_calibration'); 
save(filename, 'A', 'calibEnergy', 'tof_peak', 'wavelength'); 
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
% twoOmega_location = 130; % MAKE THIS AUTOMATICALLY DETECTED
twoOmega_signal = squeeze(tmp(:,130,:)); % + conj(squeeze(mean(tmp(:,94,:),2)));

twoOmega_nosum = twoOmega_signal; 

%% drift compensation
twoOmega_signal = twoOmega_nosum; 
sideband_list = [16]*1240/wavelength-IP(1); % select 16th harmonic of lowest ionization state 
% phase_shift = zeros([1 size(twoOmega_signal,2)]); 

% figure; hold on; 
peak_phase = 0; 
for ii=1:1:length(sideband_list)
    [~, index] = min(abs(E - sideband_list(ii)));
    window_center = index; 
%     window_center = 1301; 
    window = 5; 
    histogram_windows = twoOmega_signal((window_center-window):(window_center+window),:); 
    % integrate over window
    peak_phase = peak_phase + squeeze(sum(histogram_windows, 1))./squeeze(sum(sum(E_SpectraArray,1),2))'; 
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

signal = squeeze(sum(twoOmega_signal,2)); 
% signal = XUV_only; 

% Ar
% region = [2.5 2.8]; % sideband 12
% region = [5.3 6]; % sideband 14
% region = [8 9.5]; % sideband 16
region = [11 12.3]; % sideband 18
% region = [2.4 2.9]; % sideband 12
% region = [5.3 6]; % sideband 14
% region = [8.5 9.2]; % sideband 16
% region = [11 12.3]; % sideband 18

% H2
% region = [1.6 2.8]; % sideband 12
% region = [4.7 6.15]; % sideband 14
% region = [7.7548 9.18]; % sideband 16
% region = [7.72 9.18]; % sideband 16
% region = [10.8258 12.2]; % sideband 18
% region = [0.2 1.8]; % harmonic 11

% CO2
% region = [2.8 4.1419]; 
% region = [4.3742 5.7677]; 
% region = [5.8452 7.2387]; 
% region = [7.4 8.7]; 
% region = [8.9 10.08]; 
% region = [10.2 11.7]; 
% region = [11.8 15.5]; 

% Kr
% region = [3.4 4.6]; 
% region = [12.5 14]; 
% region = [8 9.2]; 

tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

% subtract_floor = 9.7; 
subtract_floor = 0; 
if subtract_floor == 0
    tmp_signal = signal; 
else
    tmp = [abs(signal), angle(signal)]';
    tmp(1,:) = tmp(1,:) - 0.9*abs(signal(find(abs(E-subtract_floor)<tolerance, 1, 'last'))); 
%     tmp(1,:) = tmp(1,:) - 0; 
    tmp_signal = tmp(1,:) .* exp(1j*tmp(2,:));
    tmp_signal = tmp_signal.'; 
end


[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), abs(tmp_signal(start:stop)), tmp_signal(start:stop), 1, 1); 
% save as labeled variables
paramout_original = paramout; 
fval_original = fval

%% CASE RESAMPLING (BOOTSTRAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trials = 100; 
numfiles = size(E_SpectraArray, 3); 
% array to save fit parameters
paramout_array = cat(2, zeros([size(paramout_gauss), trials]), zeros([size(paramout_original), trials])); 
fval_array = zeros([1, trials]); 
sample_2w_array = zeros([size(twoOmega_signal,1), trials]); 
sample_ES_array = zeros([size(twoOmega_signal,1), trials]); 

for nn=1:1:trials
%     resample files with replacement
    sample_list = randsample(numfiles, numfiles, true); 
    sample_ESpectra = E_SpectraArray(:,:,sample_list); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(sum(tmp(:,130,:),2));
    for mm=1:1:length(phase_shift)
        sample_twoOmega(:,mm) = sample_twoOmega(:,mm).*exp(-1j*phase_shift(sample_list(mm))); 
    end
    sample_twoOmega = sum(sample_twoOmega, 2); 
    
    if subtract_floor == 0
        tmp_signal = sample_twoOmega; 
    else
        tmp = [abs(sample_twoOmega), angle(sample_twoOmega)]';
        tmp(1,:) = tmp(1,:) - 0.9*abs(sample_twoOmega(find(abs(E-subtract_floor)<tolerance, 1, 'last'))); 
%         tmp(1,:) = tmp(1,:)-0; 
        tmp_signal = tmp(1,:) .* exp(1j*tmp(2,:));
        tmp_signal = tmp_signal.'; 
    end
    
    sample_2w_array(:,nn) = squeeze(sample_twoOmega); 
    sample_ES_array(:,nn) = squeeze(sum(sum(sample_ESpectra,2),3)); 

    % fit data
    [paramout, fval] = complexfit_section_bootstrap(wavelength, E(start:stop), tmp_signal(start:stop), paramout_gauss, paramout_original, 0); 
    paramout_array(:,:,nn) = paramout; 
    fval_array(nn) = fval; 

end 
check_if_done = 'done!'

%% save
data = paramout_array; 

bootstrap_phase = mean(data(:,4,:),3);
bootstrap_slope = mean(data(:,5,:),3); 

bootstrap_phase_std = sqrt(trials/(trials-1))*std(data(:,4,:),0,3);
bootstrap_slope_std = sqrt(trials/(trials-1))*std(data(:,5,:),0,3);

Ar_SB18_paramout = data; 
Ar_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
Ar_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
% H2_SB12_paramout = data; 
% H2_SB12_phase = [bootstrap_phase, bootstrap_phase_std]; 
% H2_SB12_slope = [bootstrap_slope, bootstrap_slope_std]; 

%% check fit convergence on resampled sets

xin = E(start:stop); 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    p = plot(xin, abs(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))), 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plot(xin, abs(Spectrum(xin, paramout_gauss, paramout_original)), 'Color', 'b', 'LineWidth', 2); 
% scatter(xin, abs(twoOmega_signal(start:stop))/sum(abs(twoOmega_signal(start:stop))), 'b'); 
title('Fit Amplitude', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Amplitude'); 
hold off; 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    plotphase = unwrap(mod(angle(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))),2*pi)); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plotphase = unwrap(mod(angle(Spectrum(xin, paramout_gauss, paramout_original)),2*pi)); 
plotphase(1)
% plotphase = plotphase - plotphase(1); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
scatter(xin, unwrap(mod(angle(twoOmega_signal(start:stop)),2*pi)), 'b'); 
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








