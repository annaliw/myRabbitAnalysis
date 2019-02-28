%% load calibration data
plotting = 1; % debug setting
    
folderName = '/Users/annaliw/code/2018_07_31-16Scan/'; % debug setting
alternate = [2 2]; % debug setting
t0=0; 
wavelength=810; 
IP = [15.763]; 
IP_label = ['Argon']; 
n = 11:1:21; 
E_vec = [0 20 900]; 

[HistTot_array, stageTimes, freqAxis] = getrawdata(folderName, 1, wavelength);  
HistTot_array = HistTot_array(:,:,alternate(1):alternate(2):end); % might need to alternate files
HistTot_array = HistTot_array./sum(HistTot_array(:)); % normalize!

%% plot and find peaks

figure; plot(sum(sum(HistTot_array,2),3)); 

hw = 1240/wavelength; 
calibEnergy = n*hw - IP; 
tof_peak = fliplr([577 603 634 673 716 771 838 936 1079 1335 1962]); 

%% find energy calibration
calibType = 'Ar'; 
n = 11:1:21; 

config.calibEnergy = n*1240/wavelength - IP; 
config.tofPeaks = tof_peak;   % redo with peak finding
config.IPcal = IP; 
config.Plot = 1; 
% calibrate to Kr peaks
A = ECalibrate(t0, n, wavelength, calibType, config); % TO DO: hard set t0 into ECalibrate

%% stage drift compensation

% do energy conversion
tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); 
%convert the ToF spectrum to energy (linear energy scale)
[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be used in loop
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
% Difference signal
%     E_SpectraArray = E_SpectraArray - mean(E_SpectraArray,2); 
% fft data
E_SpectraArray = fftshift(fft(E_SpectraArray, [], 2), 2); 
twoOmega_location = 130; % MAKE THIS AUTOMATICALLY DETECTED
twoOmega_signal = squeeze(E_SpectraArray(:,twoOmega_location,:));

% only restack the 2w data
% cut out window of ~6 tof bins
% sideband_list = [(12:2:18)*1240/810 - 13.778, (12:2:18)*1240/810 - 17.706];
sideband_list = [14, 16, 18]*1240/wavelength-IP; % select 16th harmonic of lowest ionization state 
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

%% save values
filename = strcat(folderName, 'calibration/', calibType, '_calibration'); 
save(filename, 'A', 'peak_phase', 'calibEnergy', 'tof_peak', 'wavelength'); 