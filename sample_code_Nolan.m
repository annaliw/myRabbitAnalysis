% SAMPLE CODE FOR NOLAN
% load in data and sort into matrix
folderName = 'insert_folder_path/2020_08_25-16Scan/'; % Argon data folder

alternate = [1 1]; % specifies whether data alternates between two samples, which we used to do to try to keep conditions similarbetween samples throughout the experiment. Now we are more confident that our APT is constant through several days. 
wavelength=810; 
% ionization potential for the sample. 
global IP; IP = [15.7596]; % global from some messy attempts at modular code
global IP_label; IP_label = ["Ar 2P 3/2"]; 
% includes the other spin orbit state IP, only use if you're trying to fit both peaks
% global IP; IP = [15.763, 15.763+0.17749]; 
% global IP_label; IP_label = ["Ar 2P 3/2", "Ar 2P 1/2"]; % start with 2

[HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderName, 1, wavelength); % getrawdata makes the array
HistTot_array = HistTot_array(:,:,alternate(1):alternate(2):end); % might need to alternate files in some data sets
% the last scan in every file is XUV only (IR very far away) to get a pure XUV spectrum
XUV_only = sum(XUV_only(:,alternate(1):alternate(2):end),2); 

% I use this figure to find energy calibration input. For now I have just
% given you the energy calibration. 
figure; plot(sum(sum(HistTot_array,2),3)); 

%% load calibration
load(strcat(folderName, 'calibration/Ar_calibration.mat')); 
% returns vector A, which contains coefficients of the tof vs. energy fit
t0=65; % I'll explain this value later, but it's the true t0 for the detector, where a light peak shows up
%% do energy conversion
tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); % make tof array, it's just 1:8000 ns
E_vec = [0 25 900]; % tells the next program to make 0:25 eV array with 900 bins

%convert the ToF spectrum to energy (linear energy scale)
[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % bins tof into E using overlap matrix for "in-between" values
% C=conversion matrix,
% E=new energy array
% OM=overlap matrix, which I save because it's expensive to compute

% E_SpectraArray is the new HistTot_array, converted to energy
% E_Spectra array is 3d: energy x time delay x file
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
Ebins = E; % save for later in case E gets overwritten

%% fft and filter data to get 2w signal
% tmp = fftshift(fft(E_SpectraArray(:,:,:), [], 2), 2); % no mean (DC) subtraction
tmp = fftshift(fft(ifftshift(E_SpectraArray - repmat(mean(E_SpectraArray,2),[1 numel(stageTimes)]),2),[],2),2); % mean subtraction

twoOmega_ind = 133; % this is the location of the 2w bin. Right now I plot sum(abs(tmp,3)) and select it by hand. 
twoOmega_signal = squeeze(sum(tmp(:,twoOmega_ind,:),2)); 
twoOmega_nosum = twoOmega_signal; % in case I need to redo some analysis and I have overwritten twoOmega_signal

% at this point, try plotting the sum of all the twoOmega_signals. Do you
% see much? Then compare after doing the drift compensation. 

%% drift compensation
% Over the course of the experiment, t0 drifts due to beam drift. This is
% like having a global phase drift on each file. 
% The idea here is to find the global phase of each file, and then set
% all the global phases to zero. 
twoOmega_signal = twoOmega_nosum; % reset point
sideband_list = [18]*1240/wavelength-IP(1); % energy of 18th harmonic. I select 18th harmonic of lowest ionization state as global phase reference. 18th harmonic is good because it's out at high energy where there isn't a lot of crazy stuff going on.  

peak_phase = 0; 
for ii=1:1:length(sideband_list) % we can also average multiple sidebands, sometimes I try that. This is why it's a loop over SB's. 
    [~, index] = min(abs(E - sideband_list(ii))); % find index of 18th sideband
    window_center = index; 
    window = 4; 
    histogram_windows = twoOmega_signal((window_center-window):(window_center+window),:); 
    % integrate over window
    % this is the average phase of the 18th harmonic
    peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); 
end
peak_amp = abs(peak_phase); 
peak_phase = angle(peak_phase); % global phase (phase correction) of each file
phase_shift = peak_phase; % variable renamed for no good reason

for ii=1:1:length(phase_shift) % now do phase correction of each file
   twoOmega_signal(:,ii) = twoOmega_signal(:,ii).*exp(-1j*phase_shift(ii)); 
end
% sum phase matched values
twoOmega_signal = squeeze(sum(twoOmega_signal, 2)); 

% now you have a single, summed 2w signal for your data set
% You can plot its amplitude and phase
% You need to do more processing to extract a phase value for each 2w peak
% You can either do spectral averaging (amplitude weighted sum of phase),
% or use fitting to extract a phase value (more useful in
% crowded/overlapped spectra). 