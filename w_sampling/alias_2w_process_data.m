% folderName = '/Users/annaliw/code/2019_12_18-longscan/';
folderName = '/Users/annaliw/code/2019_12_19-18Scan/';

alternate = [1 1]; 
t0=0; 
wavelength=810; 
global IP; IP = [15.763]; 
global IP_label; IP_label = ['Argon']; 

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
% load('/Users/annaliw/code/2019_12_13-17Scan/calibration/Ar_calibration.mat'); 
% load('/Users/annaliw/code/KrCO2_scan/calibration/Kr_calibration.mat'); 
%% OR redo it

t0=0; 
n = 11:1:19; 
calibEnergy = n*1240/810 - IP(1); 
tof_peak = [1838 1275 1031 892 798 732 680 638 600]; 
calibType = 'none'; 

% config.calibEnergy = n*1240/wavelength - IP; 
config.calibEnergy = calibEnergy; 
config.tofPeaks = tof_peak;   % redo with peak finding
config.IPcal = 15.736; 
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
E_SpectraArray = E_SpectraArray./repmat(norm, size(E_SpectraArray, 1), 1, 1); 

%% FT, filter, process
E_SpectraArray_FFT = fftshift(fft(E_SpectraArray(:,224:end,:), [], 2), 2); 

num_timepts = size(E_SpectraArray,2)-223; 
dt = (stageTimes(225) - stageTimes(224));  
dw =  1/num_timepts/dt;
wrange = 1/dt/2;
fftaxis = (-wrange+dw):dw:wrange; % in Hz
fftaxis = fftaxis * 810E-9 / 2.998E8; % in w

%%
% find 0.5w, where we expect 2w signal to be aliased
w_center_ind = find(abs(fftaxis - 0.3746) < 0.02, 1); 
% plot SB14 w phase by file number for 10 w points around w center
E_start = find(abs(E - 5.3) < 0.02, 1); 
E_stop  = find(abs(E - 6) < 0.02, 1); 

sideband_list = [12 14 16 18]*1240/wavelength-IP(1); % select 16th harmonic of lowest ionization state 
slope_shift_array = zeros(size(-3:1:3)); 
offset_shift_array = slope_shift_array; 

figure; hold on; 
ColorSet = colormap(hsv); 
% set(gca, 'ColorOrder', ColorSet(1:15:end,:)); 
ColorSet = ColorSet(1:18:end,:); 
for ii=-3:1:3
    tmp = 0; 
    for jj = 1:numel(sideband_list)
        E_start = find(abs(E - sideband_list(jj)+0.4) < 0.02, 1); 
        E_stop  = find(abs(E - sideband_list(jj)-0.4) < 0.02, 1); 
        tmp = tmp + squeeze(sum(E_SpectraArray_FFT(E_start:E_stop,w_center_ind+ii,:),1)); 
    end
    plot(unwrap(angle(tmp)) - unwrap(angle(tmp(1))), 'o-', 'color', ColorSet(ii+4,:)); 
    ptmp = polyfit(1:size(E_SpectraArray_FFT,3), unwrap(angle(tmp))'- unwrap(angle(tmp(1))), 1); 
    slope_shift_array(ii+4) = ptmp(1); 
    offset_shift_array(ii+4) = ptmp(2); 
    plot(1:1:14, offset_shift_array(ii+4) + slope_shift_array(ii+4)*(1:1:14), '-', 'color', ColorSet(ii+4,:)); 
end
xlabel('file number'); 
ylabel('phase'); 
title('phase of different w components for SB14'); 

%% 
twoOmega_alias = squeeze(zeros([numel(E), size(E_SpectraArray_FFT,2)])); 
for ii=1:size(E_SpectraArray_FFT,2)
    for jj = 1:14
        twoOmega_alias(:,ii+4) = twoOmega_alias(:,ii+4) + ...
            E_SpectraArray_FFT(:, w_center_ind + ii, jj) .* ...
            exp(-1j*(slope_shift_array(ii+4).*jj)); 
    end
end

% E_start = find(abs(E - 5.3) < 0.02, 1); 
% E_stop  = find(abs(E - 6) < 0.02, 1); 
E_start = find(abs(E - 8.2) < 0.02, 1); 
E_stop  = find(abs(E - 9) < 0.02, 1); 
figure; hold on; 
yyaxis left; 
plot(fftaxis(w_center_ind-3:w_center_ind+3), abs(sum(twoOmega_alias(E_start:E_stop,:),1))); 
yyaxis right; 
plot(fftaxis(w_center_ind-3:w_center_ind+3), unwrap(angle(sum(twoOmega_alias(E_start:E_stop,:),1)))); 















