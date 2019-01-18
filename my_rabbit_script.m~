% load data into arrays, don't do anything else in this cell 
clear all 

% time zero (take as parameter later)
t0 = 21; 

% data location
folderString = '/Users/annaliw/code/NOscan/'; 
saveString = 'extracted_data'; 

% get data files
dataList = dir(fullfile(folderString, '*.mat'));
numberSubScans = length(dataList);

% % load a single scan for reference
% currentSubScan = dataList(1); 
% load(char(string(currentSubScan.folder) + '/' + string(currentSubScan.name))); 
% % get delay averaged histogram (just a sum of all positions)
% histSum = sum(HistTot,2)/length(stage_positions); 
% figure; hold on; 
% plot(histSum); plot(HistTot(:,end)); 
% hold off; 
% % I will track a sideband that appears at tof bin 531 in the reference scan

% load in all the files 
for ii=1:1:numberSubScans
    creationCheck = exist('HistTot'); 
    if creationCheck == 0 
        makeArrays = 1; 
    else
        makeArrays = 0; 
    end
    currentSubScan = dataList(ii); 
    load(char(string(currentSubScan.folder) + '/' + string(currentSubScan.name))); 
    %SigSum is sum of anode signal
    %HistTot is the set of histograms
    %stage_positions are the stage positions for each scan
    if makeArrays == 1  
        SigSum_array = zeros([size(SigSum), numberSubScans]); 
        HistTot_array = zeros([size(HistTot), numberSubScans]); 
        stagePositions_array = zeros([size(stage_positions), numberSubScans]); 
    end
    % saving each subscan data in an array
    SigSum_array(:,:,ii) = SigSum; 
    HistTot_array(:,:,ii) = HistTot; 
    stagePositions_array(:,:,ii) = stage_positions; 
end

%% make stage times
stage_times = (stage_positions*2*1e-3)/(2.9979e8);
stageCenter = stage_times(round(length(stage_times)/2));
stageTimes_array = (stagePositions_array*2*1e-3)/(2.9979e8) - stageCenter;

%% trim off non-overlapped data 
HistTot_early = HistTot_array(:,1:45,:); 
stageTimes_early = squeeze(stageTimes_array(:,1:45,1));

HistTot_array = HistTot_array(:,46:end-1,:); 
stageTimes = squeeze(stageTimes_array(:,46:end-1,1));
% (or not...)
% histTotSum = sum(HistTot_array, 3);
% stageTimes = squeeze(stageTimes); 


%% Convert to Energy after loading (no binning)
t0=21; 
wavelength = 810; 
%get conversion parameters
% calibEnergy = [((11:1:19)*(1239.84e-9/IRwavelength) - 14.00) ((11:1:19)*(1239.84e-9/IRwavelength) - 14.665)];
calibEnergy = [((11:1:19)*(1240/wavelength) - 14.00) ((11:1:19)*(1240/wavelength) - 14.665)];
n = 9:2:19; 
% n = 11:1:19; 
% tof_peak = [1082, 793, 687, 610, 553, 510]-t0; 
tof_peak = [1256 1030 892 804 734 684 639 604 573 1424 1118 946 838 763 707 657 619 585]-t0;
tof_fitvar = 1./tof_peak.^2; 
IP = 9.3; 

A = polyfit(tof_fitvar, calibEnergy, 1); 
xFit = (1:1:1000)*(tof_fitvar(end)-tof_fitvar(1))/1000 + tof_fitvar(1); 
yFit = polyval(A, xFit);

A = fliplr(A); 

% %convert the ToF spectrum to energy (linear energy scale)
% [C,E]=Convert_Eng_V2(1:length(HistTotSum(:,1)),HistTotSum,[t0, A] , E_vec);
% %store the converted Energy spectrum
% E_SpectraArray = C'; 

% Convert tof -> energy no rebinning (nonlinear scale)
tof = 1:1:8192; 
t = (tof-t0); %Adjust for prompt time and convert to ns
x = 1./t.^2;

E = 0;
for i = 0:numel(A)-1
    E = E + A(i+1) .* x.^i;
end
E(tof<=t0)=0; 

% compensate for change in variables 
E_SpectraArray = zeros(size(HistTotSum)); 
E_SpectraArray_early = zeros(size(HistTotSum_early)); 
for i=1:1:223
    E_SpectraArray(:,i) = -1/(2*A(2)) .* t.^3 .* HistTotSum(:,i).';  
end
for i=1:1:45
   E_SpectraArray_early(:,i) = -1/(2*A(2)) .* t.^3 .* HistTotSum_early(:,i).'; 
end
% E_SpectraArray = HistTotSum; 

% cut out zero energy region (save old...)
E_full = E; 
E_SpectraArray_full = E_SpectraArray; 
E_SpectraArray_early_full = E_SpectraArray_early; 
E = E(451:4748); 
E_SpectraArray = E_SpectraArray(451:4748,:); 
E_SpectraArray_early = E_SpectraArray_early(451:4748,:); 

%% Energy conversion (binning)
wavelength = 810; 
E_vec_max = 25; 
dE = 2*(19*(1240/wavelength)-9.3)/510; 
% E_vec_size = floor(E_vec_max/dE); 
E_vec_size = 400; 
E_vec = [0, E_vec_max, E_vec_size];

%get conversion parameters
% calibEnergy = [((11:1:19)*(1239.84e-9/IRwavelength) - 14.00) ((11:1:19)*(1239.84e-9/IRwavelength) - 14.665)];
calibEnergy = [((11:1:19)*(1240/wavelength) - 14.00) ((11:1:19)*(1240/wavelength) - 14.665)];
n = 9:2:19; 
% n = 11:1:19; 
% tof_peak = [1082, 793, 687, 610, 553, 510]-t0; 
tof_peak = [1256 1030 892 804 734 684 639 604 573 1424 1118 946 838 763 707 657 619 585]-t0;
tof_fitvar = 1./tof_peak.^2; 
IP = 9.3; 

A = polyfit(tof_fitvar, calibEnergy, 1); 
xFit = (1:1:1000)*(tof_fitvar(end)-tof_fitvar(1))/1000 + tof_fitvar(1); 
yFit = polyval(A, xFit);

A = fliplr(A); 

tobeconverted = HistTotSum; 

%convert the ToF spectrum to energy (linear energy scale)
[C,E]=Convert_Eng_20190116(1:length(tobeconverted(:,1)),tobeconverted,[t0, A] , E_vec);
% [Cf,Ef]=Convert_Eng_V2(1:length(HistTotSum_full(:,1)),HistTotSum_full,[t0, A] , E_vec);
% [Ce, Ee] = Convert_Eng_V2(1:length(HistTotSum_early(:,2)), HistTotSum_early, [t0, A], E_vec); 
%store the converted Energy spectrum
E_SpectraArray = C.'; 
% E_SpectraArray_early = Ce.'; 

% % cut out zero energy region (save old...)
% E_SpectraArray_early_full = E_SpectraArray_early; 
% E = E(451:4748); 
% E_SpectraArray = E_SpectraArray(451:4748,:); 
% E_SpectraArray_early = E_SpectraArray_early(451:4748,:); 

figure;
plot(E,abs(C(112,:)))

clear('tobeconverted', 'C', 'Ce', 'Ee'); 


%% identify and compensate for stage drift
% will compare location of peak (at tof bin 531 in first subscan) 

% cut out window of ~20 tof bins
window_center = 531; 
window = 8; 
histogram_windows = HistTot_array(window_center-window:window_center+window,:,:); 
% integrate over window
peak_vol = squeeze(sum(histogram_windows, 1)); 
% fft wrt IR delay 
twoOmega_location = 130; % from frequency axis calculated previously
peak_fft = fftshift(fft(peak_vol, [], 1), 1); 
peak_phase = angle(peak_fft(twoOmega_location, :)); 

% % reference everything to first scan. subtract off drift. 
% peak_phase_offset = peak_phase - peak_phase(1); 
% peak_phase_offset = reshape(peak_phase_offset, [1, size(peak_phase_offset)]); 
% % convert phase to time
% peak_time_offset = peak_phase_offset/(2*2.9979e8/810e-09); 
% stageTimes_array = stageTimes_array - peak_time_offset; 

HistTot_fft = fftshift(fft(HistTot_array, [], 2), 2); 
% HistTot_phase = squeeze(angle(HistTot_fft(:,130,:))); 
HistTot_early_fft = fftshift(fft(HistTot_early, [], 2), 2); 
% HistTot_early_phase = squeeze(angle(HistTot_early_fft(:,130,:))); 

for ii=1:1:length(peak_phase)
   HistTot_fft(:,:,ii) = HistTot_fft(:,:,ii).*exp(-1i*peak_phase(ii)); 
   HistTot_early_fft(:,:,ii) = HistTot_early_fft(:,:,ii).*exp(-1i*peak_phase(ii)); 
end
% sum phase matched values
HistTotSum = sum(HistTot_fft, 3); 
HistTotSum_early = sum(HistTot_early_fft, 3); 

clear('window', 'window_center','histogram_windows','histogram_windows','peak_vol','peak_fft', 'peak_phase'); 

%% plot Kr Spectrum with harmonics
fh = figure; 
% stairs(E,C);
line(E,sum(abs(E_SpectraArray),2)/224); 
xlabel('Electron Energy (eV)')

IP = [14.0, 14.665];

axl = AddHarmonicAxis(fh,IP);

axl(1).XLabel.String = 'Kr state 1';
axl(2).XLabel.String = 'Kr state 2';

for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.2903    0.99    0.0000];
end


%% grab 2w data and Plot Result with Harmonics
oneOmega_signal = E_SpectraArray(:,121); 
twoOmega_signal = E_SpectraArray(:,130); 
IP = [(9.553+9.839+10.121)/3,16.56,18.319,21.722];
% IP = [(9.262+9.553+9.839)/3,16.56,18.319,21.7];
% IP = [(9.553+9.839+10.121)/3, 16.56, 15.667, 15.8, 15.9, 16.11, 16.26]; %, 18.318, 21.722];%10.39,10.67]; %

fh = figure; 
line(E, mean(abs(E_SpectraArray),2), 'Color', 'k', 'DisplayName', 'average spectra'); 
% line(E, mean(abs(E_SpectraArray_early),2), 'Color', 'r', 'DisplayName', 'early data average'); 
% line(E, sum(abs(E_SpectraArray),2)/num_spectra - sum(abs(E_SpectraArray_early),2)/45, 'Color', 'k', 'DisplayName', 'difference spectra'); 
xlabel('Electron Energy (eV)')

axl = AddHarmonicAxis(fh,IP, wavelength);

axl(1).XLabel.String = 'X (avg v. 2-4)';
axl(2).XLabel.String = 'b ^3\Pi';
% axl(3).XLabel.String = 'A^1\Pi';
% axl(4).XLabel.String = 'B^1\Pi or c^3\Pi'; 

for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.2903    0.99    0.0000];
end

% yyaxis left
twoOmega_abs = abs(twoOmega_signal); 
twoOmega_phi = angle(twoOmega_signal); 
% oneOmega_abs = abs(oneOmega_signal); 
% oneOmega_phi = angle(oneOmega_signal); 
hold on; 
plot(E, movmean(twoOmega_abs./(mean(abs(E_SpectraArray), 2)), 3), 'b-', 'DisplayName', '2w amplitude');
% plot(E, movmean(oneOmega_abs, 3), 'r-', 'DisplayName', '1w amplitude'); 
yyaxis right
plot(E, twoOmega_phi, 'c-', 'DisplayName', '2w phase'); 
% plot(E, unwrap(oneOmega_phi), 'm-', 'DisplayName', '1w phase'); 
legend; 

xlim([0 25]); 
title('400 bins energy conversion'); 

hold off;  

%% check mystery peaks

peak1 = sum(E_SpectraArray(find(E-3.765<0.005):find(E-3.582<0.005),:),1); 
peak2 = sum(E_SpectraArray(find(E-4.05<0.005):find(E-3.87<0.005),:),1); 
peak3 = sum(E_SpectraArray(find(E-4.368<0.005):find(E-4.175<0.005),:),1); 
peak4 = sum(E_SpectraArray(find(E-1.284<0.005):find(E-1.173<0.005),:),1);

peak1_s = sum(E_SpectraArray(find(E-3.765-1.53<0.005):find(E-3.582-1.53<0.005),:),1);
peak2_s = sum(E_SpectraArray(find(E-4.05-1.53<0.005):find(E-3.87-1.53<0.005),:),1); 
peak3_s = sum(E_SpectraArray(find(E-4.368-1.53<0.005):find(E-4.175-1.53<0.005),:),1); 

peaks = [peak1; peak2; peak3; peak4; peak1_s; peak2_s; peak3_s]; 

peaks_ifft = ifft(ifftshift(peaks, 2), [], 2); 

noIR_peak = sum(E_SpectraArray_early(find(E-1.284<0.005):find(E-1.173<0.005),:),1);
noIR_ifft = ifft(ifftshift(noIR_peak,2),[],2); 

figure; hold on; 
plot(abs(peaks_ifft(1,:)), 'DisplayName', 'peak1');
plot(abs(peaks_ifft(2,:)), 'DisplayName', 'peak2');
plot(abs(peaks_ifft(3,:)), 'DisplayName', 'peak3');
plot(abs(peaks_ifft(4,:)), 'DisplayName', 'peak4');
plot(abs(noIR_ifft), 'DisplayName', 'noIR'); 
legend; 
% figure; hold on; 
% title('potential sidebands'); 
% plot(abs(peaks_ifft(4,:)), 'DisplayName', 'peak1'); 
% plot(abs(peaks_ifft(5,:)), 'DisplayName', 'peak1'); 
% plot(abs(peaks_ifft(6,:)), 'DisplayName', 'peak1'); 
% legend; 

%% early vs. overlap

IP = [(9.262+9.553+9.839)/3,16.56,18.319,21.7];

fh = figure; hold on;
% plot(E, sum(abs(E_SpectraArray),2)/sum(abs(E_SpectraArray(:))), 'k-', 'DisplayName', 'overlap'); 
% plot(E, sum(abs(E_SpectraArray_early),2)/sum(abs(E_SpectraArray_early(:))), 'r-', 'DisplayName', 'not overlap'); 
plot(E, sum(abs(E_SpectraArray),2)/sum(abs(E_SpectraArray(:))) - sum(abs(E_SpectraArray_early),2)/sum(abs(E_SpectraArray_early(:))), 'k-', 'DisplayName', 'XUV+IR - IR late'); 
% plot(E, sum(abs(E_SpectraArray),2)/sum(abs(E_SpectraArray(:))) - movmean(abs(E_SpectraArray(:,end))/sum(abs(E_SpectraArray(:,end))),3), 'r-', 'DisplayName', 'XUV+IR - IR early'); 
xlabel('Electron Energy (eV)')

axl = AddHarmonicAxis(fh,IP,wavelength);

axl(1).XLabel.String = 'X (HOMO, average v=0-2)';
axl(2).XLabel.String = 'b ^3\Pi ';
axl(3).XLabel.String = 'A^1 \Pi (v=0)';
axl(4).XLabel.String = 'B ^1\Pi';

for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.2903    0.99    0.0000];
end
 
hold off; 

%% plotting in fourier space

df = 1/(stageTimes(end)-stageTimes(1)); 
freqAxis = (-(numel(stageTimes)-1)/2:1:(numel(stageTimes)-1)/2)*df*(810*10^(-9))/(3*10^8); 




