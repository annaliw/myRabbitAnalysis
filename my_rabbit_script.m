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

%% Convert to Energy after loading 
% E_vec_max = 25; 
% dE = 2*(19*(1240/810)-9.3)/510; 
% E_vec_size = floor(E_vec_max/dE); 
% % E_vec_size = 1200; 
% E_vec = [0, E_vec_max, E_vec_size];

t0=21; 

%get conversion parameters
% calibEnergy = [((11:1:19)*(1239.84e-9/IRwavelength) - 14.00) ((11:1:19)*(1239.84e-9/IRwavelength) - 14.665)];
calibEnergy = [((11:1:19)*(1240/810) - 14.00) ((11:1:19)*(1240/810) - 14.665)];
n = 9:2:19; 
% n = 11:1:19; 
% tof_peak = [1082, 793, 687, 610, 553, 510]-t0; 
tof_peak = [1256 1030 892 804 734 684 639 604 573 1424 1118 946 838 763 707 657 619 585]-t0;
tof_fitvar = 1./tof_peak.^2; 
IP = 9.3; 
% En_X = n*(1240/810)-IP; 
% En_X = [calibEnergy n*(1240/810)-IP]; 
En_X = calibEnergy; 

A = polyfit(tof_fitvar, En_X, 1); 
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

%% Try plotting things? 
E_Spectra_dim = size(E_SpectraArray); 
num_spectra = E_Spectra_dim(2); 
E_Spectra_average = sum(abs(E_SpectraArray),2)/num_spectra; 

figure; hold on; 
plot(E, abs(E_Spectra_average)); 
xlabel('Energy (eV)'); 
ylabel('Average counts'); 
hold off; 

%% FFT plotting
IRwavelength = 810e-09; 

%figure out time and frequency axis for upcoming FFT
dt = (stageTimes(end) - stageTimes(1))/(length(stageTimes));
TimeWindow = stageTimes(end) - stageTimes(1);
sampleCount = numel(stageTimes);
FreqWindow = 1/dt;
df = 1/TimeWindow;
freqScaling = (2.9979e8/IRwavelength);
freqAxis = (-sampleCount/2:1:(sampleCount/2 - 1)) * df / freqScaling;

% tx = (stage_positions - mean(stage_positions(1:end))) * (1E3*2/0.29979); %Time Delay Label (in fs), assuming t0 is centered in the scan range. 
% %the  stage_positions(1:end-1) is because we usually take a scan where the last point is very far off from overlap to get the XUV only spectrum for reference. 
% T = tx(end) - tx(1); %Total time duration
% Np = numel(tx); % Number of delay points (excluding the early time)
% df = 1/T; %The frequency step in the fft is 1 over the total time duration
% f = ((-Np/2):1:(Np/2-1)) * df;
% f = f * (0.810/0.29979); % This changes the frequency label to multiples of the fundamental laser frequency. 
% 
% %perform the FFT
% E_Spectra_freq = fft(abs(E_SpectraArray - mean(E_SpectraArray, 2)), [], 2);
% E_Spectra_freqShifted = fftshift(E_Spectra_freq, 2);
    
figure; hold on; 
surf(freqAxis, E, abs(E_SpectraArray), 'LineStyle', 'None'); 
colorbar; 
hold off; 

%% grab 2w data and Plot Result with Harmonics
twoOmega_signal = E_SpectraArray(:,130); 
% IP = [(9.262+9.553+9.839)/3,16.56,18.319,21.722];
IP = [(9.262+9.553+9.839)/3,16.56,18.319,22.7];

fh = figure; 
line(E, mean(abs(E_SpectraArray),2), 'Color', 'k', 'DisplayName', 'average spectra'); 
% line(E, mean(abs(E_SpectraArray_early),2), 'Color', 'r', 'DisplayName', 'early data average'); 
% line(E, sum(abs(E_SpectraArray),2)/num_spectra - sum(abs(E_SpectraArray_early),2)/45, 'Color', 'k', 'DisplayName', 'difference spectra'); 
xlabel('Electron Energy (eV)')
% make sideband lines... 
m = 10:2:18; 
for i=1:1:length(IP)
    sideband = m*(1240/810)-IP(i);
    top = max(mean(abs(E_SpectraArray),2)); 
    bottom = 0; 
    for j=1:1:length(sideband)
        line([sideband(j), sideband(j)], [top, bottom], 'LineStyle', '--'); 
    end
end
% sideband = m*(1240/810)-IP(end);
%     top = 2.5e05; 
%     bottom = -5e04; 
%     for j=1:1:length(sideband)
%         line([sideband(j), sideband(j)], [top, bottom], 'LineStyle', '--'); 
%     end
axl = AddHarmonicAxis(fh,IP);

axl(1).XLabel.String = 'X (HOMO, average v=0-2)';
axl(2).XLabel.String = 'b ^3\Pi ';
axl(3).XLabel.String = 'A^1 \Pi (v=0)';
axl(4).XLabel.String = 'B ^1\Pi';

for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.2903    0.99    0.0000];
end

% yyaxis left
twoOmega_abs = abs(twoOmega_signal); 
twoOmega_phi = angle(twoOmega_signal); 
hold on; 
plot(E, twoOmega_abs, 'b-', 'DisplayName', '2w amplitude');
yyaxis right
plot(E, unwrap(twoOmega_phi), '-', 'DisplayName', '2w phase'); 
legend; 

xlim([0 25]); 

hold off;  





