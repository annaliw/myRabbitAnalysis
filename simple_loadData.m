% Simple load data does not account for stage drift and just adds together
% the histograms. 
clear all 

% time zero (take as parameter later)
t0 = 27; 

% data location
folderString = '/Users/annaliw/code/2018_07_27-17Scan/'; 
saveString = 'extracted_data'; 

% get data files
dataList = dir(fullfile(folderString, '*.mat'));
numberSubScans = length(dataList);

% instead of partitioning the data I will sum all the histograms
HistTotSum = 0; 
for ii=1:1:numberSubScans
    currentSubScan = dataList(ii); 
    load(char(string(currentSubScan.folder) + '/' + string(currentSubScan.name))); 
    %SigSum is sum of anode signal
    %HistTot is the set of histograms
    %stage_positions are the stage positions for each scan
    if(ii == 1)
        HistTotSum = zeros(size(HistTot));
    end
    HistTotSum = HistTotSum + HistTot;
end

%% split into two data sets 
HistTotSum = HistTotSum(:, 46:end-1); 
stage_positions = stage_positions(46:end-1); 
% (or not...)
% histData = HistTotSum; 
%% Convert to Energy after loading
E_vec_max = 24; 
E_vec_size = 300; 
E_vec = [0.02, E_vec_max, E_vec_size];

E_SpectraArray=zeros(E_vec(3), size(HistTotSum,2));
% convert for each stage position (dimension 2 in histogram array)
for ii = 1:size(HistTotSum,2)
    %select at a ToF spectrum
    Y = HistTotSum(:,ii);
    %convert the ToF spectrum to energy
    [C,E]=Convert_Eng_V2(1:length(Y),Y,[t0],E_vec);
    %store the converted Energy spectrum
    E_SpectraArray(:,ii) = C;
end

%% Try plotting things? 
E_Spectra_dim = size(E_SpectraArray); 
num_spectra = E_Spectra_dim(2); 
E_Spectra_average = sum(E_SpectraArray,2)/num_spectra; 

figure; hold on; 
plot(E, E_Spectra_average); 
xlabel('Energy (eV)'); 
ylabel('Average counts'); 
hold off; 

%% FFT plotting
IRwavelength = 810e-09; 
%convert stage positions to time
stage_times = (stage_positions*2*1e-3)/(2.9979e8);
stageCenter = stage_times(round(length(stage_times)/2));
stageTimes = stage_times - stageCenter;

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

%perform the FFT
E_Spectra_freq = fft(abs(E_SpectraArray - mean(E_SpectraArray, 2)), [], 2);
E_Spectra_freqShifted = fftshift(E_Spectra_freq, 2);
    
figure; hold on; 
surf(freqAxis, E, abs(E_Spectra_freqShifted), 'LineStyle', 'None'); 
colorbar; 
hold off; 
