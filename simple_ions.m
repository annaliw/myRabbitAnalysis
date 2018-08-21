% Simple load data does not account for stage drift and just adds together
% the histograms. 
clear all 

% time zero (take as parameter later)
t0 = 25; 

% data location
folderString = '/Users/annaliw/code/2018_07_19_scan/'; 
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

%% get stage times from stage positions
%convert stage positions to time
stage_times = (stage_positions*2*1e-3)/(2.9979e8);
stageCenter = stage_times(round(length(stage_times)/2));
stageTimes = stage_times - stageCenter;

%% Sum bins around peak 
peakRange = 100; 
peaks = [3202, 2556, 1933, 1675]; % CO2+, CO+, O+, C+ locations in ToF 
ionSum = zeros(length(peaks), length(stageTimes)); 

figure; hold on; 
for ii=1:1:length(peaks)
    ionSum(ii,:) = HistTotSum(peaks(ii), :); 

    plot(stageTimes(1:end-1), ionSum(ii,1:end-1)-mean(ionSum(ii,1:end-1))); 
end 

xlabel("IR delay"); 
ylabel("Ion Counts (-mean)"); 
legend("CO2+", "CO+", "O+", "C+"); 
hold off; 

%% FFT 
IRwavelength = 810e-09; 
ionSumIRfft = fftshift(fft(ionSum(:,1:end-1), [], 2), 2);

%figure out time and frequency axis for upcoming FFT
dt = (stageTimes(end-1) - stageTimes(1))/(length(stageTimes) - 2);
TimeWindow = stageTimes(end-1) - stageTimes(1);
sampleCount = numel(stageTimes)-1;

FreqWindow = 1/(2*dt);
df = 1/TimeWindow;

freqScaling = (2.9979e8/IRwavelength);
freqAxis = -FreqWindow:df:FreqWindow; 
freqAxis = freqAxis/freqScaling; 

figure; hold on; 
for ii=1:1:length(peaks)
    plot(freqAxis, abs(ionSumIRfft(ii,:))); 
end

xlabel("IR delay frequency"); 
ylabel("Ion Counts"); 
legend("CO2+", "CO+", "O+", "C+"); 
hold off;

figure; hold on; 
for ii=1:1:length(peaks)
    plot(freqAxis, unwrap(angle(ionSumIRfft(ii,:)))); 
end

xlabel("IR delay frequency"); 
ylabel("Phase Angle (of ion counts fft)"); 
legend("CO2+", "CO+", "O+", "C+"); 
hold off;

%% filter fft and transform back
filter = zeros(1, 223); 
filter(98:106) = 1; 
% filter(118:126) = 1; 
% filter(110:114) = 1; 

filteredFFTdata = filter.*ionSumIRfft; 
filteredData = ifft(ifftshift(filteredFFTdata, 2), [], 2); 

figure; hold on; 
plot(stage_times(1:end-1), abs(filteredData)-mean(abs(filteredData),2)); 
xlabel("time of flight")
ylabel("ion count (-mean)")
legend("CO2+", "CO+", "O+", "C+");
hold off; 
%% check charge/mass vs tof

mass = [12+16*2, 12+16, 16, 12]; 
chargetomass = 1./mass; 

xtest = 0:0.005:0.1; 
ytest = 1./sqrt(xtest); 

figure; hold on; 
plot(chargetomass, peaks); 
plot(xtest, ytest); 
xlabel("Charge to mass ratio"); 
ylabel("Time of flight (ns)"); 
hold off; 





