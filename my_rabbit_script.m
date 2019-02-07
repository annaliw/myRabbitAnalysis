% load data into arrays, don't do anything else in this cell 
clear all 
[HistTot_array, stageTimes, freqAxis] = getrawdata('/Users/annaliw/code/KrCO2_scan/', 1); 

%%
t0 = 75; 
wavelength = 810; 
IP = [13.778, 17.706, 18.077, 19.394]; 
IP_label = ["X", "A", "B", "C"]; 
calibType = 'Kr'; 
E_vec = [0 25 500]; % need to add no binning option to get2wdata
alternate = [2, 2]; % right now there are two optional variables in get2wdata. 
                    % if plotting is used but not alternate it will not
                    % work. Need both or neither. Fix this. 
plotting = 0; 


[E, E_SpectraArray, twoOmega_signal] = get2wdata(HistTot_array, t0, IP, calibType, E_vec, alternate, plotting); 

plotfun_rabbitspectrum(9:1:19, IP, IP_label, 810, E, twoOmega_signal, 'twoOmega');

%% redo wavelength calculation for this data (maybe different for Kr calibration set)
Xpeaks = [14, 15, 16, 17, 18, 19; 11.7225, 13.2775, 14.7869, 16.2276, 17.8971, 19.3149]; 
bpeaks = [14, 15, 16, 17, 18, 19; 4.9667, 6.4476, 7.9929, 9.4416, 10.9547, 12.5322]; 
% [harmonic_lin,~,mu] = polyfit(Xpeaks(1,:), Xpeaks(2,:), 1); 
fitx = [Xpeaks(1,:), bpeaks(1,:)]; 
fity = [Xpeaks(2,:)+IP(1), bpeaks(2,:)+IP(2)]; 
harmonic_lin = polyfit(fitx, fity, 1); 

figure; hold on; 
scatter(fitx, fity); 
plot(fitx, polyval(harmonic_lin, fitx));  
hold off;

wavelength_mod = 1240/harmonic_lin(1); 
nshift = (9:1:19)+harmonic_lin(2)/harmonic_lin(1); 

plotfun_rabbitspectrum(nshift, IP, IP_label, wavelength_mod, E, E_SpectraArray, 'twoOmega');

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


%%

signal = [twoOmega_abs_linear, twoOmega_abs_quadratic, twoOmega_linear_ex, twoOmega_quadratic_ex]; 
IP = [9.553 16.56 18.319 21.722]; 

fh = figure; 
tmp = axes; 
xlabel('Electron Energy (eV)')

axl = AddHarmonicAxis(fh,IP, wavelength);

axl(1).XLabel.String = 'X';
axl(2).XLabel.String = 'b ^3\Pi';
axl(3).XLabel.String = 'A^1\Pi';
axl(4).XLabel.String = 'B^1\Pi or c^3\Pi'; 

for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.2903    0.99    0.0000];
end
hold on; 
for i=1:1:length(signal(1,:)) 
    plot(E, abs(signal(:,i)));
end
legend('linear fit', 'quadratic fit', 'linear fit excluding point', 'quadratic fit excluding point'); 

xlim([0 25]); 
% delete(tmp); 

hold off;  

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
E_SpectraArray = zeros(size(HistTot_array)); 
for i=1:1:223
    E_SpectraArray(:,i) = -1/(2*A(2)) .* t.^3 .* HistTot_array(:,i).';  
end





