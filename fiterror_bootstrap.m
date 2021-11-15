%%%% bootstrap method for getting errorbars on fit
%%%% case resampling 

% Download the data files in the folder 2020_08_25-16scan from our Google
% Drive. 
folderName = 'myfoldername'; % location of the data you just downloaded


alternate = [1 1]; % this is used for older data files where we did alternating scans of reference and target gas. 
% It would tell which position to start in (alternate(1)) and how many to skip between (alternate(2)). We don't do this anymore. 
wavelength=810; 
% IP and IP_label are globals so I don't have to code them in as arguments
% for all of my dependent functions, particularly the ones that plot the
% expected locations of sidebands. This was not necessarily a good
% engineering decision, and we might want to change it in the future. 
global IP; IP = [15.7596];
global IP_label; IP_label = ["Ar 2P 3/2"]; 
% If you want to look at the two argon spin-orbit states: 
% global IP; IP = [15.763, 15.763+0.17749];
% global IP_label; IP_label = ["Ar 2P 3/2", "Ar 2P 1/2"]; % start with 2

[HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderName, 1, wavelength);  
HistTot_array = HistTot_array(:,:,alternate(1):alternate(2):end); 
XUV_only = sum(XUV_only(:,alternate(1):alternate(2):end),2); 
XUV_only_raw = XUV_only; 
% Use this plot to find the energy calibration peaks, usually in argon
figure; plot(sum(sum(HistTot_array,2),3)); 
%% load calibration
% If you have already done energy calibration for this experiment, it will
% be saved in the same folder. You can load it in  here, or redo it in the
% next cell. I typically do energy calibration in argon and then use that
% calibration for the adjacent molecular measurement. 

% load(strcat(folderName, 'calibration/Ar_calibration.mat')); 
%% OR redo it
t0=25; % check this value by looking for the "light peak" near time zero. It can shift around a little. 
wavelength = 810; 
n = 10:1:21; % these are the peaks you're going to use for your energy calibration. I included sidebands here. 
calibEnergy = repmat(n, [size(IP,2) size(IP,1)])*1240/wavelength - repmat(IP, [size(n,2), size(n,1)]).'; 
tof_peak = [2107 1354 1077 929 821 750 695 653 614 581 554]; 

% config parameters for ECalibrate function. Can't remember why I did it
% this way. Seems not useful compared to simpler input structure. 
config.calibEnergy = calibEnergy; 
config.tofPeaks = tof_peak;   % redo with peak finding
config.IPcal = IP; 
config.Plot = 1; 

% Find energy calibration to be used in Convert_Eng_V2.m
A = ECalibrate(t0, n, wavelength, config); 
filename = strcat(folderName, 'calibration/', '_calibration'); 
save(filename, 'A', 'calibEnergy', 'tof_peak', 'wavelength'); 
%% Fourier Analysis begins here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use to check out your experiment file and get preliminary phase values. 
% Run before doing full analysis on the original file before bootstrapping. This will give
% structures to reference for size in the loop. 

% do energy conversion
tmp = reshape(HistTot_array, size(HistTot_array,1), []);
tof = 1:size(tmp,1); 
E_vec = [0 25 900]; % range 0-25eV with 900 bins. Adjust as you want. Increasing bins can cause artifacts. 

% convert the ToF spectrum to energy (linear energy scale)
% Convert_Eng_V2.m rebins tof into energy by calculating the overlap matrix
% (OM)
[C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be reused in bootstrap loop, which is good because it takes a while to calculate it
E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
Ebins = E; % save for later in case E gets overwritten by some of my other code for theoretical Wigner delay

% normalize on each delay step (dim 2 in array). sometimes I do this,
% sometimes I don't. 
norm = sum(E_SpectraArray,1); 
% E_SpectraArray = E_SpectraArray./repmat(norm, size(E_SpectraArray, 1), 1, 1); 

%% At any point after this, if you want to look at your raw time averaged spectrum you can run
% plotfun_rabbitspectrum(n, wavelength, E, data, mode)
% n=peaks you want to see, so probably something like 11:1:21 or 12:2:20
% wavelength = 810
% E=Ebins
% data = sum(sum(abs(tmp),2),3)
% mode = 'average'
%% fft and filter

% This step is kind of annoying and I never automated it. 
% You want to run this cell and then find the 2w bin on the plot using the
% cursor. Just the positive omega side index. 
% Then rewrite the value for twoOmega_ind using what you just found. 
% Then run this cell again with the correct index value. 

tmp = fftshift(fft(ifftshift(E_SpectraArray - repmat(mean(E_SpectraArray,2),[1 numel(stageTimes)]),2),[],2),2); 

figure; hold on; 
imagesc(sum(abs(tmp),3))

twoOmega_ind = 130; % MAKE THIS AUTOMATICALLY DETECTED EVENTUALLY? 
twoOmega_signal = squeeze(sum(tmp(:,twoOmega_ind,:),2)); 

twoOmega_nosum = twoOmega_signal; % save for later use

%% 2w analysis drift comensation

twoOmega_signal = twoOmega_nosum; % reset value in case you want to run this cell multiple times
sideband_list = [18]*1240/wavelength-IP(1); % select 18th harmonic of lowest ionization state. At SB18, the electron has high kinetic energy so we can use this as a global zero phase. Actually now that I think about it maybe we want to choose a higher ionization state for our molecular targets. 

peak_phase = 0; 
for ii=1:1:length(sideband_list) % sometimes we average over multiple sidebands (which is why 18 is in brackets above). I haven't found this to make much of a difference. 
    [~, index] = min(abs(E - sideband_list(ii)));
    window_center = index; 
    window = 2; % for argon you can make this window bigger, maybe 4
    histogram_windows = twoOmega_signal((window_center-window):(window_center+window),:); 
    % integrate over window
    peak_phase = peak_phase + squeeze(sum(histogram_windows, 1)); 
end
peak_amp = abs(peak_phase); 
peak_phase = angle(peak_phase);

phase_shift = peak_phase; % too lazy to go around changing variable names, sorry

for ii=1:1:length(phase_shift)
   twoOmega_signal(:,ii) = twoOmega_signal(:,ii).*exp(-1j*phase_shift(ii)); % apply global phase offsets to twoOmega signal
end
% sum phase matched values
twoOmega_signal = squeeze(sum(twoOmega_signal, 2)); % Add together phase corrected runs

%% At any point after this, if you want to look at your 2w spectrum you can run 
% plotfun_rabbitspectrum(n, wavelength, E, data, mode)
% n=peaks you want to see, so probably something like 11:1:21 or 12:2:20
% wavelength = 810
% E=Ebins
% data = twoOmega_signal
% mode = 'twoOmega'
%% convert XUV only data
% run this if you need XUV only data for something, like making plots

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

%% CHECK FITTING
% This is the script for doing spectral congestion fitting
% If peaks are very overlapped, we can try to extract separate phase values
% for each peak by fitting a sum of complex phases weighted by mixed
% gaussian amplitudes. 
% We have to limit the number of peaks in this kind of fit, so I select
% regions with 2-4 overlapping peaks. Adjust regions as needed to get the
% best fit. 
% I like to use this even if there is only a single peak so that I am doing
% consistent analysis across experiments, particularly for comparing
% reference to target. 
% Also, argon does have two visible spin-orbit states that we can
% disentangle using this method. 

signal = twoOmega_signal; 

% Ar fitting regions
region = [2.44 2.8]; % sideband 12
% region = [5.4 5.85]; % sideband 14
% region = [8.2 9.2]; % sideband 16
% region = [11.3 12.4]; % sideband 18

tolerance = abs(E(2)-E(1)); 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

% complexfit_section_full.m pulls on some supporting functions to extract
% phase and amplitude parameters for any peaks in the selected energy
% region, determined by IP in the first cell. 
% the parameters are: paramout = [phase slope, phase offset],
% paramout_gauss = [amp, center, sigma]. Each row is a different peak,
% ordered by IP. 
[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), abs(signal(start:stop)), signal(start:stop), 1, 1); 
% save "originals" to use as initial guess in boostrapping
paramout_original = paramout;  
fval_original = fval

%% CASE RESAMPLING (BOOTSTRAP)
% Run this after running the above fitting test cell to do boostrapping
% analysis on a selected region. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trials = 100; 
numfiles = numel(phase_shift); 
% array to save fit parameters
paramout_array = cat(2, zeros([size(paramout_gauss), trials]), zeros([size(paramout_gauss,1), 2, trials])); 
fval_array = zeros([1, trials]); 
% at some point I started saving all of the boostrap generated arrays to
% make sure nothing funny is happening. 
sample_2w_array = zeros([size(twoOmega_signal,1), trials]); 
sample_ES_array = zeros([size(twoOmega_signal,1), trials]); 

for nn=1:trials
    % resample files with replacement
    sample_list = randsample(numfiles, numfiles, true); 
    sample_ESpectra = E_SpectraArray(:,:,sample_list); 
    tmp = fftshift(fft(sample_ESpectra, [], 2), 2); 
    sample_twoOmega = squeeze(sum(tmp(:,twoOmega_ind,:),2));
    for mm=1:1:numfiles
        sample_twoOmega(:,mm) = sample_twoOmega(:,mm).*exp(-1j*phase_shift(sample_list(mm))); 
    end
    sample_twoOmega = sum(sample_twoOmega, 2); 
    sample_2w_array(:,nn) = squeeze(sample_twoOmega); 
    sample_ES_array(:,nn) = squeeze(sum(sum(sample_ESpectra,2),3)); 

    % fit data
    % complexfit_section_bootstrap is slightly modified from
    % complexfit_section_full
    [paramout, fval] = complexfit_section_bootstrap(wavelength, E(start:stop), sample_twoOmega(start:stop), paramout_gauss, paramout_original, 0); 
    paramout_array(:,:,nn) = paramout; 
    fval_array(nn) = fval; 

end 
check_if_done = 'done!'

%% check fit convergence on resampled sets
% Plots all the boostrap sampled sets, raw data and spectra made from extracted
% phases and mixed gaussian model. Make sure it converged and the sampled
% datasets are useful. 

xin = E(start:stop); 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    p = plot(xin, abs(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))), 'Color', 'c', 'LineWidth', 2); 
%     p = plot(xin, abs(yfit_array(:,ii)), 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
end
plot(xin, abs(Spectrum(xin, paramout_gauss, paramout_original)), 'Color', 'b', 'LineWidth', 2); 
% plot(xin, mean(abs(yfit_array),2), 'Color', 'b', 'LineWidth', 2); 
% scatter(xin, abs(twoOmega_signal(start:stop))/sum(abs(twoOmega_signal(start:stop))), 'bo'); 
title('Fit Amplitude', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Amplitude'); 
hold off; 

figure; hold on; grid on; 
for ii=1:size(paramout_array,3)
    plotphase = mod(angle(Spectrum(xin, paramout_array(:,1:3,ii),paramout_array(:,4:end,ii))),2*pi); 
%     plotphase = unwrap(mod(angle(yfit_array(:,ii)),2*pi)); 
%     plotphase = plotphase - plotphase(1); 
    p = plot(xin, plotphase, 'Color', 'c', 'LineWidth', 2); 
    p.Color(4) = 0.1; 
    s = plot(xin, mod(angle(sample_2w_array(start:stop,ii)),2*pi), 'cs');
    s.Color(4) = 0.1; 
end
plotphase = unwrap(mod(angle(Spectrum(xin, paramout_gauss, paramout_original)),2*pi)); 
% plotphase = unwrap(mod(mean(angle(yfit_array),2),2*pi)); 
plotphase(1)
% plotphase = plotphase - plotphase(1); 
plot(xin, plotphase, 'Color', 'b', 'LineWidth', 2); 
plot(xin, unwrap(mod(angle(Spectrum(xin, mean(paramout_array(:,1:3,:),3), mean(paramout_array(:,4:end,:),3))),2*pi)), 'b--', 'Linewidth', 2);
scatter(xin, unwrap(mod(angle(twoOmega_signal(start:stop)),2*pi)), 'bo'); 
title('Fit Phase', 'FontSize', 16); 
xlabel('Photoelectron Energy (eV)'); ylabel('Phase (referenced to v=6 phase)'); 
hold off; 

%% save 
% saves in the correct format to use in my H2_plot.m files and other H2
% analysis. You will probably do something smarter. 
data = paramout_array; 

bootstrap_phase = mean(mod(data(:,4,:),2*pi),3);
bootstrap_slope = mean(data(:,5,:),3); 

bootstrap_phase_std = sqrt(trials/(trials-1))*std(mod(data(:,4,:),2*pi),0,3);
bootstrap_slope_std = sqrt(trials/(trials-1))*std(data(:,5,:),0,3);

% Ar_SB18_paramout = data; 
% Ar_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
% Ar_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
% H2_SB18_paramout = data; 
% H2_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
% H2_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
Kr_SB18_paramout = data; 
Kr_SB18_phase = [bootstrap_phase, bootstrap_phase_std]; 
Kr_SB18_slope = [bootstrap_slope, bootstrap_slope_std]; 
% tmp = [bootstrap_phase, bootstrap_phase_std]; 

% bootstrap_d1 = mean(data(:,end,:),3);
% bootstrap_d2 = mean(data(:,end-1,:),3); 
% % bootstrap_d3 = mean(data(:,end-2,:),3); 
% 
% bootstrap_d1_std = sqrt(trials/(trials-1))*std(data(:,end,:),0,3);
% bootstrap_d2_std = sqrt(trials/(trials-1))*std(data(:,end-1,:),0,3);
% % bootstrap_d3_std = sqrt(trials/(trials-1))*std(data(:,end-2,:),0,3);
% 
% Ar_SB12_paramout = data; 
% Ar_SB12_phase = cat(3, [bootstrap_d1, bootstrap_d1_std], ...
%                        [bootstrap_d2, bootstrap_d2_std]); 
% %                        [bootstrap_d3, bootstrap_d3_std]); 
% % H2_SB14_paramout = data; 
% % H2_SB14_phase = cat(3, [bootstrap_d1, bootstrap_d1_std], ...
% %                        [bootstrap_d2, bootstrap_d2_std], ...
% %                        [bootstrap_d3, bootstrap_d3_std]); 


%% JACKKNIFE
% very similar to boostrapping but no replacement in sampling
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
    sample_twoOmega = squeeze(tmp(:,twoOmega_ind,:));
 
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
        % peak_phase = angle(peak_fft(twoOmega_ind, :)); 
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








