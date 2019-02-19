function [peak_nn, phase_nn, slope_nn] = fiterror_bootstrap(folderName, alternate, t0, wavelength, IP, E_vec, trials, fitRegion, peakFlag)
    %%
    if ~exist('peakflag', 'var')
        peakflag = 0; 
    end
    if length(fitRegion)==2
        slopeflag = 0; 
        slope_guess = 0; 
        width = 0.1; % default width guess value
    elseif length(fitRegion)==3
        slopeflag = 0;
        slope_guess = 0; 
        width = fitRegion(3); 
    elseif length(fitRegion)==4
        slopeflag = 1; 
        slope_guess = fitRegion(4); 
        width = fitRegion(3); 
    else
        error('invalid fitRegion size'); 
    end
    
    %%
    plotting = 1; % debug setting
    
    folderName = '/Users/annaliw/code/2018_07_31-16Scan/'; % debug setting
%     folderName = '/Users/annaliw/code/KrCO2_scan/'; 
    alternate = [1 2]; % debug setting
    t0=51; 
    wavelength=810; 
    IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206]; 
    IP_label = ["0", "1", "2", "3", "4", "5"]; 
%     IP = [13.9, 17.706, 18.077, 19.394]; 
%     IP_label = ["X", "A", "B", "C"]; 
    E_vec = [0 20 900]; 
    trials = 5; 
    energy_range_start = 1.5355; 
    energy_range_stop = 2.8258; 
    slope_guess = 1.0978; 
    fitRegion = [energy_range_start energy_range_stop 0.08 slope_guess]; 
    width = fitRegion(3); 
    slope_guess = fitRegion(4); 
    slopeflag = 1; % debug setting 
    peakflag = 0; % debig setting
    [HistTot_array, stageTimes, freqAxis] = getrawdata(folderName, 1, wavelength);  
    HistTot_array = HistTot_array(:,:,alternate(1):alternate(2):end); % might need to alternate files
%     % do data padding for better fitting
%     pad_data = zeros([size(HistTot_array,1)*10-9, size(HistTot_array,2), size(HistTot_array,3)]); % not sure how the sizing works here
%     for ii=1:1:size(HistTot_array, 2)
%         for jj=1:1:size(HistTot_array, 3)
%             pad_data(:,ii,jj) = interp1(0:1:(size(HistTot_array,1)-1), HistTot_array(:,ii,jj), 0:0.1:(size(HistTot_array,1)-1)); 
%             pad_data(:,ii,jj) = poissrnd(pad_data(:,ii,jj)); 
%         end
%     end
%     HistTot_array = pad_data; 
    
    %%
    % calibrate to Krypton
    calibType = 'Kr'; 
    % calibrate to Kr peaks
    A = ECalibrate(t0, calibType, 0); % TO DO: hard set t0 into ECalibrate
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % need to do full analysis on the original file. this will also give
    % structures to reference for size in the loop. 
    
    % do energy conversion
    tmp = reshape(HistTot_array, size(HistTot_array,1), []);
    tof = 1:size(tmp,1); 
    %convert the ToF spectrum to energy (linear energy scale)
    [C,E,OM]=Convert_Eng_V2(tof, tmp, [t0, A] , E_vec); % OM will be used in loop
    E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
    % fft data
    E_SpectraArray = fftshift(fft(E_SpectraArray, [], 2), 2); 
    twoOmega_location = 130; % MAKE THIS AUTOMATICALLY DETECTED
    twoOmega_signal = squeeze(E_SpectraArray(:,twoOmega_location,:));
    
    % only restack the 2w data
    % cut out window of ~6 tof bins
    % sideband_list = [(12:2:18)*1240/810 - 13.778, (12:2:18)*1240/810 - 17.706];
    sideband_list = [18*1240/wavelength-IP(1)]; % select 16th harmonic of lowest ionization state 
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
    for ii=1:1:length(peak_phase)
       twoOmega_signal(:,ii) = twoOmega_signal(:,ii).*exp(-1j*peak_phase(ii)); 
    end
    % sum phase matched values
    twoOmega_signal = sum(twoOmega_signal, 2); 

    %%
    % fit data
    fitRegion = [1.6387 2.8 0.08 slope_guess]; 
    fixguess = [1,4.5,1]; 
    [paramout, fval, exitflag, output] = ...
        complexfit_section(wavelength, E, twoOmega_signal, fitRegion, 'FixWidth', 'FixGuess', 'FitSlope', 'Plot'); 
    % save as labeled variables
    paramout_allfiles = paramout; 
    fval_allfiles = fval; 
    exitflag_allfiles = exitflag; 
    output_allfiles = output; 
   
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % array to save fit parameters
    paramout_array = zeros([size(paramout_allfiles), size(HistTot_array, 3), trials]); 
    fval_array = zeros([size(HistTot_array, 3), trials]); 
    exitflag_array = zeros([size(HistTot_array, 3), trials]); 
    
    fixguess = paramout_allfiles; 
    fixguess(:, 4) = ones(size(fixguess(:,4)))*slope_guess;  
    
    % loop over each file nn
    % generate data fluctuation mm times
    for nn=1:1:size(HistTot_array,3)
        for mm=1:1:trials
            
            % add Poisson fluctuation to data
            tmp = poissrnd(HistTot_array(:,:,nn)); % make this fluctuation on HistTot_array(:,:,nn)
            tof = 1:1:size(HistTot_array, 1); 
            %Convert ToF to Energy using previously calculated Overlap Matrix (OM)
            if (numel(tof) == size(tmp,1))
                Counts = zeros( size(tmp,2), numel(E) );
                for ind = 1:size(tmp,2)
                    %Counts(ind, :) = OM * [tcounts(2:end,ind);0];
                    Counts(ind, :) = OM * circshift(tmp(:,ind),-1);
                end
            elseif (numel(tof) == size(tmp,2))
                Counts = zeros( size(tmp,1), numel(E) );
                for ind = 1:size(tmp,1)
                    %Counts(ind, :) = OM * [tcounts(ind,2:end)';0];
                    Counts(ind, :) = OM * circshift(tmp(ind,:)',-1);
                end
            else
                Counts = 0;
            end
            E_SpectraArray = reshape(Counts.', [E_vec(3) size(HistTot_array,2)]); 
            % fft data
            E_SpectraArray = fftshift(fft(E_SpectraArray, [], 2), 2); 
            twoOmega_signal = E_SpectraArray(:,130);  % hard-coded 2w location....fix this eventually
            
            % fit data
            [paramout, fval, exitflag, output] = ...
                complexfit_section(wavelength, E, twoOmega_signal, fitRegion, 'FixWidth', 'FitSlope'); 
            paramout_array(:,:,nn,mm) = paramout; 
            fval_array(nn,mm) = fval; 
            exitflag_array(nn,mm) = exitflag; 
        
        end % end mm (trials) loop
    end % end nn (files) loop
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do analysis on all of the extracted parameters. 
    
    % mean and variance of samples for each file
    peak_mm_mean  = squeeze(mean(paramout_array(:, 2, :, :), 4)); 
    peak_mm_var   = squeeze(var(paramout_array(:, 2, :, :), [], 4)); 
    phase_mm_mean = squeeze(mean(paramout_array(:, 3, :, :), 4)); 
    phase_mm_var  = squeeze(var(paramout_array(:, 3, :, :), [], 4)); 
    slope_mm_mean = squeeze(mean(paramout_array(:, 4, :,:), 4)); 
    slope_mm_var  = squeeze(var(paramout_array(:, 4, :, :), [], 4)); 
    
    % mean and variance of files from samples mean
    peak_nn  = [mean(peak_mm_mean, 2).'; var(peak_mm_mean, [], 2).']; 
    phase_nn = [mean(phase_mm_mean, 2).'; var(phase_mm_mean, [], 2).']; 
    slope_nn = [mean(slope_mm_mean, 2).'; var(slope_mm_mean, [], 2).']; 

%     if plotting ==1
%         figure; hold on; 
%         errorbar(1:1:size(HistTot_array,3), peak_mm_mean(1,:), peak_mm_var(1,:)); 
%         hold off; 
%         figure; hold on; 
%         errorbar(1:1:size(HistTot_array,3), phase_mm_mean(1,:), phase_mm_var(1,:)); 
%         hold off; 
%         figure; hold on; 
%         errorbar(1:1:size(HistTot_array,3), slope_mm_mean(1,:), slope_mm_var(1,:)); 
%         hold off;
%     end
    if plotting == 1
        figure; hold on; histogram(peak_mm_mean(1,:), 10); hold off; 
        figure; hold on; histogram(phase_mm_mean(1,:), 10); hold off; 
        figure; hold on; histogram(slope_mm_mean(1,:), 10); hold off; 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







end