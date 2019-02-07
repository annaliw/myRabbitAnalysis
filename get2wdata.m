function [E, E_SpectraArray, twoOmega_signal] = get2wdata(HistTot_array, t0, IP, calibType, E_vec, alternate, plotting)   
    
    if ~exist('alternate', 'var')
        alternate = [1, 1]; 
    end
    if ~exist('plotting', 'var')
        plotting = 0; 
    end
    
    % energy calibration based off of previously found peak values
    if strcmp(calibType, 'Kr') == 1
        % calibrate to Kr
        % expected photoelectron energies
        wavelength=810; 
        IPcal = [14, 14.665]; 
        calibEnergy = [((11:1:19)*(1240/wavelength) - IPcal(2)); ((11:1:19)*(1240/wavelength) - IPcal(1))];
        % corresponding peaks in Kr tof data (hand selected, super annoying)
        tof_peak = fliplr(flipud([573 604 639 684 735 803 893 1034 1258; ...
            585 619 657 706 762 840 946 1117 1426])); 
        % calibrate to Kr peaks
        A = ECalibrate(t0, [IPcal(2), IPcal(1)], calibEnergy, tof_peak, plotting);
        % wavelength = wavelength_mod; 
    elseif strcmp(calibType, 'NO') == 1
        % NO self calibrate
        wavelength=810; 
        IPcal = [9.553,16.56,18.319,21.722];
        calibEnergy = [((14:1:19)*(1240/wavelength)-IPcal(2)); ((14:1:19)*(1240/wavelength)-IPcal(1))];  
        tof_peak = [983 860 785 718 675 625; 646 611 579 553 530 510]; 
        A = ECalibrate(t0, [IPcal(2), IPcal(1)], calibEnergy, tof_peak, plotting); 
    else
        error("Unknown calibration type. Choose 'Kr' or 'NO'. ")
    end
    
    % do energy conversion
    tmp = reshape(HistTot_array, size(HistTot_array,1), []);
    %convert the ToF spectrum to energy (linear energy scale)
    [C,E]=Convert_Eng_V2(1:size(tmp,1), tmp, [t0, A] , E_vec);
    E_SpectraArray = reshape(C.', [E_vec(3) size(HistTot_array,2) size(HistTot_array,3)]); 
%     E_SpectraArray_noshift = E_SpectraArray; 

    % only restack the 2w data
    tmp = E_SpectraArray(:,:,alternate(1):alternate(2):end); 
    tmp = fftshift(fft(tmp, [], 2), 2); 
    twoOmega_location = 130; % MAKE THIS AUTOMATICALLY DETECTED
    twoOmega_signal = squeeze(tmp(:,twoOmega_location,:)); 
    % cut out window of ~6 tof bins
    % sideband_list = [(12:2:18)*1240/810 - 13.778, (12:2:18)*1240/810 - 17.706];
    sideband_list = [16*1240/810-IP(1)]; % select 16th harmonic of lowest ionization state 
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

end