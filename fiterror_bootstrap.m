function [] = fiterror_bootstrap('folderName', alternate, t0, wavelength, IP, E_vec, trials, fitRegion, slopeFlag, peakFlag)
    if ~exist('peakflag', 'var')
        peakflag = 0; 
    end
    if ~exist('slopeflag', 'var')
        slopeflag = 0; 
    end
    
    [HistTot_array, stageTimes, freqAxis] = getrawdata('/Users/annaliw/code/2018_07_31-16Scan/', 1, wavelength); 
    
    % calibrate to Krypton
    calibType = 'Kr'; 
    % calibrate to Kr peaks
    A = ECalibrate(t0, calibType, plotting);
    
    % do energy conversion
    tmp = HistTot_array(:,:,1); 
    %convert the ToF spectrum to energy (linear energy scale)
    [C,E,OM]=Convert_Eng_V2(1:size(tmp,1), tmp, [t0, A] , E_vec);
    
    
    
    % arrays to save fit parameters
    paramout_array = zeros([size(
    
    % loop over each file nn
    % generate data fluctuation mm times
    for nn=1:1:size(HistTot_array,3)
        for mm=1:1:trials
            
            % add Poisson fluctuation to data
            tmp = HistTot_array(:,:,nn); % make this fluctuation on HistTot_array(:,:,nn)
            
            %Convert ToF to Energy using previously calculated Overlap Matrix (OM)
            if (numel(tof) == size(tmp,1))
                Counts = zeros( size(tmp,2), numel(Energy) );
                for ind = 1:size(tmp,2)
                    %Counts(ind, :) = OM * [tcounts(2:end,ind);0];
                    Counts(ind, :) = OM * circshift(tmp(:,ind),-1);
                end
            elseif (numel(tof) == size(tmp,2))
                Counts = zeros( size(tmp,1), numel(Energy) );
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
                complexfit_section(IP, E, twoOmega_signal, fitRegion, slope_guess, slopeflag, peakflag); 
        
        end % end mm (trials) loop
    end % end nn (files) loop














end