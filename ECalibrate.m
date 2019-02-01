function [A, nshift, wavelength_mod] = ECalibrate(t0, wavelength, calibEnergy, tof_peak)
%     % debug input
%     t0 = 21; 
%     wavelength=810; 
%     folderString = '/Users/annaliw/code/KrCO2_scan/';
%     E_vec = [0, 25, 500];     

    E_vec_max = E_vec(1); 
    E_vec_size = E_vec(3); 

%%

    % prepare for linearized fit of Kr 14 and 14.6 eV peak separations
    tof_fitvar = 1./(tof_peak-t0).^2; 
    dE_over_dx = -0.665./(tof_fitvar(1,:) - tof_fitvar(2,:)); 
    midpt_x = (tof_fitvar(1,:) + tof_fitvar(2,:))/2;
    [paramout, S] = polyfit(midpt_x, dE_over_dx, 1); 
    [yout, delta] = polyval(paramout, midpt_x, S); 

%     % plot fit results
%     figure; hold on; 
%     scatter(midpt_x, dE_over_dx, 'b'); 
%     plot(midpt_x, yout, 'b'); 
%     plot(midpt_x,yout+2*delta,'c--',midpt_x,yout-2*delta,'c--'); 
%     legend('Data','Linear Fit','95% Prediction Interval'); 
%     xlabel('x-midpoint'); 
%     ylabel('\Delta E/\Delta x'); 
%     hold off; 

    % do energy conversion of Kr spectrum
    A = fliplr([paramout(1)/2 paramout(2) 0]); 
    tmp = reshape(HistTot_array, size(HistTot_array,1), []);
    [C,E]=Convert_Eng_V2(1:size(tmp,1), tmp, [t0, A] , E_vec);
    E_SpectraArray = reshape(C.', [E_vec_size size(HistTot_array,2) size(HistTot_array,3)]); 

%%
    %identify and compensate for stage drift
    % will compare location of peak (at tof bin 531 in first subscan) 
    % user input select reference peak
%     figure; plot(E_SpectraArray(:,10,1)); 
%     window_center = round(ginput); 
    window_center = 59; 
    window = 3; 
    histogram_windows = E_SpectraArray(window_center(1)-window:window_center(1)+window,:,:); 
    % integrate over window
    peak_vol = squeeze(sum(histogram_windows, 1)); 
    % fft wrt IR delay 
    twoOmega_location = 130; % from frequency axis calculated previously
    peak_fft = fftshift(fft(peak_vol, [], 1), 1); 
    peak_phase = angle(peak_fft(twoOmega_location, :)); 
    E_SpectraArray_fft = fftshift(fft(E_SpectraArray, [], 2), 2); 
    for ii=1:1:length(peak_phase)
       E_SpectraArray_fft(:,:,ii) = E_SpectraArray_fft(:,:,ii).*exp(-1i*peak_phase(ii)); 
    end
    % sum phase matched values
    E_SpectraArray = sum(E_SpectraArray_fft, 3); 
%%
    % plot Kr Spectrum with harmonics
    IP = [14.0, 14.665];
    IP_label = ["14.665eV", "14eV"]; 
    mode='average'; 
    n = 9:1:19; 
%     plotfun_rabbitspectrum(n, IP, IP_label, 813, E, E_SpectraArray, mode); 

    % diagnostic plots
    tmp1 = A(3)*midpt_x + A(2); 
    tmp2 = A(1) + A(2)*tof_fitvar + A(3)*tof_fitvar.^2; 
    photon_energy(1,:) = tmp2(1,:)+14.665; 
    photon_energy(2,:) = tmp2(2,:)+14; 

%     figure; hold on; scatter(midpt_x, tmp1); title('x_m vs. \Delta E/\Delta x'); 
% 
%     figure; hold on; 
%     s1 = scatter(tof_peak(:), tmp2(:)); s1.MarkerEdgeColor = 'r'; s1.DisplayName = 'fit'; 
%     s2 = scatter(tof_peak(1,:), calibEnergy(1,:)); s2.MarkerEdgeColor = 'b'; s2.DisplayName = '14eV peaks'; 
%     s3 = scatter(tof_peak(2,:), calibEnergy(2,:)); s3.MarkerEdgeColor = 'b'; s3.DisplayName = '14.6eV peaks'; 
%     title('1/t^2 vs. E'); legend; 
% 
%     figure; scatter(midpt_x, tmp2(1,:)-tmp2(2,:)); title('x_m vs. \Delta E'); 
    
% %% test harmonic drift compensation 
%     n = 11:1:19;
% %     figure; hold on; 
% %     s1 = scatter(n, tmp2(1,:)+14.665); s1.MarkerEdgeColor = 'b'; s1.DisplayName = 'E from fit at 14eV tof peaks';
% %     s2 = scatter(n, tmp2(2,:)+14); s2.MarkerEdgeColor = 'c'; s2.DisplayName = 'E from fit at 14.6eV tof peaks';
% %     p1 = line(n, n*1240/wavelength); p1.Color = 'r'; p1.DisplayName = 'expected energy'; 
% %     title('harmonic energy shift'); legend; 
%     photon_energy(1,:) = tmp2(1,:)+14.665; 
%     photon_energy(2,:) = tmp2(2,:)+14; 
%     pe_shift = circshift(photon_energy, -1, 2); 
%     dEbyN = pe_shift(:,1:1:(end-1))-photon_energy(:,1:1:(end-1));
%     nfit = [n(1:1:(end-1)); n(1:1:(end-1))]; 
%     
%     [harmonic_fit,~,mu] = polyfit(nfit, dEbyN, 1); 
% %     harmonic_fit = polyfit(nfit, dEbyN, 1); 
%     n = 9:1:19; 
%     wavelength_mod = 1240/harmonic_fit(2); 
%     nshift = n + n.^2*harmonic_fit(1)/(2*(harmonic_fit(2)-harmonic_fit(1))); 
% %     alpha = harmonic_fit(1)/(2*1240/wavelength); 
% %     n0 = (harmonic_fit(2) - 1240/wavelength)/harmonic_fit(1); 
% %     nshift = n - n0; 
%     
%     figure; hold on;  
%     s1 = scatter(nfit(1,:), dEbyN(1,:)); s1.MarkerEdgeColor = 'b'; s1.DisplayName = 'E(n+1)-E(n) 14eV';
%     s2 = scatter(nfit(2,:), dEbyN(2,:)); s2.MarkerEdgeColor = 'c'; s2.DisplayName = 'E(n+1)-E(n) 14.6eV';
%     p1 = line(n, zeros(size(n))+1240/wavelength); p1.LineStyle = '--'; p1.Color = 'r'; p1.DisplayName = 'expected \Delta E'; 
%     p2 = line(9:1:19, polyval(harmonic_fit, 9:1:19, [], mu)); p2.Color = 'r'; p2.DisplayName = 'fit of \Delta E vs n'; 
% %     p2 = line(0:1:20, polyval(harmonic_fit, 0:1:20)); p2.Color = 'r'; p2.DisplayName = 'fit of \Delta E vs n'; 
%     title('harmonic energy shift'); legend; 
%     xlabel('harmonic number (n)'); 
%     ylabel('\Delta E=E_{n+1}-E_n (eV)'); 
% 
%     E_hfit = [(1240/wavelength)*nshift-14 ; (1240/wavelength)*nshift-14.665]; 
% %     E_hfit = polyval([harmonic_fit(1)/2, harmonic_fit(2), 0], nshift, [], mu); 
%     plotfun_rabbitspectrum(nshift, IP, IP_label, wavelength_mod, E, E_SpectraArray, mode); 
% %     plot(E_hfit, zeros(size(E_hfit)),'o')

%% test2
    harmonic_lin = polyfit(11:1:19, photon_energy(1,:), 1); 
    harmonic_sub = polyval(harmonic_lin, 11:1:19); 
    
%     leftover(1,:) = photon_energy(1,:) - harmonic_sub; 
%     leftover(2,:) = photon_energy(2,:) - harmonic_sub; 
%     
%     leftover_fit1 = fit((11:1:19).', leftover(1,:).', 'power1'); 
%     leftover_fit2 = fit((11:1:19).', leftover(2,:).', 'power1'); 
%     
%     figure; hold on; 
% %     s1 = scatter(11:1:19, leftover(1,:)); s1.MarkerEdgeColor = 'b'; 
% %     s2 = scatter(11:1:19, leftover(2,:)); s2.MarkerEdgeColor = 'b'; 
%     p1 = plot(leftover_fit1, 11:1:19, leftover(1,:)); 
    wavelength_mod = 1240/harmonic_lin(1); 
    nshift = (9:1:19)+harmonic_lin(2)/harmonic_lin(1); 
end