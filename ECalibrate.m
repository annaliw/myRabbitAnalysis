function A = ECalibrate(t0, IP, calibEnergy, tof_peak, plotting)
%     % debug input
%     t0 = 25; 
%     wavelength=810; 
%     IP = [14 14.665]; 
%     % expected photoelectron energies
%     calibEnergy = [((11:1:19)*(1240/wavelength) - 14.00); ((11:1:19)*(1240/wavelength) - 14.665)];
%     % corresponding peaks in Kr tof data (hand selected, super annoying)
%     tof_peak = fliplr(flipud([573 604 639 684 735 803 893 1034 1258; ...
%         585 619 657 706 762 840 946 1117 1426])); 
%     plotting=1;    

    if ~exist('plotting','var')
     % third parameter does not exist, so default it to something
        plotting = 0;
    end

%%

    % prepare for linearized fit of Kr 14 and 14.6 eV peak separations
    tof_fitvar = 1./(tof_peak-t0).^2; 
    dE_over_dx = (min(IP)-max(IP))./(tof_fitvar(1,:) - tof_fitvar(2,:)); 
    midpt_x = (tof_fitvar(1,:) + tof_fitvar(2,:))/2;
    [paramout, S] = polyfit(midpt_x, dE_over_dx, 1); 
    [yout, delta] = polyval(paramout, midpt_x, S); 

    % plot fit results
    if plotting==1
        figure; hold on; 
        scatter(midpt_x, dE_over_dx, 'b'); 
        plot(midpt_x, yout, 'b'); 
        plot(midpt_x,yout+2*delta,'c--',midpt_x,yout-2*delta,'c--'); 
        legend('Data','Linear Fit','95% Prediction Interval'); 
        xlabel('x-midpoint'); 
        ylabel('\Delta E/\Delta x'); 
        hold off; 
    end

    % do energy conversion of Kr spectrum
    A = fliplr([paramout(1)/2 paramout(2) 0]); 
%     Eng = 0;
%     for i = 0:numel(A)-1
%         Eng = Eng + A(i+1) .* tof_fitvar(:).^i;
%     end 
    
    [peaksinorder, order] = sort(tof_peak(:)); 
    Eng = A(3)*tof_fitvar.^2 + A(2)*tof_fitvar + A(1); 
    Eng = Eng(:); 

    if plotting==1
        figure; hold on; 
        s1 = scatter(tof_peak(1,:), calibEnergy(1,:)); s1.MarkerEdgeColor = 'b'; s1.DisplayName = 'b'; 
        s2 = scatter(tof_peak(2,:), calibEnergy(2,:)); s2.MarkerEdgeColor = 'r'; s2.DisplayName = 'X'; 
        p1 = plot(peaksinorder, Eng(order)); p1.Color = 'c';   
        title('Energy Calibration Result'); 
        xlabel('time of flight'); 
        ylabel('photoelectron energy'); 
        legend; 
        hold off; 
    end
%     tmp = reshape(HistTot_array, size(HistTot_array,1), []);
%     [C,E]=Convert_Eng_V2(1:size(tmp,1), tmp, [t0, A] , E_vec);
%     E_SpectraArray = reshape(C.', [E_vec_size size(HistTot_array,2) size(HistTot_array,3)]); 

end