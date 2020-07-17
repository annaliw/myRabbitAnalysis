function A = ECalibrate(t0, n, wavelength, calibType, config)
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
    
%     % energy calibration based off of previously found peak values
%     if strcmp(calibType, 'Kr') == 1
%         % calibrate to Kr
%         % expected photoelectron energies
%         wavelength=810; 
%         IPcal = [14, 14.665]; 
%         calibEnergy = [((11:1:19)*(1240/wavelength) - IPcal(2)); ((11:1:19)*(1240/wavelength) - IPcal(1))];
%         % corresponding peaks in Kr tof data (hand selected, super annoying)
%         tof_peak = fliplr(flipud([573 604 639 684 735 803 893 1034 1258; ...
%             585 619 657 706 762 840 946 1117 1426])); 
%         % calibrate to Kr peaks
%         % wavelength = wavelength_mod; 
%     elseif strcmp(calibType, 'NO') == 1
%         % NO self calibrate
%         wavelength=810; 
%         IPcal = [9.553,16.56];
%         calibEnergy = [((14:1:19)*(1240/wavelength)-IPcal(2)); ((14:1:19)*(1240/wavelength)-IPcal(1))];  
%         tof_peak = [983 860 785 718 675 625; 646 611 579 553 530 510]; 
%     else
%         if (isfield(config,'tofPeaks')) % Peak positions:
%         tof_peak = config.tofPeaks; 
%         end
%         if (isfield(config,'IPcal')) % Peak positions:
%             IPcal = config.IPcal; 
%         end
%         if (isfield(config,'calibEnergy')) % Peak positions:
%             calibEnergy = config.calibEnergy; 
%         end
% 
%         if (isfield(config,'Plot')) % Peak positions:
%             plotting = config.Plot; 
%         else
%             plotting = 0; 
%         end
%         
%     end


if (isfield(config,'tofPeaks')) % Peak positions:
    tof_peak = config.tofPeaks; 
end
if (isfield(config,'IPcal')) % Peak positions:
    IPcal = config.IPcal; 
end
if (isfield(config,'calibEnergy')) % Peak positions:
    calibEnergy = config.calibEnergy; 
end

if (isfield(config,'Plot')) % Peak positions:
    plotting = config.Plot; 
else
    plotting = 0; 
end


    
% %%
%     % debug input Ar
%     t0 = 0; 
%     wavelength=810; 
%     IPcal = [15.763]; 
%     n = 11:1:21; 
% %     calibEnergy = [n*1240/wavelength - IP]; 
%     calibEnergy = n*hw - IP; 
% %     tof_peak = [576 601 636 672 719 770 840 934 1082 1339 1968];  % redo with peak finding 
%     tof_peak = fliplr([577 603 634 673 716 771 838 936 1079 1335 1962]); 
%     plotting=1; 
%     
% %     % debug input Kr
% %     IPcal = [14, 14.665]; 
% %     calibEnergy = [((11:1:19)*(1240/wavelength) - IPcal(2)); ((11:1:19)*(1240/wavelength) - IPcal(1))];
% %     % corresponding peaks in Kr tof data (hand selected, super annoying)
% %     tof_peak = fliplr(flipud([573 604 639 684 735 803 893 1034 1258; ...
% %         585 619 657 706 762 840 946 1117 1426])); 
% %     calibEnergy = calibEnergy(1,:); 
% %     tof_peak = tof_peak(1,:); 
% %     IPcal = IPcal(2); 

    %% get peak separation fit
    hw = 1240/wavelength; 
    % prepare for linearized fit of Kr 14 and 14.6 eV peak separations
    tof_fitvar = 1./(tof_peak-t0).^2; 
    tmp = abs(tof_fitvar - circshift(tof_fitvar, 1)); 
    dE_over_dx = abs((calibEnergy - circshift(calibEnergy,1)))./tmp; 
    dE_over_dx = dE_over_dx(2:end); 
    tmp = tof_fitvar + circshift(tof_fitvar, 1);
    midpt_x = (tmp(2:end))/2;
    
    if strcmp(calibType, 'Kr') == 1
        midpt_x = midpt_x(1:2:end); 
        dE_over_dx = dE_over_dx(1:2:end); 
        n = repmat(n, [size(IPcal,2) size(IPcal,1)]); 
        n = n(:); 
    end
    
    [paramout, S] = polyfit(midpt_x, dE_over_dx, 1); 
    [yout, delta] = polyval(paramout, midpt_x, S); 

    % plot fit results
    if plotting==1
        figure; hold on; 
        ax = gca; ax.FontSize=12; 
        scatter(midpt_x, dE_over_dx, 'b'); 
        plot(midpt_x, yout, 'b'); 
        plot(midpt_x,yout+2*delta,'c--',midpt_x,yout-2*delta,'c--'); 
        legend('Data','Linear Fit','95% Prediction Interval'); 
        xlabel('x-midpoint'); 
        ylabel('\Delta E/\Delta x'); 
        hold off; 
    end
    
    %% get wavelength fit
    
    A = fliplr([0.5*paramout(1) paramout(2) 0]); 
    xdata = A(3)*tof_fitvar.^2 + A(2)*tof_fitvar; 
    ydata = n*hw; 
    
    [paramout_wavelength, S] = polyfit(xdata, ydata, 1); 
    [yout, delta] = polyval(paramout_wavelength, xdata, S); 
    
%     % plot fit results
%     if plotting==1
%         figure; hold on; 
%         ax = gca; ax.FontSize=12; 
%         scatter(xdata, ydata, 'b'); 
%         plot(xdata, yout, 'b'); 
%         plot(xdata,yout+2*delta,'c--',xdata,yout-2*delta,'c--'); 
%         legend('Data','Linear Fit','95% Prediction Interval'); 
%         xlabel('xdata'); 
%         ylabel('energy'); 
%         hold off; 
%     end
    
    % calculate new wavelength
%     hw = hw/paramout_wavelength(1); 
%     wavelength = 1240/hw; 
    % calculate offset
    A(1) = paramout_wavelength(2) - mean(IPcal); 
    
    [peaksinorder, order] = sort(tof_peak(:)); 
    
    Eng = A(1) + A(2)*tof_fitvar + A(3)*tof_fitvar.^2; 

    if plotting==1
        figure; hold on; 
        ax = gca; ax.FontSize=12; 
        s1 = scatter(tof_peak, calibEnergy); s1.MarkerEdgeColor = 'b'; s1.DisplayName = 'b'; 
        p1 = plot(peaksinorder, Eng(order)); p1.Color = 'c'; p1.DisplayName = 'calibration';  
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