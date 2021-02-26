function A = ECalibrate(t0, n, wavelength, calibType, config)

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

    % get peak separation fit
    hw = 1240/wavelength; 
    tof_fitvar = 1./(tof_peak-t0).^2; 
    if size(tof_peak,1) == 1
        % dx
        tmp = abs(tof_fitvar - circshift(tof_fitvar, 1)); 
        dE_over_dx = abs((calibEnergy - circshift(calibEnergy,1)))./tmp; 
        dE_over_dx = dE_over_dx(2:end); 
        % midpt
        tmp = tof_fitvar + circshift(tof_fitvar, 1);
        midpt_x = (tmp(2:end))/2;
    elseif size(tof_peak,1) == 2 % Kr
        % dx  
        tmp = abs(tof_fitvar(1,:) - tof_fitvar(2,:)); 
        dE_over_dx = abs(calibEnergy(1,:) - calibEnergy(2,:))./tmp; 
        % midpt
        tmp = tof_fitvar(1,:) + tof_fitvar(2,:); 
        midpt_x = tmp/2; 
    else
        error('incorrect input size'); 
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
    if size(tof_peak,1) == 1
        ydata = n*hw; 
    elseif size(tof_peak,1) == 2
        ydata = [n*hw; n*hw]; 
        xdata = xdata(:).'; 
        ydata = ydata(:).'; 
    else
        error('incorrect input size'); 
    end
    
    [paramout_wavelength, S] = polyfit(xdata, ydata, 1); 
    [yout, delta] = polyval(paramout_wavelength, xdata, S); 
    
    % plot fit results
    if plotting==1
        figure; hold on; 
        ax = gca; ax.FontSize=12; 
        scatter(xdata, ydata, 'b'); 
        plot(xdata, yout, 'b'); 
        plot(xdata,yout+2*delta,'c--',xdata,yout-2*delta,'c--'); 
        legend('Data','Linear Fit','95% Prediction Interval'); 
        xlabel('xdata'); 
        ylabel('energy'); 
        hold off; 
    end
    
    % calculate new wavelength
%     hw = hw/paramout_wavelength(1); 
%     wavelength = 1240/hw; 
    % calculate offset
    A(1) = paramout_wavelength(2) - mean(IPcal); 
    
    [peaksinorder, order] = sort(tof_peak(:)); 
    
    Eng = A(1) + A(2)*tof_fitvar + A(3)*tof_fitvar.^2; 
    Eng = Eng(:).'; 

    if plotting==1
        figure; hold on; 
        ax = gca; ax.FontSize=12; 
        s1 = scatter(tof_peak(:).', calibEnergy(:).'); s1.MarkerEdgeColor = 'b'; s1.DisplayName = 'b'; 
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