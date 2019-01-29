function trash = plotfun_fit(IP, wavelength, xin, yin, peaks, paramout)
    % text and color settings
    text_size = 12; 
    line_weight = 1.5; 
    abs_color = [0 0 0.8]; 
    phi_color = [0 0.7 1]; 
    lgnd_pad = 1.4; 
    lgnd_pos = 'southeast'; 
    
    % plotting input
    yin_abs = yin(:,1); 
    yin_phi = mod(yin(:,2), 2*pi); 

    xout = linspace(xin(1,1),xin(1,end),length(xin(1,:))*100);
    yout = mydist(xout, peaks, paramout); 
    yout_abs = yout(:,1); 
    yout_phi = mod(yout(:,2), 2*pi); 

    
    fh = figure('Position', [10 600 560 400]); 
    hold on; 
    % amplitude axis
    ax1 = gca; 
    ax1_pos = ax1.Position; 
    ax1.XColor = 'k'; 
    ax1.XLabel.String = 'Photoelectron Energy (eV)'; 
    ax1.YLabel.String = 'Amplitude (arbitrary units)';
    ax1.FontSize = text_size; 
    ax1.YColor = abs_color; 
    ax1.LineWidth = line_weight*0.5; 
    hold(ax1, 'on'); 
    s1 = scatter(ax1, xin, yin_abs, 'o'); 
    s1.MarkerEdgeColor = abs_color; 
    s1.LineWidth = line_weight; 
    s1.DisplayName = 'amplitude data'; 
    l1 = line(xout, yout_abs, 'Parent', ax1, ...
        'Color', abs_color, 'LineStyle', '-', 'LineWidth', line_weight, ...
        'DisplayName', 'total amplitude fit'); 
    for i=1:1:(length(peaks))
        tmp = mydist(xout, peaks(i), paramout(i,:)); 
%         tmp_abs = abs(tmp(:,1).*exp(1j*tmp(:,2))); 
        tmp_abs = tmp(:,1); 
        line(xout, tmp_abs, 'Parent', ax1, ...
            'Color', abs_color, 'LineStyle', '--', 'LineWidth', line_weight, ...
            'HandleVisibility', 'off'); 
    end
    
    % phase axis
    ax2 = axes('Position', ax1_pos, ...
        'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
        'Color', 'none', 'LineWidth', line_weight*0.5); 
    ax2.XColor = 'none'; 
    ax2.YColor = phi_color; 
    ax2.YLabel.String = 'phase (radians)'; 
    ax2.FontSize = text_size; 
    hold(ax2, 'on'); 
    s2 = scatter(ax2, xin, yin_phi, '+'); 
    s2.MarkerEdgeColor = phi_color; 
    s2.LineWidth = line_weight; 
    s2.DisplayName = 'phase data'; 
    l2 = line(xout, yout_phi, 'Parent', ax2, ...
        'Color', phi_color, 'LineStyle', '-', 'LineWidth', line_weight, ...
        'DisplayName', 'total phase fit'); 
    for i=1:1:(length(peaks))
        tmp = mydist(xout, peaks(i), paramout(i,:));
%         tmp_phi = angle(tmp(:,1).*exp(1j*tmp(:,2))); 
        tmp_phi = mod(tmp(:,2), 2*pi); 
        line(xout, tmp_phi, 'Parent', ax2, ...
            'Color', phi_color, 'LineStyle', '--', 'LineWidth', line_weight, ...
            'HandleVisibility', 'off'); 
    end
    
    linkaxes([ax1,ax2],'x')
%     xlabel('photoelectron energy (eV)'); 


    % harmonic labels
    axl = AddHarmonicAxis(ax1,IP, wavelength);
    axl(1).XLabel.String = 'X (avg v. 0-2)';
    axl(2).XLabel.String = 'b ^3\Pi';
    axl(3).XLabel.String = 'A^1\Pi';
    axl(4).XLabel.String = 'c^3\Pi'; 

    for i = 1:numel(IP)
        axl(i).XLabel.Position = [ -1.9    0.99];
        axl(i).FontSize = text_size*0.75; 
    end
    ax1.XLim = [xout(1), xout(end)]; 
    ax1.YLim = [0, max(yout_abs)*lgnd_pad]; 
    lgnd = legend([s1, l1, s2, l2], 'Location', lgnd_pos); 
    lgnd.FontSize = text_size*0.6; 
    
    hold off; 
    
    trash=1; 
    
end