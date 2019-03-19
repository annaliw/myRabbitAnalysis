function trash = plotfun_rabbitspectrum(n, wavelength, E, data, mode)

    global IP
    global IP_label
    
    checksize = size(data); 
    % text and color settings
    text_size = 12; 
    line_weight = 1.5; 
    abs_color = [0 0 0.8]; 
    phi_color = [0 0.7 1]; 
    
%     n = 9:1:19; 

    fh = figure('Position', [10 10 400*2.5 400]); 
    hold on; 
    
    if strcmp(mode,'twoOmega')==1
        
        % plotting inputs
        if ismember(1, checksize)==0
            oneOmega_signal = data(:,121); 
            twoOmega_signal = data(:,130); 
        else
            twoOmega_signal = data; 
        end
        twoOmega_abs = abs(twoOmega_signal); 
%         twoOmega_phi = unwrap(angle(twoOmega_signal)); 
        twoOmega_phi = angle(twoOmega_signal); 
        % oneOmega_abs = abs(oneOmega_signal); 
        % oneOmega_phi = angle(oneOmega_signal); 
    
        ax1 = gca; 
        ax1_pos = ax1.Position; 
        ax1.XColor = 'k'; 
        ax1.XLabel.String = 'Photoelectron Energy (eV)'; 
        ax1.YLabel.String = 'Amplitude (arbitrary units)';
        ax1.FontSize = text_size; 
        ax1.YColor = abs_color; 
        ax1.LineWidth = line_weight*0.5; 
        line(E, twoOmega_abs, ...
            'Parent', ax1, 'Color', abs_color, 'LineWidth', line_weight, ...
            'DisplayName', '2w amplitude');
%         line(E, 14*mean(abs(E_SpectraArray), 2)./mean((mean(abs(E_SpectraArray),2))), ...
%             'Parent', ax1, 'Color', 'k', 'LineWidth', line_weight, ...
%             'DisplayName', 'time delay average'); 
        % plot(E, movmean(oneOmega_abs, 3), 'r-', 'DisplayName', '1w amplitude'); 

        ax2 = axes('Position', ax1_pos, ...
            'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
            'Color', 'none', 'LineWidth', line_weight*0.5); 
        ax2.XColor = 'none'; 
        ax2.YColor = phi_color; 
        ax2.YLabel.String = 'phase (radians)'; 
        ax2.FontSize = text_size; 
        line(E, twoOmega_phi, ...
            'Parent', ax2, 'Color', phi_color, 'LineWidth', line_weight, ...
            'DisplayName', '2w phase'); 
        % plot(E, unwrap(oneOmega_phi), 'm-', 'DisplayName', '1w phase'); 
    else
        ax1 = gca; 
        line(E, mean(abs(data), 2), 'Color', 'k'); 
        ax1.XLabel.String = 'photoelectron energy (eV)'; 
    end
    
    axl = AddHarmonicAxis(ax1, IP, wavelength);

    for i = 1:numel(IP)
        axl(i).XLabel.Position = [ -1.9    0.99];
        axl(i).FontSize = text_size*0.75; 
        axl(i).XLabel.String = IP_label(i); 
    end
%     axl(4).XLabel.Position = [-1.9 14.8]; % figure out why this one is different
%     hold on; 
    
    
%     legend({}, 'FontSize', text_size); 

    xlim([min(E(end), E(1)) max(E(end),E(1))]); 
    if strcmp(mode, 'twoOmega')==1
        linkaxes([ax1,ax2],'x')
    end

    hold off;  
    
    trash=1; 
    
end