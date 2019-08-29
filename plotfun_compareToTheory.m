function offset_array = plotfun_compareToTheory(data, error, energy, theory, name, fignum)
    % data: [x,y] by n array, x = energy in eV
    % vref is which vibrational state to reference
    % theory: [x,y] array, x = energy in eV
    % name is the theory label (string)
    
    n = size(data,3); 
    offset_array = zeros([n size(theory,1)]); 
    
    set(groot, 'defaultAxesLineWidth', 2); 
    set(groot, 'defaultAxesFontSize', 14); 
    
    h = figure(fignum); hold on; 
    set(h, 'Position', [1         399        1440         399]); 
    axh = axes; 
    
    for ii=1:1:n
        % get theory range for section
        [val1, ind1] = min(abs(energy - min(data(1,:,ii)))); 
        [val2, ind2] = min(abs(energy - max(data(1,:,ii)))); 
        figure(fignum); 
        subplot(1, n, ii); hold on; 
        %data
        errorbar(data(1,:,ii), data(2,:,ii), error(2,:,ii), 'o', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); 
        xlim([min(data(1,:,ii))*0.9, max(data(1,:,ii))*1.1]); 
        range = [mean(data(2,:,ii))-abs(max(error(2,:,ii))*5), ...
            mean(data(2,:,ii))+abs(max(error(2,:,ii))*5)]; 
        ylim(range); 
        xlabel('photoelectron energy (eV)'); ylabel('time delay (as)'); 
        for jj=1:1:size(theory, 1)
            figure(fignum); 
            subplot(1, n, ii); 
            offset_array(ii, jj) = mean(data(2,:,ii)) - mean(theory(jj,ind1:ind2)); 
            %theory
            plot(energy, theory(jj,:)+offset_array(ii, jj), 'LineWidth', 2, 'DisplayName', name(jj)); 
            
            figure(fignum+1); hold on; 
            subplot(1, size(theory,1), jj); hold on; 
            errorbar(data(1,:,ii), data(2,:,ii)-offset_array(ii,jj), error(2,:,ii), 'HandleVisibility', 'off'); 
            if ii==1
                plot(energy, theory(jj,:), 'k', 'DisplayName', name(jj)); 
                xlim([min(min(data(1,:,:))), max(max(data(1,:,:)))]); 
                legend;
                xlabel('photoelectron energy (eV)'); ylabel('delay (as)'); 
            end
        end
    end
    
    figure(fignum)
    legend; 
    hold off; 
    
    
%     h = 1; 
    

end