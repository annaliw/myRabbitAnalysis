function h = plotfun_compareToTheory(data, error, vref, energy, theory, name)
    % data: [x,y] by n array, x = energy in eV
    % vref is which vibrational state to reference
    % theory: [x,y] array, x = energy in eV
    % name is the theory label (string)
    
    n = size(data,3); 
    
    set(groot, 'defaultAxesLineWidth', 2); 
    set(groot, 'defaultAxesFontSize', 14); 
    
    h = figure; hold on; 
    set(h, 'Position', [1         399        1440         399]); 
    
    for ii=1:1:n
        [val, ind] = min(abs(energy - data(1,vref,ii))); 
        subplot(1, n, ii); hold on; 
        for jj=1:1:size(theory, 1)
            offset = data(2,vref,ii) - theory(jj,ind); 
            %data
            errorbar(data(1,:,ii), data(2,:,ii), error(2,:,ii), 'o', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); 
            %theory
            plot(energy, theory(jj,:)+offset, 'LineWidth', 2, 'DisplayName', name(jj)); 
            xlim([min(data(1,:,ii))*0.9, max(data(1,:,ii))*1.1]); 
            range = [mean(data(2,:,ii))-abs(max(error(2,:,ii))*3), ...
                mean(data(2,:,ii))+abs(max(error(2,:,ii))*3)]; 
            ylim(range); 
            xlabel('photoelectron energy (eV)'); ylabel('time delay (as)'); 
        end
        hold off; 
    end
    legend; 
    hold off; 
    
    
%     h = 1; 
    

end