function h = plotfun_compareToTheory(data, vref, theory, name)
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
        [val, ind] = min(abs(theory(1,:) - data(1,vref,ii))); 
        offset = data(2,vref,ii) - theory(2,ind); 
        subplot(1, n, ii); hold on; 
        %data
        plot(data(1,:,ii), data(2,:,ii)-offset, 'o', 'MarkerFaceColor', 'k', 'DisplayName', 'data'); 
        %theory
        plot(theory(1,:), theory(2,:), 'b', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', name); 
        xlim([min(data(1,:,ii))*0.9, max(data(1,:,ii))*1.1]); 
        xlabel('photoelectron energy (eV)'); ylabel('time delay (as)'); 
        legend; 
        hold off; 
    end
    
    hold off; 
    
    
%     h = 1; 
    

end