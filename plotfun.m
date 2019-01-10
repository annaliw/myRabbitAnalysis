function trash = plotfun(xin, yin_abs, yin_phi, x_out, yout_abs, yout_phi, peaks, paramout)

    yyaxis left
    ylabel('amplitude'); 
    scatter(xin, yin_abs, 'o', 'DisplayName', 'amplitude data'); 
    plot(x_out, yout_abs, '-', 'DisplayName', 'amplitude fit'); 
    for i=1:1:(length(peaks))
        tmp = mydist(x_out, peaks(i), paramout(i,:)); 
        tmp_abs = abs(tmp(:,1).*exp(1j*tmp(:,2))); 
        plot(x_out, tmp_abs, '--', 'HandleVisibility', 'off'); 
    end
    yyaxis right
    ylabel('phase'); 
    scatter(xin, yin_phi, '+', 'DisplayName', 'phase data'); 
    plot(x_out, yout_phi, '-', 'DisplayName', 'phase fit'); 
    for i=1:1:(length(peaks))
        tmp = mydist(x_out, peaks(i), paramout(i,:));
        tmp_phi = angle(tmp(:,1).*exp(1j*tmp(:,2))); 
        plot(x_out, unwrap(tmp_phi), '--', 'HandleVisibility', 'off'); 
    end
    
    xlim([x_out(1), x_out(end)])
    legend; 
    xlabel('photoelectron energy (eV)'); 
    
    trash=1; 
    
end