function trash = plotfun(xin, yin_abs, yin_phi, x_out, yout_abs, yout_phi, peaks_guess, paramout)

    yyaxis left
    scatter(xin, yin_abs, 'o'); 
    plot(x_out, yout_abs, '-'); 
    for i=1:1:(length(peaks_guess))
        tmp = mydist(x_out, peaks_guess(i), paramout(i,:)); 
        tmp_abs = abs(tmp(:,1) + 1j*tmp(:,2)); 
        plot(x_out, tmp_abs, '--'); 
    end
    yyaxis right
    scatter(xin, yin_phi, '+'); 
    plot(x_out, yout_phi, '-'); 
    for i=1:1:(length(peaks_guess))
        tmp = mydist(x_out, peaks_guess(i), paramout(i,:));
        tmp_phi = angle(tmp(:,1) + 1j*tmp(:,2)); 
        plot(x_out, unwrap(tmp_phi), '--'); 
    end
    xlim([x_out(1), x_out(end)])
    
    trash=1; 
    
end