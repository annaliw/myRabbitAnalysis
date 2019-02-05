function SI = spectral_integration(xin, yin_amp, yin_phi, peaks, window)

    phase = zeros([1 length(peaks)]); 
%     windowcheck = zeros([2 length(peaks)]); 
    for i=1:1:length(peaks)
        idx_minus = find(abs(xin-peaks(i)+window)<0.01,1); 
        idx_plus = find(abs(xin-peaks(i)-window)<0.01,1); 
%         windowcheck(:,i) = [idx_minus idx_plus]; 

        phase(i) = mod(angle(sum(yin_amp(idx_minus:idx_plus).*exp(1j*yin_phi(idx_minus:idx_plus)))), 2*pi()); 
    end
    
    SI = phase; 
end