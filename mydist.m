function yfit_mat = mydist(xdata, peaks, paramlist)
    yfit = 0; 
    yfit_abs = 0;
    yfit_phi = 0; 
    for i=1:1:length(peaks)
        param = paramlist(i,:);
        param(1) = abs(param(1)); 
%         param(3) = unwrap(param(3)); 
        yfit = yfit + param(1)*exp(-(xdata-peaks(i)).^2/(2*param(2)^2))*exp(1j*param(3)); 
    end
    
    yfit_mat = [abs(yfit); mod(unwrap(angle(yfit)), 2*pi)].'; 
%     yfit_mat = [real(yfit); imag(yfit)].'; 
end