function yfit_mat = mydist_fixpeaks(xdata, peaks, paramlist, slope)
    yfit = 0; 
    yfit_abs = 0;
    yfit_phi = 0; 
    if ~exist('slope','var')
     % third parameter does not exist, so default it to something
        slope = 0;
    end
    for i=1:1:length(paramlist(:,1))
        param = paramlist(i,:);
        param(1) = abs(param(1));  
        if slope==0
            yfit = yfit + param(1)*exp(-(xdata-peaks(i)).^2/(2*param(2)^2))...
                .*exp(1j*param(3)); 
        else
            yfit = yfit + param(1)*exp(-(xdata-peaks(i)).^2/(2*param(2)^2))...
                .*exp(1j*(param(3)+param(4)*(xdata-peaks(i)))); 
        end
    end
    
    yfit_mat = [abs(yfit); unwrap(angle(yfit))].'; 
%     yfit_mat = [real(yfit); imag(yfit)].'; 
end