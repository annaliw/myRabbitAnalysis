function leastsq = complex_lsq(data, fix, param, slope, peakflag)
    xin = data(1,:); 
    yin = data(2,:); 

    if ~exist('slope', 'var')
        slope = 0; 
    end
    if ~exist('peakflag', 'var')
        peakflag = 0; 
    end

    if length(fix) > 1
        yout = mydist_fixpeaks(xin, fix, param, slope); 
    else
        if peakflag==1
            yout = mydist_fixpeaks(xin, fix, param, slope); 
        else
            yout = mydist_fixwidth(xin, fix, param, slope); 
        end
    end
    
    yout = yout(:,1).*exp(1j*yout(:,2)); 
    yout = yout.'; 

%     leastsq = sum((yin-yout).*conj(yin-yout)./abs(yout)); 
    leastsq = sum((real(yin)-real(yout)).^2) + sum((imag(yin)-imag(yout)).^2); 

end