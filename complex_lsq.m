function leastsq = complex_lsq(xin, yin, fix, param, slope, peakflag)

    if ~exist(slope)
        slope = 0; 
    end
    if ~exist(peakflag)
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

    leastsq = sum((yin-yout).*conj(yin-yout)./abs(yout)); 

end