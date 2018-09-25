function chi2=myfit(xin,yin,paramlist)
    
    yfit = mydist(xin, paramlist); 
    
    chi2_r = (real(yin) - real(yfit)).^2/real(yin).^2; 
    chi2_i = (imag(yin) - imag(yfit)).^2/imag(yin).^2; 

    chi2 = chi2_r + chi2_i; 
    
end