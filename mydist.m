function yfit_mat = mydist(xdata, width, paramlist)
    yfit = 0; 
    yfit_abs = 0;
    yfit_phi = 0; 
    for i=1:1:length(paramlist(:,1))
        param = paramlist(i,:);
        param(1) = abs(param(1));  
        yfit = yfit + param(1)*exp(-(xdata-param(2)).^2/(2*width^2))*exp(1j*param(3)); 
    end
    
    yfit_mat = [abs(yfit); unwrap(angle(yfit))].'; 
%     yfit_mat = [real(yfit); imag(yfit)].'; 
end