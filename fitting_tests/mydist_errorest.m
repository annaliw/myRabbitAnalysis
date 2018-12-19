function yfit = mydist(xin, x0, paramlist)
    yfit = 0; 
    for i=1:1:length(paramlist(:,1))
        param = paramlist(i,:); 
        yfit = yfit + abs(param(1))*exp(-(xin-x0).^2/param(3)^2)*exp(1j*+ param(4)); 
    end
end