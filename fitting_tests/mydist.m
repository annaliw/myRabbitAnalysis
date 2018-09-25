function yfit = mydist(xin, paramlist)
    yfit = 0; 
    for i=1:1:length(paramlist(:,1))
        param = paramlist(i,:); 
        yfit = yfit + param(1)*exp(-(xin-param(2)).^2/param(3)^2)*exp(1j*param(4)); 
    end
end