function yout = gaussian_mixture_fixsigma(x, xdata)
    
%     if neq(size(x,2), 3)
%         error('Incorrect parameter matrix dimensions. There are 3 parameters.')
%     end
    
    x(:,3) = x(1,3); 
    yout = sum(x(:,1).*exp(-(xdata-x(:,2)).^2./(2.*x(:,3)).^2), 1);


end