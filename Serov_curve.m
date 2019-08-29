function [curve, wignerFactor] = Serov_curve(Z, E)
    
    E_AU = 27.2114; 
    wavelength = 810; 
    tmpE = E/E_AU; 
    tmpl = wavelength/0.0529; 
    tmpw = 2*pi*137/tmpl; 
    ge = 0.5772;
    
    k = sqrt(2*tmpE); 
    a = Z*tmpw./k.^3; 
    alpha = pi*a/2 - (3*a.^2/2).*(log(2./a) - 1/6 - ge); 
    
    curve = -Z./(1 + alpha)./k.^3 .* (log(2*k.^2/tmpw) - 1 - ge + 3*pi*tmpw*Z/4./k.^3); 
%     curve = -Z./(1 + alpha)./k.^3 .* (log(2*k.^2/tmpw) - 1 - ge + 3*pi*tmpw*Z/4./k.^3); 

    wignerFactor = 1./(1 + alpha); 
%     wignerFactor = 1 - alpha.*(tmpw*Z./k.^3); 
    
end