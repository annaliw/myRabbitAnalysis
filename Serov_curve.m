function curve = Serov_curve(alpha, E)
    
    E_AU = 27.2114; 
    wavelength = 810; 
    tmpE = E/E_AU; 
    tmpl = wavelength/0.0529; 
    tmpw = 2*pi*137/tmpl; 
    ge = 0.5772;

    curve = -(1./(1 + alpha*tmpw./((2*tmpE).^(3/2)))).*(1./(2*tmpE).^(3/2)).*(log(4*tmpE/tmpw) - ge - 1 + 3*pi*tmpw./(4*(2*tmpE).^(3/2))); 
end