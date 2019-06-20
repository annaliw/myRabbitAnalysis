function curve = Ivanov_curve(alpha, E)
    
    E_AU = 27.2114; 
    wavelength = 810; 
    tmpE = E/E_AU; 
    tmpl = wavelength/0.0529; 
    tmpw = 2*pi*137/tmpl; 
    ge = 0.5772;

    curve = -(1./(2*tmpE).^(3/2)).*(log(alpha*2*tmpE/tmpw) - ge + pi*tmpw./(2*alpha*2*tmpE)); 
end