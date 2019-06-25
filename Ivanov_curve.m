function curve = Ivanov_curve(l, Z, E)
    
    E_AU = 27.2114; 
    wavelength = 810; 
    tmpE = E/E_AU; 
    tmpl = wavelength/0.0529; 
    tmpw = 2*pi*137/tmpl; 
    ge = 0.5772; 
    
    k = sqrt(2*tmpE); 
 
    [sigma, delay] = coulombScatteringPhase(l, Z, E); 
    sigma = unwrap(sigma); 
    a = 2*exp(-2*sigma.*k); 

    curve = -(1./((2*tmpE).^(3/2)+tmpw*pi/2)).*(log(a.*2.*tmpE/tmpw) - ge + pi*tmpw./(2*a.*2.*tmpE)); 
end