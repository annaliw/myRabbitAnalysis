function phi_cc = continuumphase(k2, k1, Z)

    phi_cc = (2*k1)^(1j*Z*k1)/(2*k2)^(1j*Z*k2); 
    phi_cc = phi_cc * gamma(2 + 1j*Z*(1/k1 - 1/k2))/(k1 - k2)^(1j*Z*(1/k1 - 1/k2)); 
    
    return phi_cc

end