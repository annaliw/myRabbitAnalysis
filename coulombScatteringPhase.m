function [sigma, delay] = coulombScatteringPhase(l, Z, E)

    E_AU = 27.2114; 
    
    k = sqrt(2*E/E_AU); % au
    
%     sigma_n = unwrap(angle(igamma(l + 1 - 1j./kplus, 0))) - unwrap(angle(igamma(l + 1 - 1j./kminus, 0))); 
    sigma = unwrap(angle(igamma(l + 1 - 1j*Z./k, 0))); 
    delay = 24.2 * diff(sigma)./diff(E/E_AU); 

end