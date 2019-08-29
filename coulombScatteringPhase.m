function [sigma, delay] = coulombScatteringPhase(l, Z, E)

    E_AU = 27.2114; 
    
    k = sqrt(2*E/E_AU); % au
    
%     sigma_n = unwrap(angle(igamma(l + 1 - 1j./kplus, 0))) - unwrap(angle(igamma(l + 1 - 1j./kminus, 0))); 
%     sigma = unwrap(angle(igamma(l + 1 + 1j*Z./k, 0)./igamma(l + 1 - 1j*Z./k, 0)))/2; 
    sigma = angle(igamma(l + 1 - 1j*Z./k, 0)); 
    sigma = fliplr(unwrap(fliplr(sigma))); 

%     sigma = (l+0.5)*atan(Z./k/(l+1)) + Z./k .* (log(sqrt((l+1)^2 + (Z./k).^2))-1); 
%     sigma = sigma - (Z./k)./(12*((l+1)^2+(Z./k).^2)); 

    delay = 24.2 * diff(sigma)./diff(E/E_AU); 

end