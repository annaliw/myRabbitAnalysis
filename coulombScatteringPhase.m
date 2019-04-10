function [E, k, sigma_n] = coulombScatteringPhase(n, wavelength, IP)

    hw = 1240/wavelength; % eV
    E_AU = 27.2114; 
    
    E      = repmat((n)*hw, [numel(IP),1])   - repmat(IP, [numel(n),1])';  % eV
    Eplus  = repmat((n+1)*hw, [numel(IP),1]) - repmat(IP, [numel(n),1])';  % eV
    Eminus = repmat((n-1)*hw, [numel(IP),1]) - repmat(IP, [numel(n),1])';  % eV
    
    k      = sqrt(2*E/E_AU); % au
    kplus  = sqrt(2*Eplus/E_AU); % au
    kminus = sqrt(2*Eminus/E_AU); % au
    
    l = 1; 
    sigma_n = unwrap(angle(igamma(l + 1 - 1j./kplus, 0))) - unwrap(angle(igamma(l + 1 - 1j./kminus, 0))); 
    
    E = E(:)'/E_AU; % au 
    k = k(:)'; % au
    sigma_n = sigma_n(:)'; 
end