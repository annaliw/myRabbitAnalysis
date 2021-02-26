function times = HHG_trajectories(wavelength, fwhm, intensity)

    times = []; % initialize list of time values to save
    c = 0.2998;  % speed of light um/fs
    wavelength = wavelength * 10^-3; % convert nm to um
    sigma = fwhm ./ (2*sqrt(2*log(2))); 
    
    % make laser field
    t = -50:0.01:50; % time array in fs
    w = 2*pi * c / wavelength; 
    E_field = sqrt(intensity).*sin(w.*t) .* exp(-t.^2./(2.*sigma.^2)); 
    p = cumtrapz(t, E_field); % momentum (v) = integral of field (acceleration)
    figure; hold on; 
    plot(t, E_field, 'b-'); 
    plot(t, p, 'r-'); 
    
    % correct momentum values to boundary conditions
    % p(ti) = 0
    % x(ti) = x(tr) = 0
    for ii=1:numel(t)
        
        ti = t(ii); 
        p_tmp = p - p(ii); % solve boundary condition p(ti) = 0
        
        x = cumtrapz(t, p_tmp); % integrate for position
        x = x - x(ii); % solve boundary condition x(ti) = 0
        
        for jj=1:numel(t)
            
            tr = t(jj); 
            if abs(x(jj)) < 0.01
                times = [times; ti tr]; % save on condition x(tr) = 0
            end
            
        end
        
    end






end