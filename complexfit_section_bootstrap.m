function [paramout, fval] = complexfit_section_full(wavelength, xin, yin, paramout_gauss, paramout_original, plotting)
    global IP; 
    global IP_label; 
    
%     region = fitRegion; 
    n = 9:1:19; 
    tolerance = 0.02; 

%     signal = data'; 
%     x = E; 
%     y = signal; 
% 
%     % fit section set-up
%     start = find(abs(x-region(1))<tolerance, 1, 'last'); 
%     stop = find(abs(x-region(2))<tolerance, 1, 'first'); 
% 
%     xin = x(start:stop); 
%     yin = y(start:stop); 
    yin_abs = abs(yin)./sum(abs(yin)); 
    yin_phi = mod(unwrap(angle(yin)),2*pi); 
     

    % find new peak center indices
    % paramout_gauss(:,2) = peaks_guess'; 
    peak_ind = 1:1:length(paramout_gauss(:,2)); 
    for i=1:1:length(peak_ind)
        peak_ind(i) = find(abs(xin - paramout_gauss(i,2)) < tolerance, 1);  
    end

    %% fit complex 2w data using gaussians found in above cell

    % fit new gaussian amplitudes
    yin = yin_abs; 
    guess = ones([1 size(paramout_gauss,1)]); 
    [newamp, fval] = fitGaussianAmplitudes(xin, yin, paramout_gauss, guess, 0); 
    paramout_gauss(:,1) = newamp; 
    
    
    
    yin = [yin_abs; yin_phi];
    guess = paramout_original; 
    [paramout, fval] = fit2OmegaSum(xin, yin, paramout_gauss, guess, plotting); 
    % paramout(:,2) = mod(paramout(:,2), 2*pi); 

    % plotfun_fit(n, 810, xin, yin, fix, paramout, slope, peakflag)
    
    paramout = cat(2, paramout_gauss, paramout); 
    

end