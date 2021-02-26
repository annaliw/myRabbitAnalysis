function [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E, avdata, data, slope, plotting)
    global IP; 
    global IP_label; 
    %% fit 2w abs data for region
    
%     region = fitRegion; 
    n = 1:1:19; 
    tolerance = abs(E(1)-E(2)); 
%     n = 0; 
%     tolerance = abs(E(1)-E(2))*2; 
    

    % form peaks guess
    peaks = zeros(length(IP), length(n)); 
    for i=1:1:length(IP)
        peaks(i,:) = n*(1240/wavelength)-IP(i); 
    end

    signal = avdata'; 
    xin = E; 
    yin = signal; 
    peaks = zeros(length(IP), length(n)); 
    for i=1:1:length(IP)
        peaks(i,:) = n*(1240/wavelength)-IP(i); 
    end
    peaks = peaks(:).'; 

%     % fit section set-up
%     start = find(abs(x-region(1))<tolerance, 1, 'last'); 
%     stop = find(abs(x-region(2))<tolerance, 1, 'first'); 
% 
%     xin = x(start:stop); 
%     yin = y(start:stop); 
    yin_abs = abs(yin); %./sum(abs(yin)); 
    yin = yin_abs; 

    % find peaks and their indices
    peaks_guess = peaks; 
    remove = [find(peaks_guess < xin(1)) find(peaks_guess > xin(end))]; 
    peaks_guess(remove) = []; 
    peak_ind = 1:1:length(peaks_guess); 
    for i=1:1:length(peaks_guess)
        peak_ind(i) = find(abs(xin - peaks_guess(i)) < tolerance, 1);  
    end
    amp_guess = yin(peak_ind); % form amplitude guess
    sig_guess = ones(size(amp_guess))*0.1; % form width guess
%     % special, remove later:
%     sig_guess(end) = 0.1; 
    alpha_guess = ones(size(amp_guess))*(3); 
%     guess = [amp_guess; peaks_guess; sig_guess; alpha_guess].'; % full guess matrix
    guess = [amp_guess; peaks_guess; sig_guess].'; 

    [paramout_gauss, fval_gauss] = fitGaussianSum(xin, yin, guess, 0); 
    % find new peak center indices
    % paramout_gauss(:,2) = peaks_guess'; 
    peak_ind = 1:1:length(paramout_gauss(:,2)); 
    for i=1:1:length(peak_ind)
%         peak_ind(i) = find(abs(xin - paramout_gauss(i,2)) < tolerance, 1);  
        [~,peak_ind(i)] = min(abs(xin-paramout_gauss(i,2))); 
    end

    %% fit complex 2w data using gaussians found in above cell

    % IP = IP_full; 

    xin = E; 
    yin = data.'; 
    yin_abs = abs(yin); %./sum(abs(yin)); 
%     yin_phi = angle(yin); 
    yin_phi = mod(unwrap(angle(yin)),2*pi); 
%     yin = yin_abs; 
    yin = [yin_abs; yin_phi]; 

    % paramout_gauss(:,2) = peaks_guess'; 

    b_guess = yin(2, peak_ind); 
    c_guess = ones(size(b_guess)).*2; 
%     B_guess = ones(size(b_guess)).*0.5; 
    if slope==0
        guess = b_guess'; 
    else
        guess = [b_guess; c_guess]'; 
    end

    [paramout, fval] = fit2OmegaSum(xin, yin, paramout_gauss, guess, plotting); 
    % paramout(:,2) = mod(paramout(:,2), 2*pi); 

    % plotfun_fit(n, 810, xin, yin, fix, paramout, slope, peakflag)
    

end