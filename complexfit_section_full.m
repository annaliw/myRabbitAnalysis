function [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E, data, fitRegion, plotting)
    global IP; 
    global IP_label; 
    %% fit 2w abs data for region
    
    region = fitRegion; 
    n = 9:1:19; 
    tolerance = 0.02; 

    % form peaks guess
    peaks = zeros(length(IP), length(n)); 
    for i=1:1:length(IP)
        peaks(i,:) = n*(1240/wavelength)-IP(i); 
    end

    signal = data'; 
    x = E; 
    y = signal; 
    peaks = zeros(length(IP), length(n)); 
    for i=1:1:length(IP)
        peaks(i,:) = n*(1240/wavelength)-IP(i); 
    end
    peaks = peaks(:).'; 

    % fit section set-up
    start = find(abs(x-region(1))<tolerance, 1, 'last'); 
    stop = find(abs(x-region(2))<tolerance, 1, 'first'); 

    xin = x(start:stop); 
    yin = y(start:stop); 
    yin_abs = abs(yin)./sum(abs(yin)); 
    yin_phi = mod(unwrap(angle(yin)),2*pi); 
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
    guess = [amp_guess; peaks_guess; sig_guess].'; % full guess matrix

    [paramout_gauss, fval_gauss] = fitGaussianSum(xin, yin, guess, plotting); 
    % find new peak center indices
    % paramout_gauss(:,2) = peaks_guess'; 
    peak_ind = 1:1:length(paramout_gauss(:,2)); 
    for i=1:1:length(peak_ind)
        peak_ind(i) = find(abs(xin - paramout_gauss(i,2)) < tolerance, 1);  
    end

    %% fit complex 2w data using gaussians found in above cell

    % IP = IP_full; 

    xin = x(start:stop); 
    yin = [yin_abs; yin_phi]; 

    % paramout_gauss(:,2) = peaks_guess'; 

    % a_guess = abs(twoOmega_signal(peak_ind));
    a_guess = ones([1 size(paramout_gauss,1)]); 
    % a_guess = ((paramout_gauss(:,1)-abs(twoOmega_signal(peak_ind)))./paramout_gauss(:,1)).'; 
    % b_guess = [yin(2, peak_ind), yin(2, peak_ind(end))]; 
    b_guess = yin(2, peak_ind); 
    c_guess = ones(size(a_guess)).*2; 
    guess = [a_guess; b_guess; c_guess]'; 
    % guess(end,3) = -1; 
%     guess = [a_guess; b_guess].'; 

    [paramout, fval] = fit2OmegaSum(xin, yin, paramout_gauss, guess, plotting); 
    % paramout(:,2) = mod(paramout(:,2), 2*pi); 

    % plotfun_fit(n, 810, xin, yin, fix, paramout, slope, peakflag)
    

end