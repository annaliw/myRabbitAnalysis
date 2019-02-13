function [paramout_fixwidth, fval, exitflag, output] = complexfit_section(IP, E, data, fitRegion, slope_guess, slopeflag, peakflag)
    % ONLY FIXED WIDTH VERSION OF FITTING, fminsearch on lsq
    if ~exist('peakflag', 'var')
        peakflag = 0; 
    end
    if ~exist('slopeflag', 'var')
        slopeflag = 0; 
    end
    
    % full data set-up
    x = E; 
    y_abs = abs(data); 
    y_phi = mod(angle(data), 2*pi); 
    n=9:1:19; 
    peaks = zeros(length(IP), length(n)); 
    for i=1:1:length(IP)
        peaks(i,:) = n*(1240/wavelength)-IP(i); 
    end
    peaks = peaks(:).'; 
    
    % fit section set-up
    start = find(abs(x-fitRegion(1))<0.05, 1); 
    stop = find(abs(x1-fitRegion(2))<0.05, 1); 
    xin = x(start:stop); 
    yin_abs = y_abs(start:stop); 
    yin_phi = unwrap(y_phi(start:stop)); 
    yin = [yin_abs ; yin_phi].'; 

    % find peaks and their indices
    peaks_guess = peaks; 
    remove = [find(peaks_guess < test_x(1)) find(peaks_guess > test_x(end))]; 
    peaks_guess(remove) = []; 
    peak_ind = 1:1:length(peaks_guess); 
    for i=1:1:length(peaks_guess)
        peak_ind(i) = find(abs(test_x - peaks_guess(i)) < 0.03, 1);  
    end
    amp_guess = yin_abs(peak_ind); % form amplitude guess
    phi_guess = yin_phi(peak_ind); % form phase guess
    slope_guess = zeros(size(peak_ind)) + slope_guess; % form slope guess
    guess_fixwidth = [amp_guess; peaks_guess; phi_guess; slope_guess].'; % full guess matrix
    if slopeflag==0
        guess = guess_fixwidth(:,1:3); 
    end
    
    % do fit
    width = fitRegion(3); 
    fun = @(guess, data) complex_lsq(data, width, guess, slope, peakflag); 
    data = [xin; yin_abs.*exp(1j*yin_phi)]; 
    options = optimset('MaxFunEvals', 200000, 'MaxIter', 100000); 
    [paramout_fixwidth, fval, exitflag, output] = fminsearch(@(x) complex_lsq(data, width, x, slope, 0), guess, options); 
    paramout_fixwidth(:,3) = mod(paramout_fixwidth(:,3), 2*pi); 
    
    
    
    

end