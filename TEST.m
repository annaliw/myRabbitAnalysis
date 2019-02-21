%% fit delay integrated data for region
region = [1.744 2.7742]; 
n = 9:1:19; 

% form peaks guess
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
    peaks(i,:) = n*(1240/wavelength)-IP(i); 
end

x = E; 
y = sum(abs(E_SpectraArray),3);
y = squeeze(mean(y,2)); 
norm = mean(abs(E_SpectraArray),2); 
norm = sum(sum(norm, 3),1); 
y = y./norm; 
n=9:1:19; 
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
    peaks(i,:) = n*(1240/wavelength)-IP(i); 
end
peaks = peaks(:).'; 

% fit section set-up
start = find(abs(x-region(1))<0.05, 1); 
stop = find(abs(x-region(2))<0.05, 1); 
xin = x(start:stop); 
yin = y(start:stop).'; 

% find peaks and their indices
peaks_guess = peaks; 
remove = [find(peaks_guess < xin(1)) find(peaks_guess > xin(end))]; 
peaks_guess(remove) = []; 
peak_ind = 1:1:length(peaks_guess); 
for i=1:1:length(peaks_guess)
    peak_ind(i) = find(abs(xin - peaks_guess(i)) < 0.03, 1);  
end
amp_guess = yin(peak_ind); % form amplitude guess
sig_guess = ones(size(amp_guess))*0.1; % form width guess
guess = [amp_guess; peaks_guess; sig_guess].'; % full guess matrix

[paramout_gauss, fval_gauss] = fitGaussianSum(xin, yin, guess); 

%% fit complex 2w data using gaussians found in above cell
xin = x(start:stop); 
yin = twoOmega_signal(start:stop);
% yin = yin./squeeze(sum(sum(abs(E_SpectraArray(start:stop,:,:)),2),3)); 
yin = yin./sum(abs(twoOmega_signal(:))); 
yin = yin.'; 

% a_guess = abs(twoOmega_signal(peak_ind));
a_guess = ones([1 size(paramout_gauss,1)]); 
b_guess = angle(twoOmega_signal(peak_ind)).'; 
guess = [a_guess; b_guess].'; 

[paramout, fval] = fit2OmegaSum(xin, yin, paramout_gauss, guess); 
paramout(:,2) = mod(paramout(:,2), 2*pi); 

% plotfun_fit(n, 810, xin, yin, fix, paramout, slope, peakflag)






