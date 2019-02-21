%% fit delay integrated data for region
region = [1.744 2.7742]; 
n = 9:1:19; 

% form peaks guess
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
    peaks(i,:) = n*(1240/wavelength)-IP(i); 
end

x = E; 
y = squeeze(sum(sum(abs(E_SpectraArray),2),3));
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
yin = twoOmega_signal(start:stop).'; 

a_guess = abs(twoOmega_signal(peak_ind)); 
b_guess = angle(twoOmega_signal(peak_ind)); 
guess = [a_guess; b_guess].'; 

[paramout, fval] = fit2OmegaSum(xin, yin, paramout_gauss, guess); 






