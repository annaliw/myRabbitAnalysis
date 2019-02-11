%% import H2 data
clear all
load('HydrogenSpectra.mat'); 
% load('ArgonSpectra.mat'); 
%%
wavelength = 810; 
slope = 0; 
x1 = E; 
% y1 = (twoOmega_abs.')./((mean(abs(E_SpectraArray), 2)).'); 
y1 = abs(restackedSpectra); 
y3 = mod(angle(restackedSpectra), 2*pi); 

%% expected peak positions 
n=9:1:19; 
IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206 16.7461]; 
IP_label = ["0", "1", "2", "3", "4", "5", "5.5"]; 
% IP = [15.76 27.63]; 
% IP_label = ["S0", "P3/2"]; 

% widths = zeros(length(IP), length(n)); 
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
%     widths(i,:) = state_widths(i); 
    peaks(i,:) = n*(1240/wavelength)-IP(i); 
end
peaks = peaks(:).'; 
% widths = widths(:).'; 
plotfun_rabbitspectrum(n, IP, IP_label, wavelength, E, restackedSpectra, 'twoOmega')
%% fit section
start = find(abs(x1-in1(1))<0.05, 1); 
stop = find(abs(x1-in2(1))<0.05, 1); 
test_x = x1(start:stop); 
% test_yamp = y1(start:stop)/sum(y1(start:stop));
test_yamp = y1(start:stop); 
test_ypha = y3(start:stop); 
test_ycom = test_yamp .* exp(1j*test_ypha);

peaks_guess = peaks; 
remove = [find(peaks_guess < test_x(1)) find(peaks_guess > test_x(end))]; 
peaks_guess(remove) = []; 

peak_ind = 1:1:length(peaks_guess); 
for i=1:1:length(peaks_guess)
    peak_ind(i) = find(abs(test_x - peaks_guess(i)) < 0.03, 1);  
end

% form amplitude guess
amp_guess = test_yamp(peak_ind); 

% form peak width guess
sig_guess = ones(1, length(peaks_guess))*0.1; 

% form phase guess
pha_guess = test_ypha(peak_ind); 

% form slope guess (optional)
slope_guess = zeros(size(peak_ind)) + 1; 
%     slope_guess(end) = -0.4; 
guess_fixwidth = [amp_guess; peaks_guess; pha_guess; slope_guess].'; 
guess_fixpeaks = [amp_guess; sig_guess; pha_guess; slope_guess].'; 

width=0.1; 
xin = test_x;  
yin_abs = test_yamp; 
yin_phi = unwrap(test_ypha); 
% yin_phi = test_ypha; 
yin = [yin_abs ; yin_phi].'; 


%% fit using fixed peaks to get width estimate

slope = 1; 
variance = 0.02; 
meanValue = 0.01; 
guess = guess_fixpeaks; 
lb = [guess(:,1)*0, guess(:,2)-0.1, guess(:,3)-0.1, guess(:,4)*0]; 
ub = [guess(:,1)*10, guess(:,2)+0.1, guess(:,3)+0.1, guess(:,4)*10]; 
if slope==0
    guess = guess_fixwidth(:,1:3); 
    lb = lb(:,1:3); 
    ub = ub(:,1:3); 
end
guess = guess + sqrt(variance)*randn(size(guess));
 
fun = @(guess,xdata) mydist_fixpeaks(xdata, peaks_guess, guess, slope);
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Display', 'off', ...
    'MaxFunctionEvaluations', 200000, 'MaxIterations', 10000);
% lb = []; 
% ub = []; 
[paramout_fixpeaks, resnorm]=lsqcurvefit(fun,guess,xin,yin,lb,ub,options); 
% paramout_fixpeaks(:,3) = mod(paramout_fixpeaks(:,3), 2*pi); 

plotfun_fit(n, IP, IP_label, wavelength, xin, yin, peaks_guess, paramout_fixpeaks, slope, 1); 

%% fit using fixed width and variable peak
% width = mean(paramout_fixpeaks(1:end,2)); 
% yin(:,2) = unwrap(yin(:,2)); 
width = 0.08; 

slope = 1; 
variance = 0.02; 
meanValue = 0.01; 
guess = guess_fixwidth; 
lb = [guess(:,1)*0, guess(:,2)-0.1, guess(:,3)-0.1, -guess(:,4)*0]; 
ub = [guess(:,1)*5, guess(:,2)+0.1, guess(:,3)+0.1, guess(:,4)*10]; 
if slope==0
    guess = guess_fixwidth(:,1:3); 
    lb = lb(:,1:3); 
    ub = ub(:,1:3); 
end
% guess = guess + sqrt(variance)*randn(size(guess));
guess(:,3) = mod(guess(:,3), 2*pi); 
% guess(:,4) = guess(:,4)*20; 
 
fun = @(guess,xdata) mydist_fixwidth(xdata, width, guess, slope);
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Display', 'off', ...
    'MaxFunctionEvaluations', 200000, 'MaxIterations', 10000);

% lb = []; 
% ub = []; 
[paramout_fixwidth, resnorm]=lsqcurvefit(fun,guess,xin,yin,lb,ub,options); 
paramout_fixwidth(:,3) = mod(paramout_fixwidth(:,3), 2*pi); 
% paramout_fixwidth(:,4) = paramout_fixwidth(1,4); 

plotfun_fit(n, IP, IP_label, wavelength, xin, yin, width, paramout_fixwidth, slope, 0); 
[paramout_fixwidth(:,3) paramout_fixwidth(:,4)]

%%
width_list = mean(paramout_fixpeaks(1:end,2))*(0:0.05:1)*4; 
samples = 5; 
slope = 1; 
meanvar = zeros([2 size(guess_fixpeaks) length(width_list)]); 
resnorm_ii = zeros([2 size(width_list)]); 
resnorm_jj = zeros([1 samples]); 
tmp = zeros([size(guess_fixpeaks) samples]); 

for ii=1:1:length(width_list)
    width = width_list(ii); 
    for jj=1:1:samples
        guess = guess_fixwidth; 
        lb = [guess(:,1)*0, guess(:,2)-0.1, guess(:,3)-0.1, guess(:,4)*0]; 
        ub = [guess(:,1)*5, guess(:,2)+0.1, guess(:,3)+0.1, guess(:,4)*10]; 
        if slope==0
            guess = guess_fixwidth(:,1:3); 
            lb = lb(:,1:3); 
            ub = ub(:,1:3); 
        end
        guess(:,2) = guess(:,2) + sqrt(variance)*randn(size(guess(:,2)));
        guess(:,3) = mod(-guess(:,3), 2*pi); 
        % guess(:,4) = guess(:,4); 

        fun = @(guess,xdata) mydist_fixwidth(xdata, width, guess, slope);
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', ...
            'MaxFunctionEvaluations', 200000, 'MaxIterations', 10000);
        [paramout_fixwidth, resnorm]=lsqcurvefit(fun,guess,xin,yin,lb,ub,options); 
        paramout_fixwidth(:,3) = mod(paramout_fixwidth(:,3), 2*pi); 
        paramout_fixwidth(:,4) = paramout_fixwidth(1,4); 
        tmp(:,:,jj) = paramout_fixwidth;
        resnorm_jj(jj) = resnorm; 
    end
    
    meanvar(1, :, :, ii) = mean(tmp, 3); 
    meanvar(2, :, :, ii) = var(tmp, 0, 3); 
    resnorm_ii(1, ii) = mean(resnorm_jj); 
    resnorm_ii(2, ii) = var(resnorm_jj); 

end

figure; hold on; 
scatter(width_list, resnorm_ii(1,:)); 
hold off; 
%%
peakno = 2; 

figure; hold on; 
% plot peak locations vs width
errorbar(width_list, squeeze(meanvar(1, peakno, 2, :)), squeeze(meanvar(2, peakno, 2, :)), 'o', 'DisplayName', 'peak center'); 
line(width_list, ones(size(width_list))*guess_fixwidth(peakno, 2))
xlabel('peak width (eV)'); 
ylabel('peak center (eV)'); 
yyaxis right; 
plot(width_list, squeeze(meanvar(2, peakno, 2, :)), 'o-')
% legend; 
hold off; 
figure; hold on; 
% plot peak locations vs width
errorbar(width_list, squeeze(meanvar(1, peakno, 3, :)), squeeze(meanvar(2, peakno, 3, :)), 'o', 'DisplayName', 'phase'); 
line(width_list, ones(size(width_list))*guess_fixwidth(peakno, 3))
xlabel('peak width (eV)'); 
ylabel('phase (radians)'); 
yyaxis right; 
plot(width_list, squeeze(meanvar(2, peakno, 3, :)), 'o-')
% legend; 
hold off; 
figure; hold on; 
% plot peak locations vs width
errorbar(width_list, squeeze(meanvar(1, peakno, 4, :)), squeeze(meanvar(2, peakno, 4, :)), 'o', 'DisplayName', 'phase slope'); 
line(width_list, ones(size(width_list))*guess_fixwidth(peakno, 4))
xlabel('peak width (eV)'); 
ylabel('phase slope (radians/eV)'); 
yyaxis right; 
plot(width_list, squeeze(meanvar(2, peakno, 4, :)), 'o-')
% legend; 
hold off; 




