%% import H2 data
clear all
% load('ArgonSpectra.mat'); 
load('HydrogenSpectra.mat'); 
wavelength = 810; 
slope = 0; 

x1 = E; 
% y1 = (twoOmega_abs.')./((mean(abs(E_SpectraArray), 2)).'); 
y1 = abs(restackedSpectra); 
y3 = mod(angle(restackedSpectra), 2*pi); 

n=9:1:19; 
% IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206 16.7461]; 
% IP_label = ["0", "1", "2", "3", "4", "5", "5.5"]; 
IP = [15.38174 15.65097 15.90469 16.16865 16.39351 16.62206]; 
IP_label = ["0", "1", "2", "3", "4", "5"]; 
% IP = [15.76 27.63]; 
% IP_label = ["S0", "P3/2"]; 

% widths = zeros(length(IP), length(n)); 
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
%     widths(i,:) = state_widths(i); 
    peaks(i,:) = n*(1240/wavelength)-IP(i); 
end
peaks = peaks(:).'; 
plotfun_rabbitspectrum(n, IP, IP_label, wavelength, E, restackedSpectra, 'twoOmega'); 
%% fit section
start = find(abs(x1-in1(1))<0.05, 1); 
stop = find(abs(x1-in2(1))<0.05, 1); 
test_x = x1(start:stop); 
% test_yamp = y1(start:stop)/sum(y1(start:stop));
test_yamp = y1(start:stop); 
test_ypha = unwrap(y3(start:stop)); 
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
slope_guess = zeros(size(peak_ind)) + Ar_phase(3,4); 
%     slope_guess(end) = -0.4; 
guess_fixwidth = [amp_guess; peaks_guess; pha_guess; slope_guess].'; 
guess_fixpeaks = [amp_guess; sig_guess; pha_guess; slope_guess].'; 

width=0.1; 
xin = test_x;  
yin_abs = test_yamp; 
yin_phi = unwrap(test_ypha); 
% yin_phi = test_ypha; 
yin = [yin_abs ; yin_phi].'; 

%% fit using fixed width and variable peak
% width = mean(paramout_fixpeaks(1:end,2)); 
% yin(:,2) = unwrap(yin(:,2)); 
width = 0.11; 

slope = 1; 
variance = 0.2; 
meanValue = 0.01; 
guess = guess_fixwidth; 
lb = [guess(:,1)*0, guess(:,2)-0.1, guess(:,3)-0.5, -guess(:,4)*2]; 
ub = [guess(:,1)*5, guess(:,2)+0.1, guess(:,3)+0.5, guess(:,4)*2]; 
if slope==0
    guess = guess_fixwidth(:,1:3); 
    lb = lb(:,1:3); 
    ub = ub(:,1:3); 
end
% guess = guess + sqrt(variance)*randn(size(guess));
% guess(:,3) = mod(guess(:,3), 2*pi); 
 
% fun = @(guess,xdata) mydist_fixwidth(xdata, width, guess, slope);
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', ...
%     'MaxFunctionEvaluations', 200000, 'MaxIterations', 10000); 
fun = @(guess, data) complex_lsq(data, width, guess, slope, peakflag); 
data = [xin; yin_abs.*exp(1j*yin_phi)]; 
options = optimset('MaxFunEvals', 200000, 'MaxIter', 100000); 

lb = []; 
ub = []; 
% [paramout_fixwidth, resnorm, residual, exitflag, output]=lsqcurvefit(fun,guess,xin,yin,lb,ub,options); 
[paramout_fixwidth, fval, exitflag, output] = fminsearch(@(x) complex_lsq(data, width, x, slope, 0), guess, options); 
paramout_fixwidth(:,3) = mod(paramout_fixwidth(:,3), 2*pi); 
% paramout_fixwidth(:,4) = paramout_fixwidth(1,4); 

plotfun_fit(n, IP, IP_label, wavelength, xin, yin, width, paramout_fixwidth, slope, 0); 
[paramout_fixwidth(:,3) paramout_fixwidth(:,4)]

%% save Argon results, compare to H2 results
Ar_phase = [12 14 16 18; 5.1850 5.8276 0.3966 1.4859; 1.0978 1.037 0.9039 0.4485]; 
% H2_12_phase = [1 2 3 4 5; 4.7957 4.8650 4.6856 4.7727 4.6740; 0.0074 1.3941 0.0064 0.0210 0.7634]; 
% H2_14_phase = [1 2 3 4 5; 5.4877 5.4769 5.4476 5.3375 5.3559; 1.6247 0.1434 -0.9406 -1.0431 0.2741]; 
% H2_16_phase = [1 2 3 4 5; 6.2092 6.1293 6.2322 0.0076 6.2271; 0.1813 1.5851 2.4411 0.6244 -2.7230]; 
% H2_18_phase = [1 2 3 4 5; 1.2008 1.0910 1.1045 1.0012 1.1391; 1.1231 0.3029 -0.2744 -0.0902 -0.3323]; 
H2_12_phase = [1 2 3 4 5; 4.8853 4.7794 4.6710 4.8169 4.6888; 3.5050 1.1076 1.1484 1.7174 1.2690]; 
H2_14_phase = [1 2 3 4 5; 5.2711 5.2488 5.6277 5.5031 5.4515; 2.6573 6.2246 6.2737 3.0227 1.6551]; 
H2_16_phase = [1 2 3 4 5; 6.1560 6.1905 6.2043 6.2292 6.2342; 0.7868 1.7929 1.5691 0.8916 -1.4538]; 
H2_18_phase = [1 2 3 4 5; 1.2024 1.0990 1.1144 1.0048 1.1342; 0.9216 -0.2096 -0.8959 -1.1378 -1.2672]; 

figure; hold on; 
scatter(H2_12_phase(1,:), mod(H2_12_phase(2,:)-Ar_phase(2,1), 2*pi), 'MarkerEdgeColor', 'b', 'DisplayName', 'sideband 12'); 
plot(H2_12_phase(1,:), mod(ones(size(H2_12_phase(1,:)))*Ar_phase(2,1), 2*pi), 'Color', 'b', 'HandleVisibility', 'off'); 
scatter(H2_14_phase(1,:), mod(H2_14_phase(2,:)-Ar_phase(2,2), 2*pi), 'MarkerEdgeColor', 'g', 'DisplayName', 'sideband 14'); 
plot(H2_12_phase(1,:), mod(ones(size(H2_12_phase(1,:)))*Ar_phase(2,2), 2*pi), 'Color', 'g', 'HandleVisibility', 'off'); 
scatter(H2_16_phase(1,:), mod(H2_16_phase(2,:)-Ar_phase(2,3), 2*pi), 'MarkerEdgeColor', 'r', 'DisplayName', 'sideband 16'); 
plot(H2_12_phase(1,:), mod(ones(size(H2_12_phase(1,:)))*Ar_phase(2,3), 2*pi), 'Color', 'r', 'HandleVisibility', 'off'); 
scatter(H2_18_phase(1,:), mod(H2_18_phase(2,:)-Ar_phase(2,4), 2*pi), 'MarkerEdgeColor', 'k', 'DisplayName', 'sideband 18'); 
plot(H2_12_phase(1,:), mod(ones(size(H2_12_phase(1,:)))*Ar_phase(2,4), 2*pi), 'Color', 'k', 'HandleVisibility', 'off'); 
xlabel('H_2 vibrational state'); 
ylabel('phase (radians mod 2\pi)'); 
title('H_2 phase referenced to Ar')
legend; 
hold off; 

figure; hold on; 
scatter(H2_12_phase(1,:), H2_12_phase(3,:)-Ar_phase(3,1), 'MarkerEdgeColor', 'b', 'DisplayName', 'sideband 12'); 
plot(H2_12_phase(1,:), ones(size(H2_12_phase(1,:)))*Ar_phase(3,1), 'Color', 'b', 'HandleVisibility', 'off'); 
scatter(H2_14_phase(1,:), H2_14_phase(3,:)-Ar_phase(3,2), 'MarkerEdgeColor', 'g', 'DisplayName', 'sideband 14'); 
plot(H2_14_phase(1,:), ones(size(H2_12_phase(1,:)))*Ar_phase(3,2), 'Color', 'g', 'HandleVisibility', 'off'); 
scatter(H2_16_phase(1,:), H2_16_phase(3,:)-Ar_phase(3,3), 'MarkerEdgeColor', 'r', 'DisplayName', 'sideband 16'); 
plot(H2_16_phase(1,:), ones(size(H2_12_phase(1,:)))*Ar_phase(3,3), 'Color', 'r', 'HandleVisibility', 'off'); 
scatter(H2_18_phase(1,:), H2_18_phase(3,:)-Ar_phase(3,4), 'MarkerEdgeColor', 'k', 'DisplayName', 'sideband 18'); 
plot(H2_18_phase(1,:), ones(size(H2_12_phase(1,:)))*Ar_phase(3,4), 'Color', 'k', 'HandleVisibility', 'off'); 
xlabel('H_2 vibrational state'); 
ylabel('phase slope (unitless)'); 
title('H_2 phase slope referenced to Ar')
legend; 
hold off; 
%%
% width_list = mean(paramout_fixpeaks(1:end,2))*(0:0.05:1)*4; 
width_list = 0.1*(0:0.05:1)*4; 
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

%% data fluctuation for error bars

samples = 10; 

width = 0.11; 
slope = 1; 
guess = guess_fixwidth; 
fun = @(guess, data) complex_lsq(data, width, guess, slope, peakflag); 
options = optimset('MaxFunEvals', 200000, 'MaxIter', 100000); 

param_array = zeros([size(guess), samples]); 

for ii=1:1:samples
    data = [xin; yin_abs.*exp(1j*yin_phi)]; 
    [paramout_fixwidth, fval, exitflag, output] = fminsearch(@(x) complex_lsq(data, width, x, slope, 0), guess, options); 
    paramout_fixwidth(:,3) = mod(paramout_fixwidth(:,3), 2*pi); 
    
    param_array(:,:,ii) = paramout_fixwidth; 
end


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
