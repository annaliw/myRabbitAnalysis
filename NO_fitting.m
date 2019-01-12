% import NO workspace
x1 = fliplr(E); 
y1 = fliplr(twoOmega_abs.')./fliplr((mean(abs(E_SpectraArray), 2)).'); 
% y1 = fliplr(twoOmega_abs.'); 
y3 = mod(fliplr(twoOmega_phi.'), 2*pi()); 
% y3 = fliplr(twoOmega_phi.'); 

%% expected peak positions 
n = 9:1:19; 
% IP = [(9.553+9.839+10.121)/3, 16.56, 15.667, 15.8, 15.9, 16.11, 16.26];% n, IP should already be
% defined in NO workspace
IP = [(9.553+9.839+10.121)/3, 16.56, 18.318, 21.722]; 
% w_slopes = [0.45, 0.0889, 0.03, 0.03, 0.03, 0.03, 0.03];
% w_slopes = [(0.33-0.18)/6 
% w_offsets = [

% IP = (9.553+9.839+10.121)/3; 

% widths = zeros(length(IP), length(n)); 
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
%     widths(i,:) = state_widths(i); 
    peaks(i,:) = n*(1240/810)-IP(i); 
end
peaks = peaks(:).'; 
% widths = widths(:).'; 

%% fit section
start = find(abs(x1-10.65)<0.02, 1); 
stop = find(abs(x1-20.16)<0.02, 1); 
test_x = x1(start:stop); 
test_yamp = y1(start:stop); 
test_ypha = y3(start:stop); 
test_ycom = test_yamp .* exp(1j*test_ypha);

% [val, peaks_guess] = findpeaks(test_yamp, test_x, 'MinPeakDistance', 0.1); 
peaks_guess = peaks; 
% widths_guess = widths; 
remove = [find(peaks_guess < test_x(1)) find(peaks_guess > test_x(end))]; 
peaks_guess(remove) = []; 
% widths_guess(remove) = []; 

peak_ind = 1:1:length(peaks_guess); 
for i=1:1:length(peaks_guess)
    peak_ind(i) = find(abs(test_x - peaks_guess(i)) < 0.03, 1);  
end

% form amplitude guess
amp_guess = test_yamp(peak_ind); 

% form peak width guess
sig_guess = ones(1, length(peaks_guess))*0.3; 

% form phase guess
pha_guess = test_ypha(peak_ind); 

% guess = [-amp_guess; peaks_guess; sig_guess; -pha_guess].'; 
guess = [amp_guess; sig_guess; pha_guess].'; 



%% fmincon
xin = test_x;  
yin_abs = test_yamp;
yin_phi = test_ypha; 
yin = [yin_abs ; yin_phi].'; 
% yin = [real(test_ycom) ; imag(test_ycom)].'; 

% lb = ones([1 length(guess(1,:))])*0.15e-03; 
% ub = ones([1 length(guess(1,:))])*0.5e-03;

% [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess, [], [], [], [], lb, ub); 
% [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess(1:(end-1),:)); 
% [paramout, fval] = fmincon(@(x) myfit(xin, yin, peaks_guess, x), guess); 
% fun_abs = @(x,xdata) abs(mydist(xdata, peaks_guess, x));
% fun_phi = @(x,xdata) angle(mydist(xdata, peaks_guess, x)); 
fun = @(guess,xdata) mydist(xdata, peaks_guess, guess);
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', ...
    'MaxFunctionEvaluations', 200000, 'MaxIterations', 10000);
% lb = zeros(size(guess));
% ub = zeros(size(guess)); 
% for i=1:1:length(peaks_guess)
%     lb(i,:) = [0, 0, -2*pi()]; 
%     ub(i,:) = [amp_guess(i)*1.2, amp_guess(i)*0.1, 2*pi()];
% end
lb = []; 
ub = []; 
[paramout, resnorm]=lsqcurvefit(fun,guess,xin,yin,lb,ub,options); 
paramout(:,3) = mod(paramout(:,3), 2*pi()); 

% x = x_abs.*exp(1j*x_phi); 

x_out = linspace(xin(1,1),xin(1,end),length(xin(1,:))*100);
y_out = mydist(x_out, peaks_guess, paramout); 

% figure; hold on; 
% yyaxis('left')
% scatter(xin, abs(yin), 'o'); 
% plot(x_out, abs(y_out), '--'); 
% plot(x_out, abs(mydist(x_out, guess)), '-'); 
% yyaxis('right')
% scatter(xin, angle(yin), '+'); 
% plot(x_out, angle(y_out)); 
% plot(x_out, angle(mydist(x_out, guess)), '-'); 
% hold off; 

%%

yfit_abs = y_out(:,1); 
yfit_phi = y_out(:,2); 
% yfit = y_out(:,1) + 1j*y_out(:,2); 
% yfit_abs = abs(yfit); 
% yfit_phi = -angle(yfit); 

fh = figure;
tmp = plot(xin, yin_abs); 
axl = AddHarmonicAxis(fh,IP,810);

% axl(1).XLabel.String = 'X';
% axl(2).XLabel.String = 'A';
% axl(3).XLabel.String = 'B';
% axl(4).XLabel.String = 'C';

hold on; 
plotfun(xin, yin_abs, mod(yin_phi, 2*pi()), x_out, yfit_abs, mod(yfit_phi, 2*pi()), peaks_guess, paramout); 

delete(tmp)

% xlim([xin(1), xin(end)]); 
% xlabel('photoelectron energy (eV)'); 

hold off; 

%%
figure; hold on; 
plotfun(xin, yin(:,1), yin(:,2), x_out, y_out(:,1), y_out(:,2), peaks_guess, paramout); 


%% window variation test

stpt = 400+5; 
edpt = 493-5; 
nsteps = 20; 
inc = 1; 

lslist = zeros([1 nsteps]); 
szlist = ((1:nsteps)-1)*2 + (edpt-stpt); 

% options = optimoptions('fmincon','MaxFunctionEvaluations', 100000, 'MaxIterations', 50001, 'OptimalityTolerance',1e-6);
% make window wider at each iteration
for i=1:nsteps
% i=1; 
    stpt = stpt - (i-1)*inc; 
    edpt = edpt + (i-1)*inc; 
    test_x = x1(stpt:edpt); 
    test_yamp = y1(stpt:edpt); 
    test_ypha = y3(stpt:edpt); 
    test_ycom = test_yamp .* exp(-1j*test_ypha); 
    
    xin = test_x; 
    yin = test_ycom; 

    peak_ind = 1:1:length(peaks_guess); 
    for j=1:1:length(peaks_guess)
        peak_ind(j) = find(abs(test_x - peaks_guess(j)) < 0.1, 1); 
    end

    % form amplitude guess
    amp_guess = test_yamp(peak_ind); 

    % form peak width guess
    sig_guess = ones(1, length(peaks_guess))*0.4; 

    % form phase guess
    pha_guess = test_ypha(peak_ind); 

    guess = [amp_guess; peaks_guess; sig_guess; -pha_guess].'; 
    add = 0; 
    
%     [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess,[],[],[],[],[],[],[],options); 
    [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess); 
    lslist(i) = fval; 
    
    figure;
    hold on; 
    yyaxis('left')
    scatter(xin, abs(yin), 'o'); 
    plot(x_out, abs(y_out), '-'); 
    for i=1:1:(length(peaks_guess + add))
        plot(x_out, abs(mydist(x_out, paramout(i,:))), '--'); 
    end
    yyaxis('right')
    scatter(xin, angle(yin), '+'); 
    plot(x_out, angle(y_out), '-'); 
    for i=1:1:(length(peaks_guess + add))
        plot(x_out, angle(mydist(x_out, paramout(i,:))), '--'); 
    end
    
    xlim([x_out(1), x_out(end)])
    hold off; 
end


% plist = (1:1:nsteps) .* start_guess(1,p)*(dist/2)*nsteps; 
figure; hold on; 
scatter(szlist, lslist); 
hold off; 

%% full spectrum line matching 
fh = figure; 
plot(x2, abs(y2)); 
IP = [13.778, 17.706, 18.077, 19.394]; 
axl = AddHarmonicAxis(fh,IP,810);

axl(1).XLabel.String = 'X';
axl(2).XLabel.String = 'A';
axl(3).XLabel.String = 'B';
axl(4).XLabel.String = 'C';

hold on; 
plot(x1, abs(y1), 'r'); 

hold off; 


%% plot phases

Xevens = [[-0.763840955470096] [-0.0981492795683128] [0.738121001907595] [1.82037310275875]]; 
Aevens = [[NaN] [-0.126313874443204] [0.483006605435046] [1.49339930039400]]; 
Bevens = [[NaN] [0.178167832356136] [0.696763479492293] [1.62675425154048]]; 
Cevens = [[NaN] [NaN] [1.19304163180125] [-4.16464371299721 + 2*pi]]; 

figure; hold on; 
scatter([12:2:18], Xevens, 'r', 'DisplayName', 'X');
scatter([12:2:18], Aevens, 'b', 'DisplayName', 'A');
scatter([12:2:18], Bevens, 'g', 'DisplayName', 'B');
scatter([12:2:18], Cevens, 'm', 'DisplayName', 'C');

xlabel('harmonic number')
ylabel('phase')
legend; 
hold off; 


%% Andrei's spectral integration
peaklist = [4.521 7.612 10.634 13.735 3.839 6.941 9.922 3.363 6.425 9.526]; 
window = 0.05; 

phase = zeros([1 length(peaklist)]); 
windowcheck = zeros([2 length(peaklist)]); 
for i=1:1:length(peaklist)
    idx_minus = find(abs(x1-peaklist(i)+window)<0.01,1); 
    idx_plus = find(abs(x1-peaklist(i)-window)<0.01,1); 
    windowcheck(:,i) = [idx_minus idx_plus]; 
    
    phase(i) = sum(y1(idx_minus:idx_plus).*exp(1j*y3(idx_minus:idx_plus))); 
end

Xphase_SI = angle(phase(1:4)); 
Aphase_SI = [NaN angle(phase(5:7))]; 
Bphase_SI = [NaN angle(phase(8:10))]; 

figure; hold on; 
scatter([12:2:18], -Xevens, 'r', 'o', 'DisplayName', 'X complex fit');
scatter([12:2:18], Xphase_SI, 'r', 'd', 'DisplayName', 'X spectral int');
scatter([12:2:18], -Aevens, 'b', 'o', 'DisplayName', 'A complex fit');
scatter([12:2:18], Aphase_SI, 'b', 'd', 'DisplayName', 'A spectral int');
scatter([12:2:18], -Bevens, 'g', 'o', 'DisplayName', 'B complex fit');
scatter([12:2:18], Bphase_SI, 'g', 'd', 'DisplayName', 'B spectral int');
xlabel('harmonic'); 
ylabel('phase'); 
legend; 
hold off; 


