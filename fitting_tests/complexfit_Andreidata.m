open('plotforAnna.fig'); 
h = gcf; 

%%
h = findobj(gca,'Type','line'); 
x=get(h,'Xdata'); 
y=get(h,'Ydata'); 

x1 = cell2mat(x(1)); 
x2 = cell2mat(x(2)); 
x3 = cell2mat(x(34)); 

y1 = cell2mat(y(1)); 
y2 = cell2mat(y(2)); 
y3 = cell2mat(y(34)); 

figure; hold on; 
plot(x1, y1, 'DisplayName', 'x1'); 
plot(x2, y2, 'DisplayName', 'x2'); 
yyaxis('right')
plot(x3, y3, 'DisplayName', 'x3'); 
legend; 

%% 331-389 triple peak 
stpt = 331; 
edpt = 389; 
test_x = x1(stpt:edpt); 
test_yamp = y1(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 

% four known peaks
guess = [2e-03 x2(360) 0.2 -2 ; 4e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2 ; 1.5e-03 8.1 0.15 -3]; 
ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1]; 
lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1];  
% double center peaks
% guess = [2e-03 x2(357) 0.2 -2 ; 2e-03 x2(364) 0.2 -2.2]; 
% ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1]; 
% lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1];  
% three peaks
% guess = [2e-03 x2(360) 0.2 -2 ; 4e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2]; 
% lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
% ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
% two peaks
% guess = [4e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2]; 
% lb = guess - guess.*[0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
% ub = guess + guess.*[0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1]; 

%% 210-257 triple peak
stpt = 210; 
edpt = 257; 
test_x = x1(stpt:edpt); 
test_yamp = y1(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 

% three peaks
guess = [3e-03 4.878 0.2 -2.5 ; 1.5e-03 5.18 0.2 2.5 ; 2e-03 5.36 0.2 -2.9]; 
lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1]; 

%% 400-493 triple peak
stpt = 400; 
edpt = 493; 
test_x = x1(stpt:edpt); 
test_yamp = y1(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 

% two peaks
guess = [2.5e-03 9.144 0.2 -2.9 ; 1.5e-03 9.5 0.2 2 ; 1.5e-03 9.656 0.1 1.5 ; 0.6e-03 9.767 0.2 1.4 ; 1e-03 10.61 0.2 0.7]; 
lb = guess - guess.*[1 0.01 0.1 -2 ; 1 0.01 0.1 2 ; 1 0.01 0.1 2 ; 1 0.01 0.1 2 ; 1 0.01 0.1 2];
ub = guess + guess.*[5 0.01 0.1 -2 ; 5 0.01 0.1 2 ; 5 0.01 0.1 2 ; 5 0.01 0.1 2 ; 5 0.01 0.1 2]; 


%% expected peak positions 
n = 9:1:19; 
IP = [13.778, 17.706, 18.077, 19.394]; 
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
    peaks(i,:) = n*(1240/810)-IP(i); 
end
peaks = peaks(:)'; 
peaks_guess = peaks; 
remove = [find(peaks_guess < test_x(1)) find(peaks_guess > test_x(end))]; 
peaks_guess(remove) = []; 

%% full data set fit

% test_x = x1; 
% test_yamp = y1; 
% test_ypha = y3; 
% test_ycom = test_yamp .* exp(-1j*test_ypha); 

% try splitting into high/low energy regions
% 0-3 eV
start = find(abs(x1-8.875)<0.02, 1); 
stop = find(abs(x1-10.95)<0.02, 1); 
test_x = x1(start:stop); 
test_yamp = y1(start:stop); 
test_ypha = y3(start:stop); 
test_ycom = test_yamp .* exp(-1j*test_ypha);

% [val, peaks_guess] = findpeaks(test_yamp, test_x, 'MinPeakDistance', 0.1); 


peak_ind = 1:1:length(peaks_guess); 
for i=1:1:length(peaks_guess)
    peak_ind(i) = find(abs(test_x - peaks_guess(i)) < 0.1, 1); 
end

% form amplitude guess
amp_guess = test_yamp(peak_ind); 

% form peak width guess
sig_guess = ones(1, length(peaks_guess))*0.4; 

% form phase guess
pha_guess = test_ypha(peak_ind); 

guess = [amp_guess; peaks_guess; sig_guess; -pha_guess].'; 
add = 0; 

% % add any extra peaks by hand
% guess = [guess ; 0.0025 8.5 0.2 0]; 
% add = 1; 


%% fmincon
xin = test_x; 
yin = test_ycom; 

lb = ones([1 length(guess(1,:))])*0.15e-03; 
ub = ones([1 length(guess(1,:))])*0.5e-03;

% [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess, [], [], [], [], lb, ub); 
% [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess(1:(end-1),:)); 
[paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess); 

x_out = linspace(xin(1),xin(end),length(xin)*100);
y_out = mydist(x_out, paramout); 

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

fh = figure;
tmp = plot(x_out, abs(y_out)); 
IP = [13.778, 17.706, 18.077, 19.394]; 
axl = AddHarmonicAxis(fh,IP,810);

axl(1).XLabel.String = 'X';
axl(2).XLabel.String = 'A';
axl(3).XLabel.String = 'B';
axl(4).XLabel.String = 'C';

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

delete(tmp); 
xlim([x_out(1), x_out(end)])
hold off; 


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



