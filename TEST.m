%% fit 2w abs data for region

% H2
region = [1.45 3.055]; % sideband 12
% region = [4.94 5.85]; % sideband 14
% region = [7.7548 9.1]; % sideband 16
% region = [10.8258 12.2]; % sideband 18

% CO2
% region = [2.9032 4.1419]; 
% region = [4.3742 5.7677]; 
% region = [5.8452 7.2387]; 
% region = [7.4 8.6]; 
% region = [8.8387 10.1548]; 
% region = [10.2 11]; 
% region = [11.8 15.5]; 
n = 9:1:19; 
tolerance = 0.02; 

% form peaks guess
peaks = zeros(length(IP), length(n)); 
for i=1:1:length(IP)
    peaks(i,:) = n*(1240/wavelength)-IP(i); 
end

signal = twoOmega_signal'; 
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
yin_abs = abs(signal)./mean(abs(signal(:))); 
yin_phi = mod(unwrap(angle(signal)),2*pi); 
yin = [yin_abs; yin_phi]; 
yin = yin(1,start:stop); 

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

[paramout_gauss, fval_gauss] = fitGaussianSum(xin, yin, guess); 
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
yin = yin(:,start:stop); 

% paramout_gauss(:,2) = peaks_guess'; 

% a_guess = abs(twoOmega_signal(peak_ind));
a_guess = ones([1 size(paramout_gauss,1)]); 
% a_guess = ((paramout_gauss(:,1)-abs(twoOmega_signal(peak_ind)))./paramout_gauss(:,1)).'; 
% b_guess = [yin(2, peak_ind), yin(2, peak_ind(end))]; 
b_guess = yin(2, peak_ind); 
c_guess = ones(size(a_guess)).*2; 
guess = [a_guess; b_guess; c_guess]'; 
% guess(end,3) = -1; 
% guess = [a_guess; b_guess].'; 

[paramout, fval] = fit2OmegaSum(xin, yin, paramout_gauss, guess); 
% paramout(:,2) = mod(paramout(:,2), 2*pi); 

% plotfun_fit(n, 810, xin, yin, fix, paramout, slope, peakflag)


%% CO2 results 
% harmonic
% peak center
% peak width
% phase
% phase slope

X_results = [11, 12, 13, 14, 15, 16, 17, 18, 19; ...
            2.9956, 4.5065, 6.0653, 7.5986, 9.0850, 10.5987, 12.1161, 13.6634, 15.2298; ...
            0.0695, 0.1077, 0.1516, 0.1407, 0.1955, 0.1905, 0.2409, 0.2386, 0.2916; ...
            4.4835, 1.3540, 4.1266, 0.6879, 3.4762, 6.1762, 2.7116, 4.9663, 2.0523; ...
            -0.0741, 0.7597, 0.6612, 0.5042, 0.3341, -0.3805, -0.9494, -0.6384, -0.5772];      
            
% B_results = [14, 15, 16, 17, 18; ...
%             3.8423, 5.3913, 6.9457, 8.4019, 9.9176; ...
%             0.1867, 0.2117, 0.2044, 0.2354, 0.2442; ...
%             0.7351, 3.5664, 0.0587, 2.9627, 5.3287; ...
%             -0.8187, -0.8329, -0.9433, -0.8979, 0.1308]; 
B_results = [14, 15, 16, 17, 18; ...
            3.8423, 5.3913, 6.9457, 8.4019, 9.9176; ...
            0.1867, 0.2117, 0.2044, 0.2354, 0.2442; ...
            0.7351, 3.5664, 0.0587, 2.9627, 5.3805; ...
            -0.8187, -0.8329, -0.9433, -0.8979, -0.4069]; 
            
A_results = [14, 15, 16, 17, 18; ...
            3.3142, 4.8688, 6.4449, 7.9802, 9.5040; ...
            0.0830, 0.0868, 0.0812, 0.0703, 0.0898; ...
            6.6779, 3.2225, 6.0585, 2.5194, 5.1312, ; ...
            -1.2448, -0.5061, -2.0478, -2.4706, -1.7024]; 
            
C_results = [15, 16, 17, 18, 19; ...
            3.5690, 5.0998, 6.6307, 8.1616, 9.6924; ...
            0.0500, 0.1000, 0.1000, 0.1000, 0.1000; ...
            -0.1395, 6.0643, 4.9288, 1.6014, 7.2075; ...
            -2.6401, -9.7198, 3.7609, 0.4263, -4.5575]; 
        
X_SI = [11, 12, 13, 14, 15, 16, 17, 18, 19; 4.4319, 1.3415, 4.0996, 0.6469, 3.4144, 6.1870, 2.7260, 4.9770, 2.0459]; 
B_SI = [14, 15, 16, 17, 18; 0.7148, 3.4663, 0.0389, 2.9907, 5.3571]; 
A_SI = [14, 15, 16, 17, 18; 0.4174, 3.2235, 6.0770, 0.6669, 4.9407]; 
        
% figure; hold on; 
% scatter(X_SI(1,2:2:end)+0.1, unwrap(X_SI(2,2:2:end)), 'd', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', ...
%     'DisplayName', 'X, spectral int.'); 
% scatter(B_SI(1,1:2:end)+0.1, unwrap(B_SI(2,1:2:end)), 'd', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
%     'DisplayName', 'A, spectral int.'); 
% scatter(A_SI(1,1:2:end)+0.1, [A_SI(2,1)+2*pi, A_SI(2,3:2:end)], 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', ...
%     'DisplayName', 'B, spectral int.'); 
% 
% scatter(X_results(1, 2:2:end), unwrap(X_results(4, 2:2:end)), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', ...
%     'DisplayName', 'X, 2-part fit'); 
% scatter(B_results(1, 1:2:end), unwrap(B_results(4, 1:2:end)), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
%     'DisplayName', 'A, 2-part fit'); 
% scatter(A_results(1, 1:2:end), unwrap(A_results(4, 1:2:end)), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', ...
%     'DisplayName', 'B, 2-part fit'); 

X_percent_diff = 100*abs((unwrap(X_SI(2,2:2:end)) - unwrap(X_results(4, 2:2:end)))./unwrap(X_SI(2,2:2:end))); 
B_percent_diff = 100*abs((unwrap(B_SI(2,1:2:end)) - unwrap(B_results(4, 1:2:end)))./unwrap(B_SI(2,1:2:end))); 
A_percent_diff = 100*abs(([A_SI(2,1)+2*pi, A_SI(2,3:2:end)] - unwrap(A_results(4, 1:2:end)))./[A_SI(2,1)+2*pi, A_SI(2,3:2:end)]); 

figure; hold on; 
scatter(X_SI(1,2:2:end), X_percent_diff, 'r', 'filled', 'DisplayName', 'X'); 
scatter(A_SI(1,1:2:end), A_percent_diff, 'b', 'filled', 'DisplayName', 'B'); 
scatter(B_SI(1,1:2:end), B_percent_diff, 'k', 'filled', 'DisplayName', 'A'); 
xlabel('sideband number'); ylabel('percent difference'); 
title('Percent difference for phase values (SI vs. 2-part fit)'); 
set(gca, 'FontSize', 14); 



%% H2 plotting
sidebands = 12:2:18; 
vstates = 0:1:5; 
include = 2:length(vstates); 

figure; hold on; 

% sideband 12
subplot(4,1,1); 
yyaxis left; 
plot(vstates(include), sb12_param(include,2), 'o-'); 
ylabel('phase'); 
yyaxis right; 
plot(vstates(include), sb12_param(include,3), 'o-'); 
ylabel('phase slope'); 
title('sideband 12'); 

% sideband 14
subplot(4,1,2); 
yyaxis left; 
plot(vstates(include), sb14_param(include,2), 'o-'); 
ylabel('phase'); 
yyaxis right; 
plot(vstates(include), sb14_param(include,3), 'o-'); 
ylabel('phase slope'); 
title('sideband 14'); 

% sideband 16
subplot(4,1,3); 
yyaxis left; 
plot(vstates(include), sb16_param(include,2), 'o-'); 
ylabel('phase'); 
yyaxis right; 
plot(vstates(include), sb16_param(include,3), 'o-'); 
ylabel('phase slope'); 
title('sideband 16'); 

% sideband 18
subplot(4,1,4); 
yyaxis left; 
plot(vstates(include), sb18_param(include,2), 'o-'); 
ylabel('phase'); 
yyaxis right; 
plot(vstates(include), sb18_param(include,3), 'o-'); 
ylabel('phase slope'); 
title('sideband 18'); 

xlabel('vibrational state'); 

hold off; 

%% sideband 12 phase plot
text_size = 14; 
phase_color = [1,0,0]; 
slope_color = [1,0.5,0]; 
line_weight = 1.5; 

figure; hold on; 

ax1 = gca; 
ax1_pos = ax1.Position; 
ax1.XColor = 'k'; 
ax1.XLabel.String = 'Vibrational State Number'; 
ax1.XTick = vstates(include);
ax1.YLabel.String = 'Phase (radians)';
ax1.FontSize = text_size; 
ax1.YColor = phase_color; 
ax1.LineWidth = line_weight*0.5; 
hold(ax1, 'on'); 
s1 = plot(vstates(include), sb12_param(include,2), 'o-'); 
s1.Color = phase_color; 
s1.LineWidth = line_weight; 

ax2 = axes('Position', ax1_pos, ...
    'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
    'Color', 'none', 'LineWidth', line_weight*0.5); 
ax2.XColor = 'none'; 
ax2.YColor = slope_color; 
ax2.YLabel.String = 'Phase Slope'; 
ax2.FontSize = text_size; 
hold(ax2, 'on'); 
s2 = plot(vstates(include), sb12_param(include,3), 'o-'); 
s2.Color = slope_color; 
s2.LineWidth = line_weight; 

% text('sideband 12'); 

%% XUV only subtract

sub_signal = mean(mean(tmp1,2),3)-XUV_only./sum(XUV_only(:));
plotfun_rabbitspectrum(9:1:19, wavelength, E, sub_signal, 'average'); 





