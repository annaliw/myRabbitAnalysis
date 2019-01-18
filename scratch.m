%%% rough extraction of b3Pi and X phases
% first import NO_analysis.mat

wavelength=810; 

harmonics = 12:1:19; 
b3Pi_peaks =   [find(abs(E-1.932)<0.001) find(abs(E-1.658)<0.001); 
                find(abs(E-3.481)<0.001) find(abs(E-3.263)<0.001);
                find(abs(E-4.995)<0.001) find(abs(E-4.762)<0.001);
                find(abs(E-6.526)<0.005) find(abs(E-6.264)<0.001);
                find(abs(E-8.081)<0.001) find(abs(E-7.76)<0.001); 
                find(abs(E-9.526)<0.001) find(abs(E-9.332)<0.001); 
                find(abs(E-11.1)<0.005) find(abs(E-10.75)<0.005);
                find(abs(E-12.6)<0.005) find(abs(E-12.35)<0.005)]; 
X_peaks =      [find(abs(E-8.283)<0.001) find(abs(E-9.092)<0.001); 
                find(abs(E-10.56)<0.005) find(abs(E-10.11)<0.005);
                find(abs(E-11.91)<0.005) find(abs(E-11.46)<0.005);
                find(abs(E-13.4)<0.005) find(abs(E-12.91)<0.005); 
                find(abs(E-15.13)<0.005) find(abs(E-14.33)<0.005); 
                find(abs(E-16.62)<0.005) find(abs(E-15.88)<0.005); 
                find(abs(E-18.19)<0.005) find(abs(E-17.48)<0.005); 
                find(abs(E-19.66)<0.005) find(abs(E-19.02)<0.005)]; 
            
A1Pi_peaks =   []; 

% twoOmega_abs_norm = twoOmega_abs./sum(twoOmega_abs); 

b3Pi_phase = zeros([1 length(b3Pi_peaks)]); 
X_phase = zeros([1 length(X_peaks)]); 
% for i=1:length(b3Pi_peaks)
%     b3Pi_phase(i) = sum(twoOmega_abs(b3Pi_peaks(i,1):b3Pi_peaks(i,2)).*...
%     twoOmega_phi(b3Pi_peaks(i,1):b3Pi_peaks(i,2)))/...
%     sum(twoOmega_abs(b3Pi_peaks(i,1):b3Pi_peaks(i,2))); 
% 
%     X_phase(i) = sum(twoOmega_abs(X_peaks(i,1):X_peaks(i,2)).*...
%     twoOmega_phi(X_peaks(i,1):X_peaks(i,2)))/...
%     sum(twoOmega_abs(X_peaks(i,1):X_peaks(i,2))); 
% end
for i=1:length(b3Pi_peaks)
    b3Pi_phase(i) = sum(twoOmega_abs(b3Pi_peaks(i,1):b3Pi_peaks(i,2)).*...
    exp(1j*twoOmega_phi(b3Pi_peaks(i,1):b3Pi_peaks(i,2)))); 

    X_phase(i) = sum(twoOmega_abs(X_peaks(i,1):X_peaks(i,2)).*...
    exp(1j*twoOmega_phi(X_peaks(i,1):X_peaks(i,2)))); 
end

figure; hold on; 
scatter(harmonics(1:2:7), angle(b3Pi_phase(1:2:7))); 
scatter(harmonics(1:2:7), angle(X_phase(1:2:7))); 

%% compare to fit extraction NO

window = 3; 
peaks = peaks(:); 
phase_si = spectral_integration(x1, y1, y3, peaks, window); 
phase_si = reshape(phase_si, [length(IP), length(n)]); 
peaks = reshape(peaks, [length(IP), length(n)]); 

% % unbinned energy conversion
% X_fit = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19; 175.9033, 43.9565, 44.1536, 11.8212, 373.9577, 373.9577, 373.9577, 373.9577, 373.9577, 373.9577;...
%     4.5433, 4.8752, 0.3235, 3.8694, 0.9177, 4.7954, 1.9589, 5.7624, 3.1644, 5.4831]; 
% b_fit = [12, 13, 14, 15, 16, 17, 18, 19; 713.9706, 91.7449, 53.8667, 43.9565, 44.1536, 11.8212, 373.9577, 373.9577; ...
%     6.2832, 4.0629, 1.2348, 2.3675, 1.6701, 5.7794, 3.1728, 0.5818]; % harmonic, resnorm, phase in columns
% A_fit = [14, 16, 17, 18, 19; 91.7449, 29.2274, 44.1536, 11.8212, 373.9577; ...
%     2.0391, 1.8822, 2.0423, 4.9629, 3.3323]; % 16 ?
% c_fit = [16, 18; 91.7449, 29.2274; 1.0662, 2.5632]; % 18 ?

% binned energy conversion
X_fit = [10, 11, 12, 13, 14, 15, ; 5.3907 1.0682e+03, 152.4920, 3.0539, 1.8723e+03, 1.8723e+03, ; ...
    2.3774, 5.4546, 3.5500, 6.1692, 4.3168, 1.4265, 3.4755, ]; 
b_fit = [12, 14, 15, 16, 17, 18, ; 5.5969, 11.7024, 5.3907, 152.4920, 3.0539, 2.9140, ; ...
    3.0342, 4.2788, 1.4772, 4.9925, 2.3021, 6.1194, ]; % harmonic, resnorm, phase in columns
A_fit = [13, 15, 16, 17, 18, 19; 5.5969, 11.7024, 5.3907, 152.4920, 3.0539, 2.9140; ...
    0.8636, 1.6607, 5.1315, 4.2886, 1.4322, 0.7538]; 
c_fit = [15, 17, 18, 19; 5.5969, 11.7024, 5.3907, 152.4920; 1.4832, 1.3597, 5.5888, 1.7090]; 

    
omega = 3*10^8/(810*10^(-9)); 

figure; hold on; 
scatter(X_fit(1,1:2:end), mod(X_fit(3,1:2:end), 2*pi()), 'r', '+', 'LineWidth',1.5); 
scatter(X_fit(1,1:2:end), mod(phase_si(1,2:2:end), 2*pi()), 'r', 'o', 'LineWidth',1.5); 
scatter(b_fit(1,1:2:end), mod(b_fit(3,1:2:end), 2*pi()), 'b', '+', 'LineWidth',1.5); 
scatter(b_fit(1,1:2:end), mod(phase_si(2,4:2:end), 2*pi()), 'b', 'o', 'LineWidth',1.5); 
scatter(A_fit(1,[1,2,4]), mod(A_fit(3,[1,2,4]), 2*pi()), 'g', '+', 'LineWidth',1.5); 
scatter(A_fit(1,[1,2,4]), mod(phase_si(3,6:2:end), 2*pi()), 'g', 'o', 'LineWidth',1.5); 
scatter(c_fit(1,:), mod(c_fit(3,:), 2*pi()), 'k', '+', 'LineWidth',1.5); 
scatter(c_fit(1,:), mod(phase_si(4,10:2:end), 2*pi()), 'k', 'o', 'LineWidth',1.5); 
xlabel('harmonic number'); 
% scatter(peaks(1,2:2:end), mod(X_fit(3,1:2:end), 2*pi()), 'r', '+', 'LineWidth',2); 
% scatter(peaks(2,4:2:end), mod(b_fit(3,1:2:end), 2*pi()), 'b', '+', 'LineWidth',2); 
% scatter(peaks(3,6:2:end), mod(A_fit(3,[1,2,4]), 2*pi()), 'g', '+', 'LineWidth',2); 
% scatter(peaks(4,8:2:end), mod(c_fit(3,:), 2*pi()), 'k', '+', 'LineWidth',2); 
% xlabel('photoelectron energy'); 
ylabel('phase (radians)'); 
legend('X', 'b', 'A', 'c'); 
% xlim([10 20]); 
hold off; 

%% CO fit vs. SI
window = 5; 
phase_si = spectral_integration(x1, y1, y3, peaks, window); 
phase_si = reshape(phase_si, [length(IP), length(n)]); 

X_fit = [11, 12, 13, 14, 15, 16, 17, 18, 19; 476.8661, 3.6433, 5.2459, 2.3867, 0.7674, 1.9197, NaN, 1.2200, 3.1016; ... 
    4.1997, 0.7058, 3.7485, 6.2599, 2.8230, 5.5860, NaN, 4.3930, 1.4975]; 
b_fit = [14, 15, 16, 17, 18; 476.8661, 3.6433, 5.2459, 2.3867, 0.7674; ...
    6.0278, 2.9922, 5.8367, 2.3067, 4.9908]; 
A_fit = [14, 15, 16, 17, 18; 476.8661, 3.6433, 5.2459, 2.3867, 0.7674; ...
    0.0566, 2.6320, 5.6394, 2.5555, 4.5583]; 
c_fit = [15, 16, 17, 18, 19; 476.8661, 3.6433, 5.2459, 2.3867, 0.7674; ...
    4.7477, 5.7650, 2.1313, 5.2021, 5.3181]; 

figure; hold on; 
scatter(n(4:2:end), phase_si(1, 4:2:end), 'b', 'o', 'Linewidth', 2); 
scatter(n(4:2:end), X_fit(3, 2:2:end), 'b', '+', 'Linewidth', 2); 
scatter(n(6:2:end), phase_si(2, 6:2:end), 'r', 'o', 'Linewidth', 2); 
scatter(n(6:2:end), b_fit(3, 1:2:end), 'r', '+', 'Linewidth', 2); 
scatter(n(6:2:end), phase_si(3, 6:2:end), 'g', 'o', 'Linewidth', 2); 
scatter(n(6:2:end), A_fit(3, 1:2:end), 'g', '+', 'Linewidth', 2); 
scatter(n(8:2:end), phase_si(4, 8:2:end), 'k', 'o', 'Linewidth', 2);
scatter(n(8:2:end), c_fit(3, 2:2:end), 'k', '+', 'Linewidth', 2);


%% plot time delays
figure; hold on; 
scatter(harmonics(1:2:7), -(b_fit-X_fit)*10^18/(2*omega), 'b', 'o'); 
scatter(harmonics(1:2:7), -(angle(b3Pi_phase(1:2:7))-angle(X_phase(1:2:7)))*10^18/(2*omega), 'r', 'o'); 
xlabel('harmonic number'); 
ylabel('time delay (as)'); 
legend('complex fitting', 'spectral integration'); 
title('b - X time delays'); 
xlim([10 20]); 
hold off; 
