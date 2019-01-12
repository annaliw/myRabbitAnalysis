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

%% compare to fit extraction
% X_fit = [0.6912 1.0608 1.8472 3.1419]; 
% X_fit = [0.7 1.0608 1.8472 3.1419]; 
% b_fit = [0.1053 1.1498 1.7822 2.8808]; 
X_fit = [10, 11, 12, 13, 14, ; 175.9033, 43.9565, 44.1536, 11.8212, 373.9577;...
    4.5433, 4.8752, 0.3235, 3.8694, ]; 
b_fit = [12, 13, 14, 15, 16, 17, 18, 19, ; 713.9706, 91.7449, 53.8667, 43.9565, 44.1536, 11.8212, 373.9577, ; ...
    6.2832, 4.0629, 1.2348, 2.3675, 1.6701, 5.7794, ]; % harmonic, resnorm, phase in columns
A_fit = [14, 16, 17, 18, 19; 91.7449, 29.2274, 44.1536, 11.8212, 373.9577; 2.0391, 1.8822, 2.0423, 4.9629, 3.3323]; % 16 ?
c_fit = [16, 18, ; 91.7449, 29.2274, ; 1.0662, 2.5632, ]; % 18 ?
    
omega = 3*10^8/(810*10^(-9)); 

figure; hold on; 
% scatter(harmonics(1:2:7), angle(b3Pi_phase(1:2:7)), 'b', 'o'); 
% scatter(harmonics(1:2:7), angle(X_phase(1:2:7)), 'r', 'o');
scatter(X_fit(1,1:2:end), mod(X_fit(3,1:2:end), 2*pi()), 'r', '+', 'LineWidth',2); 
scatter(b_fit(1,1:2:end), mod(b_fit(3,1:2:end), 2*pi()), 'b', '+', 'LineWidth',2); 
scatter(A_fit(1,1:2:end), mod(A_fit(3,1:2:end), 2*pi()), 'g', '+', 'LineWidth',2); 
scatter(c_fit(1,1:2:end), mod(c_fit(3,1:2:end), 2*pi()), 'k', '+', 'LineWidth',2); 
xlabel('harmonic number'); 
ylabel('phase (radians)'); 
% legend('b SI', 'X SI', 'X fit', 'b fit'); 
legend('X', 'b', 'A', 'c'); 
xlim([10 20]); 
hold off; 

%% 
figure; hold on; 
scatter(harmonics(1:2:7), -(b_fit-X_fit)*10^18/(2*omega), 'b', 'o'); 
scatter(harmonics(1:2:7), -(angle(b3Pi_phase(1:2:7))-angle(X_phase(1:2:7)))*10^18/(2*omega), 'r', 'o'); 
xlabel('harmonic number'); 
ylabel('time delay (as)'); 
legend('complex fitting', 'spectral integration'); 
title('b - X time delays'); 
xlim([10 20]); 
hold off; 
