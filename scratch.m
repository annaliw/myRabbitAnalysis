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
X_fit = [0.7 1.0608 1.8472 3.1419]; 
b_fit = [0.1053 1.1498 1.7822 2.8808]; 
omega = 3*10^8/(810*10^(-9)); 

figure; hold on; 
scatter(harmonics(1:2:7), angle(b3Pi_phase(1:2:7)), 'b', 'o'); 
scatter(harmonics(1:2:7), angle(X_phase(1:2:7)), 'r', 'o');
scatter(harmonics(1:2:7), X_fit, 'r', '+'); 
scatter(harmonics(1:2:7), b_fit, 'b', '+'); 
xlabel('harmonic number'); 
ylabel('phase (radians)'); 
legend('b SI', 'X SI', 'X fit', 'b fit'); 
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
