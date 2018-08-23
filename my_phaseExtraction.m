%% get energy vs delay time
X_Eng = fliplr([17.82, 14.79, 11.68, 8.734]); 
b_Eng = [1.825, 4.852, 7.958, 10.99]; 
A_Eng = [3.067, 6.172, 9.2]; 
B_Eng = [2.756, 5.939]; 

window = 2; 
X_phi = X_Eng; 
X_phi_list = zeros(5, length(X_Eng)); 
for j = 1:1:5 
    subtract = j-3; 
    for i=1:1:length(X_Eng)
        ind = find(abs(E-b_Eng(i))<0.005)-subtract;  
        X_phi(i) = sum(twoOmega_phi(ind-window:ind+window).*twoOmega_abs(ind-window:ind+window))/sum(twoOmega_abs(ind-window:ind+window)); 
    end
    X_phi_list(j, :) = X_phi; 
end
% b_phi = b_Eng; 
% for i=1:1:length(b_Eng)
%     ind = find(abs(E-b_Eng(i))<0.005);  
%     b_phi(i) = sum(twoOmega_phi(ind-window:ind+window).*twoOmega_abs(ind-window:ind+window))/sum(twoOmega_abs(ind-window:ind+window)); 
% end

sideband = [12, 14, 16, 18]; 
% X_tau = X_phi/(2*2.9979e8/810e-09); 
% b_tau = b_phi/(2*2.9979e8/810e-09); 
figure; hold on; 
% scatter(sideband, X_tau);
% scatter(sideband, b_tau);
for j=1:1:5
    scatter(sideband, X_phi_list(j,:)); 
end
legend; 
hold off; 

%% try fitting to complex FFT results

% engWindow = [12.15 13.78]; 
% engWindow = [10.62 12.15]; 
% engWindow = [11.28 12.78]; 
engWindow = [9.014 9.726]; 
% % engWindow = [7.616 8.125]; 
% % engWindow = [8.637 8.988]; 
indWindow = fliplr([find(abs(E-engWindow(1))<0.005) find(abs(E-engWindow(2))<0.005)]); 
xin = E(indWindow(1)+3:indWindow(2)); 
yin = twoOmega_signal(indWindow(1)+3:indWindow(2)).'; 
% guess = [3e04 3e04 0.4 0.6 -pi/2 0 -0.6 -1.3 0]; 
% guess = [8e03 8e03 0.4 0.6 0 0 3 1.1 0]; 
% guess = [8e03 8e03 0.6 0.4 0 -0.6 1.22 0 0];
guess = [4e03 8e03 0.07 0.4 0 0 -1.6 -1.2 0]; 
% % guess = [0.4e05 2.5e05 0.005 0.2 1.2 1.9 0]; 
% % guess = [3e04 3e04 0.25 0.25 0.59 1.01 0]; 
% peaks = [12.44 13.1]; 
% peaks = [10.9 11.6]; 
% peaks = [11.72 12.48]; 
peaks = [9.2 9.5]; 
% % peaks = [7.698 7.95]; 
% % peaks = [8.662 8.835]; 

[paramout, fval, x_out, y_out] = my2wfit(xin, yin, peaks, guess); 
A1 = paramout(1); 
A2 = paramout(2); 
s1 = paramout(3);
s2 = paramout(4); 
% a1 = paramout(5); 
% a2 = paramout(6); 
a1 = 0; 
a2 = 0; 
b1 = paramout(7); 
b2 = paramout(8); 
o  = paramout(9); 

x1 = peaks(1); 
x2 = peaks(2); 

y1 = abs(A1) * exp(1i*(a1*(x_out-x1)+b1)) .* exp(-(x_out-x1).^2/(2*s1)); 
y2 = abs(A2) * exp(1i*(a2*(x_out-x2)+b2)) .* exp(-(x_out-x2).^2/(2*s2));

% y1 = abs(A1) * exp(1i*(b1)) .* exp(-(x_out-x1).^2/(2*s1)); 
% y2 = abs(A2) * exp(1i*(-pi/2*(x_out-x2)+b2)) .* exp(-(x_out-x2).^2/(2*s2));
% y_out = y1+y2; 

figure; 
hold on; 
yyaxis right; 
scatter(xin, abs(yin), 'o'); plot(x_out, abs(y_out), '-'); 
plot(x_out, abs(y1), '--'); plot(x_out, abs(y2), '-.'); 
ylabel('amplitude'); 
yyaxis left; 
scatter(xin, angle(yin), '+'); plot(x_out, angle(y_out), '-'); 
plot(x_out, angle(y1), '--'); plot(x_out, angle(y2), '-.'); 
ylabel('phase'); 
xlabel('energy'); 
hold off; 

%% plotting harmonic vs. phase 

X18_phase = pi; % read off plot
X16_phase = 1.879; % read off plot
X14_phase = mod((param_X14b19(8) + param_X14b18(7))/2, 2*pi); % fit extracted
X12_phase = 1.105; % read off plot, not good data
X10_phase = 1.942; % read off plot, not good data

b18_phase = mod(param_X14b18(8), 2*pi);
b16_phase = 1.797; %read off of plot
b14_phase = 1.19; 
b12_phase = 0.56; 

X_phase_list = [X12_phase X14_phase X16_phase X18_phase]; 
b_phase_list = [b12_phase b14_phase b16_phase b18_phase]; 
harmonics = [12 14 16 18]; 

figure; hold on; 
scatter(harmonics, b_phase_list-X_phase_list); 




