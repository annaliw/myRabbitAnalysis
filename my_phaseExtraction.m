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

engWindow = [12.15 13.78]; 
% engWindow = [10.62 12.15]; 
indWindow = fliplr([find(abs(E-engWindow(1))<0.005) find(abs(E-engWindow(2))<0.005)]); 
xin = E(indWindow(1):indWindow(2)); 
yin = twoOmega_signal(indWindow(1):indWindow(2)).'; 
% guess = [8e03 8e03 10.9 11.6 0.4 0.6 3 1.1 0]; 
guess = [4e03 4e03 12.44 13.1 0.4 0.6 -0.7 -1.3 0]; 

[paramout, fval, x_out, y_out] = my2wfit(xin, yin, guess); 
A1 = paramout(1); 
A2 = paramout(2); 
x1 = paramout(3); 
x2 = paramout(4); 
s1 = paramout(5);
s2 = paramout(6); 
b1 = paramout(7); 
b2 = paramout(8); 
o  = paramout(9); 

% y1_r = abs(A1) * cos(b1) .* exp(-(x_out-x1).^2/(2*s1)) + o; 
% y1_i = abs(A1) * sin(b1) .* exp(-(x_out-x1).^2/(2*s1)); 
% y1 = y1_r + 1i*y1_i; 
% 
% y2_r = abs(A2) * cos(b2) .* exp(-(x_out-x2).^2/(2*s2)) + o; 
% y2_i = abs(A2) * sin(b2) .* exp(-(x_out-x2).^2/(2*s2)); 
% y2 = y2_r + 1i*y2_i; 

y1 = abs(A1) * exp(1i*b1) .* exp(-(x_out-x1).^2/(2*s1)); 
y2 = abs(A2) * exp(1i*b2) .* exp(-(x_out-x2).^2/(2*s2)); 

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