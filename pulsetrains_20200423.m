% make a pulse train

centers = -15:1:15; 
% centers = 3*(0.2:0.4:4).^2; 
% centers = cat(2,-fliplr(centers(2:end)), centers); 
% centers = (0.2.^(1:1:30)); 
phases = zeros(size(centers)); 
% phases = -centers.^2/100; 
% phases = -(1:numel(centers))/20
A = 1; 
b = 0.05; 

xtrain = -20:0.1:20; 
ytrain = zeros(size(xtrain)); 
for ii=1:2:numel(centers)
    ytrain = ytrain + A*exp(-(xtrain-centers(ii)).^2/(2*b^2)) .* exp(-1j*phases(ii)); 
end
for ii=2:2:numel(centers)
    ytrain = ytrain + A*exp(-(xtrain-centers(ii)).^2/(2*b^2)) .* exp(-1j*phases(ii)) .* exp(-1j*pi); 
end
% add envelope
ytrain = ytrain .* exp(-(xtrain).^2/(2*((centers(end)-centers(1))/10).^2)); 
% % add "resonance" event
% ytrain = ytrain .* exp(1j .* exp(-(xtrain-5).^2/(2*0.05^2))); 
% % add phase that I fit...
% % p = [4.7004e-06 -7.6579e-08 -5.4032e-06 5.8933e-07 4.2396e-04 -1.7636e-06 -0.0180]; 
% ptrain = polyval(p, xtrain); 
% ytrain = ytrain .* exp(1j * ptrain); 
% % add modified 1/x phase
% x = xtrain + 23; 
% distortion = 1./x; 
% ytrain = ytrain .* exp(1j * distortion); 

distortion = 0.2*atan(4*(xtrain-1.8)); 
ytrain = ytrain .* exp(1j * distortion); 


figure; hold on; 
subplot(1,2,1); 
yyaxis left; 
plot(xtrain, abs(ytrain)); 
xlim([min(centers) max(centers)]); 
ylim([0 1.5*max(abs(ytrain))]); 
ylabel('amplitude'); 
yyaxis right; 
plot(xtrain, unwrap(angle(ytrain))); 
xlabel('time'); 
title('attosecond pulse train'); 
goodplot()

% fourier transform
dw = 1/(numel(xtrain))/(xtrain(2)-xtrain(1)); 
wrange = 1/(xtrain(2)-xtrain(1)); 
xcomb = ((-wrange/2):dw:(wrange/2 - dw/2)) + dw/2; 
ycomb = ifftshift(ifft(fftshift(ytrain))); 

subplot(1,2,2); 
yyaxis left; 
plot(xcomb, abs(ycomb)); 
ylim([0 1.5*max(abs(ycomb))]); 
yyaxis right; 
plot(xcomb, unwrap(angle(ycomb))); 
ylabel('phase'); 
xlabel('frequency'); 
title('frequency comb'); 
goodplot()

%% try going the other way? 

centers = -7:2:7; 
% centers = (0:1:4).^2/2; 
% centers = cat(2,-fliplr(centers), centers); 
phases = centers/3; 
% phases = -centers.^2/100; 
A = 1; 
b = 0.1; 

xcomb = -10:0.01:10; 
ycomb = zeros(size(xcomb)); 
for ii=1:numel(centers)
    ycomb = ycomb + A*exp(-(xcomb-centers(ii)).^2/(2*b^2)); 
end
% add envelope
ycomb = ycomb .* exp(-xcomb.^2/(2*((centers(end)-centers(1))/5).^2)); 
% % add "resonance" event
% ycomb = ycomb .* exp(1j .* exp(-(xcomb-5).^2/(2*0.05^2))); 
% add sloped phase to entire comb
ycomb = ycomb .* exp(-1j*xcomb/60); 
% add polynomail distortion to one comb tooth
[~,x0] = min(abs(xcomb-centers(3))); 
x = xcomb - xcomb(x0); 
% distortion = 10*((1000*(4*x).^3/factorial(3)- 1000) .* 0.5.*exp(-(xcomb-xcomb(x0)).^2/(2*0.3^2)))/1000; 
distortion = (-(4*x).^2 - 4) .* 0.5.*exp(-x.^2/(2*0.3^2)); 
% ycomb = ycomb .* exp(-1j*distortion/20); 



figure; hold on; 

subplot(1,2,1); 
yyaxis left; 
plot(xcomb, abs(ycomb .* exp(-1j*distortion/20))); 
% xlim([min(centers) max(centers)]); 
ylim([0 1.5*max(abs(ycomb .* exp(-1j*distortion/20)))]); 
ylabel('amplitude'); 
yyaxis right; 
plot(xcomb, unwrap(angle(ycomb .* exp(-1j*distortion/20)))); 
% ylim([-0.06 -0.04]); 
xlabel('frequency'); 
title('frequency comb'); 
goodplot()

% subplot(1,3,2); 
% yyaxis left; 
% plot(xcomb, abs(ycomb)); 
% % xlim([min(centers) max(centers)]); 
% ylim([0 1.5*max(abs(ycomb))]); 
% ylabel('amplitude'); 
% yyaxis right; 
% plot(xcomb, angle(ycomb)); 
% % ylim([-0.06 -0.04]); 
% xlim([xcomb(x0)-0.3 xcomb(x0)+0.3]); 
% xlabel('frequency'); 
% title('(zoom in on comb)'); 
% goodplot()

% fourier transform
dt = 1/(numel(xcomb))/(xcomb(2)-xcomb(1)); 
trange = 1/(xcomb(2)-xcomb(1)); 
xtrain = ((-trange/2):dt:(trange/2 - dt/2)) + dt/2; 
ytrain_1 = fftshift(fft(ifftshift(ycomb))); 
ytrain_2 = fftshift(fft(ifftshift(ycomb .* exp(-1j*distortion/20)))); 

subplot(1,2,2); hold on; 
% yyaxis left; 
plot(xtrain, abs(ytrain_1));
plot(xtrain, abs(ytrain_2)); 
ylim([0 1.5*max(real(ytrain_1))]); 
% yyaxis right; 
% [pks, locs] = findpeaks(abs(ytrain)); 
% plot(xtrain(locs(1:2:end)), angle(ytrain(locs(1:2:end))), 'o-'); 
% plot(xtrain(locs(2:2:end)), angle(ytrain(locs(2:2:end)))+pi, 'x-'); 
xlim([-3.5 3.5]); 
% ylim([-900 -800]); 
ylabel('phase'); 
xlabel('time'); 
title('pulse train'); 
goodplot()

% this made a 6th degree phase polynomial. plug that into first cell. 

%% try to get phase polynomial from data? 

[~,i1] = min(abs(E-2.1));  
[~,i2] = min(abs(E-2.22)); 
xval = E(i1:i2); 
yval = twoOmega_signal(i1:i2)'; 

[~,i1] = min(abs(E-5.2));  
[~,i2] = min(abs(E-5.35)); 
xval = cat(2, xval, E(i1:i2)); 
yval = cat(2, yval, twoOmega_signal(i1:i2)'); 

[~,i1] = min(abs(E-8.25));  
[~,i2] = min(abs(E-8.45)); 
xval = cat(2, xval, E(i1:i2)); 
yval = cat(2, yval, twoOmega_signal(i1:i2)'); 

figure; hold on; 
plot(xval, angle(yval), 'o'); 

p = polyfit(xval, angle(yval), 5); 
x = min(xval):0.1:max(xval); 
plot(x, polyval(p,x), 'k-'); 




