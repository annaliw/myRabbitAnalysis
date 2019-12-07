% 1000 points
% 2 period range
periods = 10; 
dx = periods*2*pi/1000; 
x = dx:dx:(2*periods*pi); 
y = cos(x); 

y_ft = fftshift(fft(y)); 
dx_ft = 1/x(end); 
x_ft = ((1:1000)-500)*dx_ft - dx_ft; 

figure; plot(x, y, 'o'); 
figure; plot(x_ft*2*pi, abs(y_ft)); 

%% use specific sample rate
% periods = 10; 
% xrate = 35; 
periods = 50; 
xrate = 0.4; 
dx = 1/xrate; 
x = dx:dx:(2*pi*periods); 
y = cos(x) + cos(2*x); 
numpts = numel(x); 

y_ft = fftshift(fft(y)); 
dx_ft = 1/(x(end) - x(1)); 
% x_ft = ((1:numpts)-numpts/2)*dx_ft - dx_ft; 
x_ft = ((1:numpts)-1 - numpts/2)*xrate/numpts; 

figure; plot(x, y, '.'); 
figure; plot(x_ft*2*pi, abs(y_ft)); 

%%
c = 2.998e8; 
l = 810e-9; 
f = c/l; 
w = 2*pi*f; 

stage = 4e-6; 
distance = 2*stage; 
time = distance/c; 
oscillations = time*2*f

dt = 2*0.018e-6/c
rate = oscillations / dt / 2/w

%% check phase retention

% create desired spectrum w/ phase
dw = 0.05; 
wlim = 100; 
warray = -wlim:dw:wlim; 
% bandwidth 0.05w
% flat phase for w, sloped phase for 2w
bw = 0.05; 
% spectrum_1 = exp(-(warray - 1).^2/(2*bw^2)) .* exp(-1j*pi) ...
%            + exp(-(warray + 1).^2/(2*bw^2)) .* exp(1j*pi) ...
%            + exp(-(warray-2).^2/(2*bw^2)) .* exp(-1j*(warray-2)) ...
%            + exp(-(warray+2).^2/(2*bw^2)) .* exp(1j*(warray+2)); 

spectrum_1 = exp(-(warray - 1).^2/(2*bw^2)) ...
           + exp(-(warray + 1).^2/(2*bw^2));
%                 + exp(-(warray-2).^2/(2*bw^2)) ...
%            + exp(-(warray+2).^2/(2*bw^2)); 

figure; hold on; 
yyaxis right
plot(warray, abs(spectrum_1)); 
ylabel('amplitude'); 
yyaxis left
plot(warray, unwrap(angle(spectrum_1))); 
ylabel('phase'); 
xlabel('w');
title('original spectrum'); 
goodplot(); 

tarray = (-1/dw):(1/wlim):(1/dw); 
spectrum_ft_1 = ifft(spectrum_1); 

figure; hold on; 
plot(tarray, ifftshift(real(spectrum_ft_1))); 

%%
oscillations = 1/dw; 
sample_rate_1 = numel(tarray)/oscillations; 

sample_rate_2 = floor(sample_rate_1)/2; 
numpts_2 = sample_rate_2 * oscillations; 
skip = floor(numel(tarray)/numpts_2); 
tarray_2 = tarray(1:skip:end); 
spectrum_ft_2 = spectrum_ft_1(1:skip:end); 

figure; hold on; 
plot(tarray_2(1:end), fftshift(real(spectrum_ft_2))); 

dt_2 = tarray_2(2)-tarray_2(1); 
warray_2 = (-1/dt_2):(2/(tarray_2(end)-tarray_2(1))):(1/dt_2); 
spectrum_2 = fft(real(spectrum_ft_2)); 

figure; hold on; 
yyaxis right; 
plot(warray_2, fftshift(abs(spectrum_2))); 
ylabel('amplitude'); 
yyaxis left; 
plot(warray_2, fftshift(unwrap(angle(spectrum_2)))); 
ylabel('phase'); 
xlabel('w'); 



