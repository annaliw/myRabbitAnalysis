% make desired w distribution
dw_1 = 1E-2; 
wrange_1 = 5; 
w_1 = -wrange_1:dw_1:wrange_1; 

bandwidth = 3*dw_1; 
signal_w_1 =  exp(-(w_1-1).^2/(2*bandwidth.^2))  + exp(-(w_1+1).^2/(2*bandwidth.^2)) ...
            + exp(-(w_1-2).^2/(2*bandwidth.^2)).*exp(-1j*(w_1-2))  + exp(-(w_1+2).^2/(2*bandwidth.^2)).*exp(1j*(w_1+2)) ...
            + 50*exp(-w_1.^2/(2*bandwidth.^2)); 
% signal_w_1 =  exp(-(w_1-1).^2/(2*bandwidth.^2))  + exp(-(w_1+1).^2/(2*bandwidth.^2)) ...
%             + exp(-(w_1-2).^2/(2*bandwidth.^2))  + exp(-(w_1+2).^2/(2*bandwidth.^2)); 

        
figure; hold on; 

subplot(3,1,1); 
yyaxis left
plot(w_1, abs(signal_w_1)); 
ylabel('amplitude'); 
yyaxis right
plot(w_1, angle(signal_w_1)); 
ylabel('phase'); 
xlabel('w'); 

% FT
signal_t_1 = fftshift(ifft(ifftshift(signal_w_1))); 
dt_1 = 1/(numel(w_1))/dw_1; 
trange_1 = 1/dw_1; 
t_1 = ((-trange_1/2):dt_1:(trange_1/2 - dt_1/2)) + dt_1/2; 
fs_1 = 1/dt_1; 

% resample
f0 = 1; 
fs_2 = 2.4*f0; 
dt_2 = 1/fs_2; 
% t_2 = (-trange_1/2:dt_2:(trange_1/2 - dt_1/2 + dt_2)); 
t_2 = ((-trange_1/2):dt_2:(trange_1/2 - dt_1/2)) + dt_1/2; %+ (t_1(end)-t_2(end))/2*dt_1; 
signal_t_2 = interp1(t_1, signal_t_1, t_2, 'spline'); 
% signal_t_2 = signal_t_2(2:end); 
% t_2 = t_2(2:end); 

% % subtract t0 offset
% t0_1 = t_1( find( abs(t_1-0) < dt_1, 1) ); 
% t0_2 = t_2( find( abs(t_2-0) < dt_2, 1) ); 
% t_2 = 

subplot(3,1,2); hold on; 
plot(t_1, real(signal_t_1)); 
plot(t_2, real(signal_t_2), 'o'); 
xlabel('t'); 

% FT back
signal_w_2 = fftshift(fft(ifftshift(signal_t_2))); 
dw_2 =  1/(numel(t_2))/dt_2;
wrange_2 = 1/dt_2/2;
w_2 = (-wrange_2+dw_2):dw_2:wrange_2; 

subplot(3,1,3); hold on; 
yyaxis left
plot(w_2, abs(signal_w_2)); 
ylabel('amplitude'); 
yyaxis right
plot(w_2, angle(signal_w_2)); 
ylabel('phase'); 

%% check effect of sample rate on phase slope

w1_slope_list = zeros([1 10]); 
w1_phase_list = w1_slope_list; 
w2_slope_list = w1_slope_list; 
w2_phase_list = w1_slope_list; 

numpts = 20; 
center = 2.5; 
step = 0.02; 
start_rate = center - step*numpts; 
for ii=1:numpts
    
    % resample
    f0 = 1; 
    fs_2 = (start_rate+step*ii)*f0; 
    dt_2 = 1/fs_2; 
    t_2 = ((-trange_1/2):dt_2:(trange_1/2 - dt_1/2)) + dt_1/2; %+ (t_1(end)-t_2(end))/2*dt_1; 
    signal_t_2 = interp1(t_1, signal_t_1, t_2, 'spline'); 
    % FT back
    signal_w_2 = fftshift(fft(ifftshift(signal_t_2))); 
    dw_2 =  1/(numel(t_2))/dt_2;
    wrange_2 = 1/dt_2/2;
    w_2 = (-wrange_2+dw_2):dw_2:wrange_2; 
    
    % fit to find slopes
    start = find(abs(w_2 - 0.93) < dw_2); 
    stop  = find(abs(w_2 - 1.07) < dw_2); 
    p = polyfit(w_2(start:stop), angle(signal_w_2(start:stop)), 1); 
    w1_slope_list(ii) = p(1); 
    w1_phase_list(ii) = p(2); 
    start = find(abs(w_2 - 1.93 + start_rate+step*ii) < dw_2); 
    stop  = find(abs(w_2 - 2.07 + start_rate+step*ii) < dw_2); 
    p = polyfit(w_2(start:stop), angle(signal_w_2(start:stop)), 1); 
    w2_slope_list(ii) = p(1); 
    w2_phase_list(ii) = p(2); 
    
end

figure; hold on; 
plot(start_rate+step*(1:numpts), w1_slope_list, 'o-', 'DisplayName', 'w1 slope'); 
plot(start_rate+step*(1:numpts), w2_slope_list, 'o-', 'DisplayName', 'w2 slope'); 
xlabel('sample rate'); 
ylabel('phase slope'); 
legend; 

figure; hold on; 
plot(start_rate+step*(1:numpts), w1_phase_list, 'o-', 'DisplayName', 'w1 phase'); 
plot(start_rate+step*(1:numpts), w2_phase_list, 'o-', 'DisplayName', 'w2 phase'); 
xlabel('sample rate'); 
ylabel('phase offset'); 
legend; 

%% add noise before resample 

multiplier = 1000/max(abs(signal_t_1)); 
numsim = 1000; 
% figure; hold on; 
% plot(t_1, abs(signal_t_1*multiplier)); 
% figure; hold on; 
% plot(t_1, abs(signal_t_1*multiplier) + shot_noise); 

twoOmega_values = zeros([numsim 5]); 
twoOmega_traces = zeros([numsim numel(t_2)]); 

for ii=1:numsim
    % resample
    f0 = 1; 
    fs_3 = 2.4*f0; 
    dt_3 = 1/fs_2; 
    t_3 = ((-trange_1/2):dt_2:(trange_1/2 - dt_1/2)) + dt_1/2; 
    signal_t_3_full = interp1(t_1, signal_t_1*multiplier, t_3, 'spline'); 
    signal_t_3 = poissrnd(abs(signal_t_3_full)); 
    % FT back
    signal_w_3 = fftshift(fft(ifftshift(signal_t_3))); 
    dw_3 =  1/(numel(t_3))/dt_3;
    wrange_3 = 1/dt_3/2;
    w_3 = (-wrange_3+dw_3):dw_3:wrange_3;  

    % figure; hold on; 
    % plot(w_3, abs(signal_w_3)); 
    [pks, loc] = findpeaks(abs(signal_w_3)); 
    [pks, tmp] = sort(pks); loc = loc(tmp); 
    pks = fliplr(pks); loc = fliplr(loc); 
    % plot(w_3(loc(1:5)), abs(signal_w_3(loc(1:5))), 'ro'); 

    tmp = sort(loc(1:5)); 
    twoOmega_center = w_3(tmp(4)); 
    twoOmega_amp = abs(signal_w_3(tmp(4))); 
    DC_amp = abs(signal_w_3(tmp(3))); 
    p = polyfit(w_3((tmp(4)-3):(tmp(4)+3))-twoOmega_center, angle(signal_w_3((tmp(4)-3):(tmp(4)+3))), 1); 
    twoOmega_phase = p(2); twoOmega_slope = p(1); 
    
    twoOmega_values(ii,:) = [twoOmega_center, twoOmega_amp, DC_amp, twoOmega_phase, twoOmega_slope]; 
    twoOmega_traces(ii,:) = signal_t_3; 
end



%%
figure; hold on; 
for ii=1:100:numsim
    p = plot(t_3, signal_t_3, 'Color', 'c', 'LineWidth', 1); 
    p.Color(4) = 0.1; 
end

%% make figure to compare

multiplier = 1000/max(abs(signal_t_1)); 

f0 = 1; 
fs_3 = 2.4*f0; 
dt_3 = 1/fs_2; 
t_3 = ((-trange_1/2):dt_2:(trange_1/2 - dt_1/2)) + dt_1/2; 
signal_t_3_full = interp1(t_1, signal_t_1*multiplier, t_3, 'spline'); 
signal_t_3 = poissrnd(abs(signal_t_3_full)); 
% FT back
signal_w_3 = fftshift(fft(ifftshift(signal_t_3))); 
dw_3 =  1/(numel(t_3))/dt_3;
wrange_3 = 1/dt_3/2;
w_3 = (-wrange_3+dw_3):dw_3:wrange_3;  
    
figure; hold on; 

subplot(3,1,1); 
yyaxis left
plot(w_1, abs(signal_w_1*multiplier)); 
ylabel('amplitude'); 
yyaxis right
plot(w_1, angle(signal_w_1)); 
ylabel('phase'); 
xlabel('w'); 
goodplot()

subplot(3,1,2); hold on; 
plot(t_1, real(signal_t_1*multiplier)); 
plot(t_2, real(signal_t_2*multiplier), 'o'); 
plot(t_3, real(signal_t_3), '*'); 
xlabel('t');
goodplot()

subplot(3,1,3); hold on; 
yyaxis left
plot(w_2, abs(signal_w_2*multiplier)); 
plot(w_3, abs(signal_w_3)); 
ylabel('amplitude'); 
yyaxis right
plot(w_2, angle(signal_w_2*multiplier)); 
plot(w_3, angle(signal_w_3)); 
ylabel('phase'); 
xlabel('w'); 
goodplot()


















