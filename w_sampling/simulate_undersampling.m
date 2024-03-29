% make desired w distribution
% dw_1 = 3E-2; 
dw_1 = 5E-2; 
wrange_1 = 10; 
w_1 = -wrange_1:dw_1:wrange_1; 
slope_2w = 20; 
env_size = 10; 

truncate_1 = 0; 
truncate_2 = 0; 
norm = 5E6; 
numsim = 500; 

bandwidth = 0.03; 
signal_w_1 =  exp(-(w_1-1).^2/(2*bandwidth.^2))  + exp(-(w_1+1).^2/(2*bandwidth.^2)) ...
            + exp(-(w_1-2).^2/(2*bandwidth.^2)).*exp(1j*slope_2w*(w_1-2))  ...
            + exp(-(w_1+2).^2/(2*bandwidth.^2)).*exp(1j*slope_2w*(w_1+2)) ...
            + env_size*exp(-w_1.^2/(2*bandwidth.^2)); 
% FT
signal_t_1 = fftshift(ifft(ifftshift(signal_w_1))); 
signal_t_1 = signal_t_1 * norm/sum(abs(signal_t_1)); 
dt_1 = 1/(numel(w_1))/dw_1; 
trange_1 = 1/dw_1; 
t_1 = ((-trange_1/2):dt_1:(trange_1/2 - dt_1/2)) + dt_1/2; 
fs_1 = 1/dt_1; 
% resample
f0 = 1; 
fs_2 = 2.5*f0; 
dt_2 = 1/fs_2; 
t_2 = ((-trange_1/2):dt_2:(trange_1/2 - dt_1/2)) + dt_1/2;
signal_t_2 = interp1(t_1, signal_t_1, t_2, 'spline'); 

% truncate dense sample to get same number of points between two cases
w_1_full = w_1; 
dw_1_full = dw_1; 
signal_t_1_full = signal_t_1; 
t_1_full = t_1; 
if truncate_1 == 1
    numpts = numel(t_2); 
    midpt_1 = ceil(numel(t_1)/2); 
    if mod(numpts,2)==0
        ind1 = midpt_1 - floor(numpts/2); 
        ind2 = midpt_1 + floor(numpts/2) - 1; 
    else
        ind1 = midpt_1 - floor(numpts/2); 
        ind2 = midpt_1 + floor(numpts/2); 
    end
    signal_t_1_full = signal_t_1; 
    t_1_full = t_1; 
    signal_t_1 = signal_t_1(ind1:ind2); 
    t_1 = t_1(ind1:ind2); 
    w_1_full = w_1; 
    dw_1_full = dw_1; 
    dw_1 =  1/(numel(t_1))/dt_1;
    wrange_1 = 1/dt_1/2;
    w_1 = (-wrange_1):dw_1:(wrange_1-dw_1); 
end
if truncate_1==1 & truncate_2~=0
    midpt_1 = ceil(numel(t_1)/2); 
    ind1 = midpt_1 - floor(truncate_2/2) + 1; 
    ind2 = midpt_1 + floor(truncate_2/2) + 1; 
    signal_t_1 = signal_t_1(ind1:ind2); 
    t_1 = t_1(ind1:ind2); 
    dw_1 =  1/(numel(t_1))/dt_1;
    wrange_1 = 1/dt_1/2;
    w_1 = (-wrange_1):dw_1:(wrange_1-dw_1); 
    
    midpt_2 = ceil(numel(t_2)/2); 
    ind1 = midpt_2 - floor(truncate_2/2) + 1; 
    ind2 = midpt_2 + floor(truncate_2/2) + 1; 
    signal_t_2_full = signal_t_2; 
    signal_t_2 = signal_t_2(ind1:ind2); 
    t_2_full = t_2; 
    t_2 = t_2(ind1:ind2); 
    dw_2 =  1/(numel(t_2))/dt_2;
    wrange_2 = 1/dt_2/2;
    w_2 = (-wrange_2+dw_2):dw_2:wrange_2;    
end


signal_w_1ft = fftshift(fft(ifftshift(signal_t_1))); 
signal_w_2ft = fftshift(fft(ifftshift(signal_t_2))); 
dw_2 =  1/(numel(t_2))/dt_2;
wrange_2 = 1/dt_2/2;
w_2 = (-wrange_2+dw_2):dw_2:wrange_2; 

figure; hold on; 

subplot(3,1,1); hold on; 
yyaxis left; 
plot(w_1_full, abs(signal_w_1)); 
ylabel('amplitude'); 
yyaxis right; 
plot(w_1_full, unwrap(angle(signal_w_1))); 
ylabel('phase'); 
xlabel('w'); 
xlim([-5 5]); 
title('high sample rate'); 
goodplot(); 

subplot(3,1,2); hold on; 
plot(t_1, real(signal_t_1), 'k-', 'DisplayName', 'high sample rate'); 
plot(t_2, real(signal_t_2), 'ro', 'DisplayName', 'low sample rate'); 
xlabel('t'); 
xlim([-20 20]); 
goodplot(); 

subplot(3,1,3); hold on; 
yyaxis left; 
plot(w_2, abs(signal_w_2ft)); 
ylabel('amplitude'); 
yyaxis right; 
plot(w_2, unwrap(angle(signal_w_2ft))); 
ylabel('phase'); 
xlabel('w'); 
title('low sample rate'); 
goodplot()


%%
% scale up to replicate experimental normalization conditions
% signal_t_1 = signal_t_1 * norm/sum(abs(signal_t_1)); 
% signal_t_2 = signal_t_2 * norm/sum(abs(signal_t_2)); 

% make many experiments with poisson noise
signal_t_1_noisy = zeros([numsim numel(t_1)]); 
signal_t_2_noisy = zeros([numsim numel(t_2)]); 
% add noise
for ii=1:numsim
    signal_t_1_noisy(ii,:) = poissrnd(abs(signal_t_1)); 
    signal_t_2_noisy(ii,:) = poissrnd(abs(signal_t_2)); 
end
% FT back
signal_w_1ft = fftshift(fft(ifftshift(signal_t_1))); 
signal_w_1ft_noisy = fftshift(fft(ifftshift(signal_t_1_noisy,2),[],2),2); 
signal_w_2ft = fftshift(fft(ifftshift(signal_t_2))); 
signal_w_2ft_noisy = fftshift(fft(ifftshift(signal_t_2_noisy,2),[],2),2); 
dw_2 =  1/(numel(t_2))/dt_2;
wrange_2 = 1/dt_2/2;
w_2 = (-wrange_2+dw_2):dw_2:wrange_2;  

% analyze phase
phase_w_1 = angle(signal_w_1ft); 
phase_w_1_noisy = angle(signal_w_1ft_noisy); 
for ii=1:numel(w_1)
    if abs(signal_w_1ft(ii)) < max(abs(signal_w_1ft))/100
        phase_w_1(ii) = NaN; 
    end
    for jj=1:numsim
        if abs(signal_w_1ft_noisy(jj,ii)) < max(abs(signal_w_1ft_noisy(jj,:)))/100
            phase_w_1_noisy(jj,ii) = NaN; 
        end
    end
end
phase_w_2 = angle(signal_w_2ft); 
phase_w_2_noisy = angle(signal_w_2ft_noisy); 
for ii=1:numel(w_2)
    if abs(signal_w_2ft(ii)) < max(abs(signal_w_2ft))/100
        phase_w_2(ii) = NaN; 
    end
    for jj=1:numsim
        if abs(signal_w_2ft_noisy(jj,ii)) < max(abs(signal_w_2ft_noisy(jj,:)))/100
            phase_w_2_noisy(jj,ii) = NaN; 
        end
    end
end

%% plot results
figure; hold on; 

subplot(2,1,1); hold on; 
yyaxis left; 
plot(w_1, abs(signal_w_1ft)); 
for ii=1:numsim
    p = plot(w_1, abs(signal_w_1ft_noisy(ii,:)), 'c-'); 
    p.Color(4) = 0.2; 
end
ylabel('amplitude'); 
yyaxis right; 
plot(w_1, phase_w_1); 
for ii=1:numsim
    p = plot(w_1, phase_w_1_noisy(ii,:), 'Color', [1 0.7 0], 'LineStyle', '-', 'Marker', 'none'); 
    p.Color(4) = 0.2; 
end
ylabel('phase'); 
xlabel('w'); 
title('high sample rate'); 
goodplot(); 

subplot(2,1,2); hold on; 
yyaxis left; 
plot(w_2, abs(signal_w_2ft)); 
for ii=1:numsim
    p = plot(w_2, abs(signal_w_2ft_noisy(ii,:)), 'c-'); 
    p.Color(4) = 0.1; 
end
ylabel('amplitude'); 
yyaxis right; 
plot(w_2, phase_w_2); 
for ii=1:numsim
    p = plot(w_2, phase_w_2_noisy(ii,:), 'Color', [1 0.7 0], 'LineStyle', '-', 'Marker', 'none'); 
    p.Color(4) = 0.1; 
end
ylabel('phase'); 
xlabel('w'); 
title('low sample rate'); 
goodplot()

%% extract slope values

slope_list_1w_1 = zeros([1 numsim]); 
phase_list_1w_1 = slope_list_1w_1; 
slope_list_1w_2 = slope_list_1w_1; 
phase_list_1w_2 = slope_list_1w_1; 
slope_list_2w_1 = slope_list_1w_1; 
phase_list_2w_1 = slope_list_1w_1; 
slope_list_2w_2 = slope_list_1w_1; 
phase_list_2w_2 = slope_list_1w_1; 

% find fit regions
window = 2; 
threshold = env_size*5; 
[~,tmp] = findpeaks(abs(signal_w_1ft), 'MinPeakProminence', max(abs(signal_w_1ft)/threshold)); 
center_1w_1 = tmp(4); 
[~,tmp] = findpeaks(abs(signal_w_1ft), 'MinPeakProminence', max(abs(signal_w_1ft)/threshold)); 
center_2w_1 = tmp(5); 
[~,tmp] = findpeaks(abs(signal_w_2ft), 'MinPeakProminence', max(abs(signal_w_2ft)/threshold)); 
center_1w_2 = tmp(5); 
[~,tmp] = findpeaks(abs(signal_w_2ft), 'MinPeakProminence', max(abs(signal_w_2ft)/threshold)); 
center_2w_2 = tmp(4); 


for ii=1:numsim
    
    start = center_1w_1 - window; 
    stop  = center_1w_1 + window; 
    x = w_1(start:stop); 
    y = signal_w_1ft_noisy(ii,start:stop);
    if size(x) == size(y)
        y = y'; 
    end
    IP = -w_1(center_1w_1); IP_label = '1w'; 
    [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, x, abs(y), y, 1, 0); 
    slope_list_1w_1(ii) = paramout(2); 
    phase_list_1w_1(ii) = mod(paramout(1)-pi,2*pi); 
    
    start = center_2w_1 - window; 
    stop  = center_2w_1 + window; 
    x = w_1(start:stop); 
    y = signal_w_1ft_noisy(ii,start:stop);
    if size(x) == size(y)
        y = y'; 
    end
    IP = -w_1(center_2w_1); IP_label = '2w'; 
    [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, x, abs(y), y, 1, 0); 
    slope_list_2w_1(ii) = paramout(2); 
    phase_list_2w_1(ii) = mod(paramout(1)-pi,2*pi); 
    
%     figure; hold on; 
%     yyaxis left; 
%     plot(x, abs(y), 'bo-'); 
%     ylabel('amplitude'); 
%     yyaxis right; 
%     plot(x, unwrap(angle(y)), 'o', 'Color', [1 0.5 0]); 
%     plot(x, polyval([paramout(2) paramout(1)], x), 'r-'); 
%     ylabel('phase'); 
%     title('rate 1 2w'); 
    
    
    start = center_1w_2 - window; 
    stop  = center_1w_2 + window; 
    x = w_2(start:stop); 
    y = signal_w_2ft_noisy(ii,start:stop);
    if size(x) == size(y)
        y = y'; 
    end
    IP = -w_2(center_1w_2); IP_label = '1w'; 
    [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, x, abs(y), y, 1, 0); 
    slope_list_1w_2(ii) = paramout(2); 
    phase_list_1w_2(ii) = mod(paramout(1)-pi,2*pi); 
    
    start = center_2w_2 - window; 
    stop  = center_2w_2 + window; 
    x = w_2(start:stop); 
    y = signal_w_2ft_noisy(ii,start:stop);
    if size(x) == size(y)
        y = y'; 
    end
    IP = -w_2(center_2w_2); IP_label = '2w'; 
    [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, x, abs(y), y, 1, 0); 
    slope_list_2w_2(ii) = paramout(2); 
    phase_list_2w_2(ii) = mod(paramout(1)-pi,2*pi); 
    
%     figure; hold on; 
%     yyaxis left; 
%     plot(x, abs(y), 'bo-'); 
%     ylabel('amplitude'); 
%     yyaxis right; 
%     plot(x, unwrap(angle(y)), 'o', 'Color', [1 0.5 0]); 
%     plot(x, polyval([paramout(2) paramout(1)], x), 'r-'); 
%     ylabel('phase'); 
%     title('rate 2 2w'); 
    
end

%%
% histogram and evaluate
found_params_array = [slope_list_1w_1; ...
                      phase_list_1w_1; ...
                      slope_list_2w_1; ...
                      phase_list_2w_1; ...
                      slope_list_1w_2; ...
                      phase_list_1w_2; ...
                      slope_list_2w_2; ...
                      phase_list_2w_2]; 
found_params_histfit = zeros([2 8]); 
for ii=1:8
    [N, edges] = histcounts(found_params_array(ii,:)); 
    f = fit(edges(2:end)', N', 'gauss1'); 
    found_params_histfit(1,ii) = f.b1; 
    found_params_histfit(2,ii) = f.c1; 
end

figure; hold on; 

subplot(4, 2, 1); 
histogram(slope_list_1w_1); 
xlabel('slope'); title('1w, rate 1'); 
text(3, 50, ...
    strcat(num2str(found_params_histfit(1,1)), ' +/- ', num2str(found_params_histfit(2,1)))); 
subplot(4, 2, 3); 
histogram(slope_list_2w_1); 
xlabel('slope'); title('2w, rate 1'); 
text(-7, 50, ...
    strcat(num2str(found_params_histfit(1,3)), ' +/- ', num2str(found_params_histfit(2,3)))); 
subplot(4, 2, 5); 
histogram(phase_list_1w_1); 
xlabel('phase'); title('1w, rate 1'); 
text(3, 50, ...
    strcat(num2str(found_params_histfit(1,2)), ' +/- ', num2str(found_params_histfit(2,2)))); 
subplot(4, 2, 7); 
histogram(phase_list_2w_1); 
xlabel('phase'); title('2w, rate 1'); 
text(3, 50, ...
    strcat(num2str(found_params_histfit(1,4)), ' +/- ', num2str(found_params_histfit(2,4)))); 

subplot(4, 2, 2); 
histogram(slope_list_1w_2); 
xlabel('slope'); title('1w, rate 2'); 
text(3, 50, ...
    strcat(num2str(found_params_histfit(1,5)), ' +/- ', num2str(found_params_histfit(2,5)))); 
subplot(4, 2, 4); 
histogram(slope_list_2w_2); 
xlabel('slope'); title('2w, rate 2'); 
text(-18, 50, ...
    strcat(num2str(found_params_histfit(1,7)), ' +/- ', num2str(found_params_histfit(2,7)))); 
subplot(4, 2, 6); 
histogram(phase_list_1w_2); 
xlabel('phase'); title('1w, rate 2'); 
text(3, 50, ...
    strcat(num2str(found_params_histfit(1,6)), ' +/- ', num2str(found_params_histfit(2,6)))); 
subplot(4, 2, 8); 
histogram(phase_list_2w_2); 
xlabel('phase'); title('2w, rate 2');  
text(3, 50, ...
    strcat(num2str(found_params_histfit(1,8)), ' +/- ', num2str(found_params_histfit(2,8)))); 













