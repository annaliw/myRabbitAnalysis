
%% from Andrei's original time series to undersampling
% length units in microns, time in s
c = 2.998E14; 
l = 0.810; 
f = c/l; 
w = 2*pi*f; 

t1 = (-5:0.018:5)*2 ./ c; 
figure; hold on; 
subplot(2,1,1)
plot(t1, cos(w*t1) + cos(w*t1*2)); 
xlabel('time'); 

dt_1 = 2*0.018/c; 
fs_1 = 1/dt_1; 
subplot(2,1,2)
df_1 = 1/(t1(end)-t1(1)); 
f1 = (-(numel(t1)-1)/2:1:(numel(t1)-1)/2)*df_1*l/c; 
plot(f1, abs(fftshift(fft(cos(w*t1) + cos(w*t1*2))))); 
xlabel('frequency'); 

fs_2 = 2.45*f; 
dt_2 = 1/(fs_2); 
t2 = (-10/c:dt_2:10/c); 
figure; hold on; 
subplot(2,1,1)
plot(t2, cos(w*t2) + cos(w*t2*2)); 
xlabel('time'); 
subplot(2,1,2)
df_2 = 1/(t2(end)-t2(1)); 
f2 = (-(numel(t2)-1)/2:1:(numel(t2)-1)/2)*df_2*l/c; 
plot(f2, abs(fftshift(fft(cos(w*t2) + cos(w*t2*2))))); 
xlabel('frequency'); 





