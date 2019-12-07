%% Argon photoionization cross section
pics_x = [15.83 15.89 16.10 16.31 16.53 16.75 16.98 17.22 17.46 17.71 17.97 ... 
    18.23 18.50 18.78 19.07 19.37 19.68 20.00 20.32 20.66 21.01 21.38 21.75 ... 
    22.14 22.54 22.96 23.39 23.84 24.31 24.80 25.30 25.83 26.38 26.95 27.55 ...
    28.18 28.83 29.52 30.24 30.99 31.79 32.63 33.51 34.44 35.42 36.46 37.57]; 
pics_y = [29.2 29.5 30.3 31.1 31.8 32.5 33.1 33.7 34.2 34.7 35.1 35.5 35.8 ...
    36.1 36.3 36.5 36.6 36.7 36.8 36.7 36.7 36.5 36.3 36.1 35.7 35.4 34.9 ...
    34.4 33.8 33.1 32.3 31.4 30.5 29.5 28.3 27.1 25.7 24.3 22.7 21.0 19.1 ...
    17.1 15.0 12.8 13.6 7.77 6.10]; 

%% Argon measured spectrum
% run CO2_cross_sections.m for Argon XUV_only data and save paramout_gauss as
% Ar_H11_param, etc. 

paramlist = [Ar_H11_param; Ar_H13_param; Ar_H15_param; Ar_H17_param; Ar_H19_param]; 
xdata = 0:0.05:20; 
ydata = Gauss(xdata, paramlist(1,1), paramlist(1,2), paramlist(1,3)) + ...
        Gauss(xdata, paramlist(2,1), paramlist(2,2), paramlist(2,3)) + ...
        Gauss(xdata, paramlist(3,1), paramlist(3,2), paramlist(3,3)) + ...
        Gauss(xdata, paramlist(4,1), paramlist(4,2), paramlist(4,3)) + ...
        Gauss(xdata, paramlist(5,1), paramlist(5,2), paramlist(5,3)); 
% xdata = E; 
% ydata = XUV_only'; 

%% deconvolve
IP = [15.763];
ygate = interp1(pics_x-IP, pics_y, xdata); 
ygate(isnan(ygate)) = 0; 
ygate = ygate / sum(ygate); 
ydata = ydata / sum(ydata); 
% [y1, y2] = deconv(ydata(4:end), ygate(4:end)); 
% y_APT = y2; 
y_APT = ifft( fft(ydata)./fft(ygate) ); 
y2 = y_APT; 

figure; plot(xdata, y2, 'DisplayName', 'APT spectrum'); 
figure; plot(xdata, abs(ygate.*y2), 'DisplayName', 'recontructed data'); 
axl = AddHarmonicAxis(gca, IP, wavelength);
for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.9    0.99];
    axl(i).FontSize = 12*0.75; 
    axl(i).XLabel.String = IP_label(i); 
end

%%
function yout = Gauss(x,A,mu,sig)
    yout = A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
end

