%  a*E = me/2 * L^2/(t-t0)^2
% E = n*hw - IP

% bohr radius 5.29177210903e−11
% l in au = l in m * 1/bohr radius

t0 = 25 * 1e-9; % convert to s
xdata = tof_peak(:).' * 1e-9; % convert to s
ydata = calibEnergy(:).' * 1.60218e-19; % convert to J

% x = [a L E0]
me = 9.10938356e-31; 
ydata = 1./ydata; 
fun = @(x,xdata) 1 ./ ((me/2 * x^2./(xdata-t0).^2)) ; 

x0 = 1; 
[paramout, resnorm] = lsqcurvefit(fun, x0, xdata, ydata); 

figure; hold on; 
plot(xdata, ydata, 'bo'); 
plot(xdata, fun(paramout, xdata), 'r-'); 