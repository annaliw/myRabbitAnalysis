[xdata, ind] = sort([H2TDSE_810_140.x_Ee; H2TDSE_810_145.x_Ee; H2TDSE_810_150.x_Ee]); 
ydata = [H2TDSE_810_140.t; H2TDSE_810_145.t; H2TDSE_810_150.t]; 
ydata = ydata(ind); 

myfun = @(x, xdata) 1./xdata .* (x(1).*log(xdata) + x(2)); 
x0 = [1 1]; 
x = lsqcurvefit(myfun, x0, xdata, ydata); 

figure; hold on; 
plot(xdata, ydata, 'ro-'); 
xtmp = interp1(1:numel(xdata), xdata, 1:0.1:numel(xdata)); 
plot(xtmp, myfun(x, xtmp), 'k-'); 











% [tmp1, ~] = coulombScatteringPhase(1,1,E+1240/810); 
% [tmp2, ~] = coulombScatteringPhase(1,1,E-1240/810); 
% cdelay = 24.2 * (tmp1 - tmp2) / (2*1240/810); 
% % plot(E, Ivanov_CC, 'b--', 'DisplayName', 'Ivanov CLC'); 
% % plot(E(2:end), delay01 + CCPAp_tmp(2:end), 'k-'); 
% % plot(E, cdelay + CCPAp_tmp, 'k-'); 
% % plot(E, cdelay + Ivanov_CC, 'k--'); 
% % plot(E, cdelay + CCP_tmp, 'k-.'); 