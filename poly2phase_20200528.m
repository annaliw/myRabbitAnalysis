% refit SB gaussians

% Argon regions
% region = [2.44 2.8]; % sideband 12
region = [5.4 5.85]; % sideband 14
% region = [8.2 9.2]; % sideband 16
% region = [11 12.3]; % sideband 18

% H2 regions
% region = [1.45 3.1]; % sideband 12
% region = [4.55 6.15]; % sideband 14
% region = [7.7548 9.18]; % sideband 16

[~,start] = min(abs(E-region(1))); 
[~,stop]  = min(abs(E-region(2))); 
allpks = repmat((9:1:19)', [1 numel(IP)]) .* 1240/810 - repmat(IP, [numel(9:1:19) 1]); 
allpks = allpks(:); 
regpks = allpks; 
for ii=1:numel(allpks)
    if regpks(ii) < E(start)
        regpks(ii) = NaN; 
    elseif regpks(ii) > E(stop)
        regpks(ii) = NaN; 
    end
end
regpks(isnan(regpks)) = []; 

xdata = E(start:stop); 
ydata = twoOmega_signal(start:stop).' ./ sum(twoOmega_signal(start:stop)); 

% sig = 0.3; % starting value for width
% fun = @(x,xdata) sum(x(:,1).*exp(-(xdata-x(:,2)).^2./(2.*x(:,3)).^2), 1); 
fun = @(x,xdata) gaussian_mixture_fixsigma(x, xdata); 
x0 = zeros([numel(regpks), 3]); 
x0(:,1) = max(abs(ydata)); 
x0(:,2) = regpks; 
x0(:,3) = 0.1; 
[paramout, resnorm] = lsqcurvefit(fun, x0, xdata, abs(ydata)); 

figure; hold on; 
plot(xdata, abs(ydata), 'bo'); 
% plot(xdata, gaussian_mixture_fixsigma(paramout, xdata), 'k-');
plot(xdata, fun(paramout, xdata), 'k-');
for ii=1:size(paramout,1)
    plot(xdata, paramout(ii,1).*exp(-(xdata-paramout(ii,2)).^2./(2.*paramout(ii,3)).^2), 'k--'); 
end


% % fit only width? 
% amps = paramout_centers(:,1); 
% centers = paramout_centers(:,2); 
% fun = @(x,xdata) sum(amps.*exp(-(xdata-centers).^2./(2.*x).^2), 1); 
% x0 = 0.3
% paramout_width = lsqcurvefit(fun, x0, xdata, abs(ydata));  
% 
% figure; hold on; 
% plot(xdata, abs(ydata), 'bo'); 
% plot(xdata, fun(paramout_width, xdata), 'k-'); 
% plot(xdata, amps(1).*exp(-(xdata-centers(1)).^2./(2.*paramout_width(1)).^2), 'k--'); 
% plot(xdata, amps(2).*exp(-(xdata-centers(2)).^2./(2.*paramout_width(1)).^2), 'k--'); 



%%
% choose polynomial degree before fitting
pdeg = 1; 

% xdata = E(start:stop); 
% ydata = tmp_signal(start:stop).'; 

figure; hold on; 
yyaxis left; ylabel('amplitude'); 
plot(xdata, abs(ydata), 'o'); 
yyaxis right; ylabel('phase'); 
plot(xdata, angle(ydata), 'x'); 

y_amp = zeros([size(paramout_gauss,1) numel(xdata)]); 
poly = zeros([size(paramout_gauss,1) pdeg+1]); 
yfit = zeros([1 numel(xdata)]); 
for ii=1:size(paramout_gauss,1)
% for ii=1:2
    y_amp(ii,:) = paramout_gauss(ii,1) .* exp(-(xdata-paramout_gauss(ii,2)).^2/(2*paramout_gauss(ii,3).^2)); 
    yyaxis left; plot(xdata, y_amp(ii,:), 'b--'); 
    
    window = 4; 
    [cval,center] = min(abs(xdata - paramout_gauss(ii,2))); 
    if center-window < 1
        ind1 = 1; 
    else
        ind1 = center-window;
    end
    if center+window > numel(xdata)
        ind2 = numel(xdata); 
    else
        ind2 = center+window; 
    end
    poly(ii,:) = polyfit(xdata(ind1:ind2)-xdata(center), angle(ydata(ind1:ind2)), pdeg); 
    yyaxis right; plot(xdata(ind1:ind2), polyval(poly(ii,:), xdata(ind1:ind2)-xdata(center)), 'k--'); 
    
    yfit = yfit + y_amp(ii,:) .* exp(1j * polyval(poly(ii,:), xdata-xdata(center))); 
end

yyaxis left; plot(xdata, abs(yfit), 'b-'); 
yyaxis right; plot(xdata, angle(yfit), 'r-'); 

xlabel('electron energy (eV)'); 
goodplot(24)



