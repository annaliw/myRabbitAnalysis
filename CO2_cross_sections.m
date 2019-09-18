% fit region
% region = [1.6 3.3]; 
% region = [4.5 6.3]; 
region = [7.6 9.5];
% region = [10.2 12.5]; 
IP = [13.9000   17.6000   18.0770];

E = XUV_old_x; 
XUV_only = XUV_old_y'; 

tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 

[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), abs(XUV_only(start:stop)), XUV_only(start:stop), 1, 0); 
% save as labeled variables
paramout_original = paramout; 
fval_original = fval; 

% plot results
% Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
Yout = 0; 

set(groot,'defaultLineLineWidth',2.0)
figure; hold on; 
plot(E(start:stop), XUV_only(start:stop)/sum(XUV_only(start:stop)), 'bo', ...
    'DisplayName', 'XUV only data'); 
for ii = 1:size(paramout_gauss,1)
    Amp = paramout_gauss(ii,1); 
%             Amp = 1; 
    E0 = paramout_gauss(ii,2); 
    wid = paramout_gauss(ii,3); 
    
    ytmp = Gauss(E,Amp,E0,wid); 
    plot(E(start:stop), ytmp(start:stop), 'b--', 'HandleVisibility', 'off'); 
    Yout = Yout + ytmp;
end
plot(E(start:stop), Yout(start:stop), 'b-', 'DisplayName', 'full gaussian fit'); 
legend; 
xlabel('electron kinetic energy (eV)'); 
ylabel('amplitude'); 
axl = AddHarmonicAxis(gca, IP, wavelength);
for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.9    0.99];
    axl(i).FontSize = 12*0.75; 
    axl(i).XLabel.String = IP_label(i); 
end
xlim([E(start), E(stop)]);  
hold off; 

%% special fit for H13

% fit region
% region = [1.6 3.3]; harmonic = 13; 
region = [7.7 9.5]; harmonic = 17; 
tolerance = 0.05; 
% fit section set-up
start = find(abs(E-region(1))<tolerance, 1, 'last'); 
stop = find(abs(E-region(2))<tolerance, 1, 'first'); 
IP = [13.78 13.95 17.32 17.45 17.59 17.72 17.86 18.1370]-0.1;

% load CO2 vpeaks parameters from previous fit
% or run above with IP = [13.78 13.95 17.32 17.45 17.59 17.72 17.86 18.1370];
Amp = paramout_gauss(end,1); 
E0 = paramout_gauss(end,2); 
wid = paramout_gauss(end,3);
% Amp = 0.02; 
% E0 = 7.9885; 
% wid = 0.1649; 
bpeak = Gauss(E(start:stop),Amp,E0,wid); 

Amp = paramout_gauss(1,1); 
E0 = paramout_gauss(1,2); 
wid = paramout_gauss(1,3);
% Amp = 0.024; 
% E0 = 9.1978; 
% wid = 0.1854; 
xpeak = Gauss(E(start:stop),Amp,E0,wid); 

% LSQCURVEFIT FIT
guess = [0.015 0 0.08]; 
fun = @(x,xdata) bpeak + xpeak + acomb(harmonic, xdata, IP(3:7), x); 
[paramout, fval] = lsqcurvefit(fun, guess, E(start:stop), (XUV_only(start:stop)/sum(XUV_only(start:stop)))'); 

% plot results
apeak = acomb(harmonic, E(start:stop), IP(3:7), paramout); 

set(groot,'defaultLineLineWidth',2.0)
figure; hold on; 
plot(E(start:stop), XUV_only(start:stop)/sum(XUV_only(start:stop)), 'bo', ...
    'DisplayName', 'XUV only data'); 
plot(E(start:stop), bpeak, 'b--', 'HandleVisibility', 'off'); 
plot(E(start:stop), xpeak, 'b--', 'HandleVisibility', 'off'); 
plot(E(start:stop), apeak, 'b--', 'HandleVisibility', 'off'); 
plot(E(start:stop), apeak + bpeak + xpeak, 'b-', 'DisplayName', 'full gaussian fit'); 
legend; 
xlabel('electron kinetic energy (eV)'); 
ylabel('amplitude'); 
axl = AddHarmonicAxis(gca, IP, wavelength);
for i = 1:numel(IP)
    axl(i).XLabel.Position = [ -1.9    0.99];
    axl(i).FontSize = 12*0.75; 
    axl(i).XLabel.String = IP_label(i); 
end
xlim([E(start), E(stop)]);  
hold off; 


%%
figure; hold on; 
for ii=3:7
    plot(E(start:stop), Gauss(E(start:stop), 1, harmonic*1240/810-IP(ii), 0.08), 'b--'); 
end
plot(E(start:stop), acomb(E(start:stop), IP(3:7), [1 0 0.08]), 'b-'); 

function yout = acomb(harmonic, x, IPlist, p)
    yout = 0; 
    A = p(1); 
    x0 = p(2); 
    wid = p(3); 
    env = Gauss(x, A*1.5, sum(harmonic*1240/810 - IPlist)/5, 0.4); % check width for envelope
    for ii=1:numel(IPlist)
        tolerance = 0.05; 
        pos = (harmonic*1240/810 - IPlist(ii)) - x0; 
        posind = find(abs(x - pos) < tolerance, 1); 
        yout = yout + Gauss(x, env(posind), pos, wid); 
    end
end

function yout = Gauss(x,A,mu,sig)
    yout = A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
end



