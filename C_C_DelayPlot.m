E_AU = 27.2114; %eV to a.u.
T_AU = 0.024189; %fs to a.u. 

w = 1239.8 / 810 / E_AU; %Photon Energy
T_L = 2*pi/w * T_AU; % laser period in fs

Z = 1; %Atomic Charge
l = 0; %Final state angular momentum

%Final Momentum of outgoing electron
% k = 0:0.001:3; 
k = logspace(-2, 0, 1500); 

%Final Energy of outgoing electron
E = k.^2/2 .* E_AU; 

%Momentum of electron before C-C transition (emission)
K = sqrt( k.^2 + 2*w );

pf = ((2.*K).^(1i*Z./K)) ./ ((2.*k).^(1i*Z./k)); %Pre-factor
d = ((K-k).^(1i*Z.*(1./K-1./k))); %Denominator

Parg = pf .* igamma( 1 + l + 1i*Z .* (1./K - 1./k), 0 ) ./ d ; %Arguemnt of angle
%igamma is upper incomplete gamma function
ePP = angle( Parg ); %C-C phase for emission in asymptotic approx. 

gkK = 1i*Z .* ((K-k).*(K.^2+k.^2))./(2.*K.^2.*k.^2) .* igamma(1 + 1i*Z.*(1./K-1./k),0);
ePPA = angle ( Parg + pf .* gkK ./ d); %C-C phase for emission in long-range approximation     

% sigma = angle(igamma(l + 1 - 1i*Z./K, 0)); 
% sigma = fliplr(unwrap(fliplr(sigma))); 
% ak = 2*exp(-2*sigma.*K); 
r0 = 1i*Z/4 .* (1./K.^2 + 1./k.^2);
% r0 = 1./(K.*ak); 
igkK = 1i*Z .* ((K-k).*(K.^2+k.^2))./(2.*K.^2.*k.^2) .* igamma(1 + 1i*Z.*(1./K-1./k),r0);

ePPAp = angle ( Parg + pf .* igkK ./ d); %C-C phase for emission in modified long-range approximation

%Momentum of electron before C-C transition (absorption)
K = sqrt( k.^2 - 2*w );

pf = ((2.*K).^(1i*Z./K)) ./ ((2.*k).^(1i*Z./k)); %prefactor
d = ((K-k).^(1i*Z.*(1./K-1./k))); %denominator

Parg = pf .* igamma( 1 + l + 1i*Z .* (1./K - 1./k), 0 ) ./ d ;
aPP = angle( Parg ); %C-C phase for absorption in asymptoic approx. 

gkK = 1i*Z .* ((K-k).*(K.^2+k.^2))./(2.*K.^2.*k.^2) .* igamma(1 + 1i*Z.*(1./K-1./k),0);
aPPA = angle ( Parg + pf .* gkK ./ d);  %C-C phase for absorption in long-range approximation     

% sigma = angle(igamma(l + 1 - 1j*Z./K, 0)); 
% sigma = fliplr(unwrap(fliplr(sigma))); 
% ak = 2*exp(-2*sigma.*K); 
r0 = 1i*Z/4 .* (1./K.^2 + 1./k.^2);
% r0 = 1./(K.*ak); 
igkK = 1i*Z .* ((K-k).*(K.^2+k.^2))./(2.*K.^2.*k.^2) .* igamma(1 + 1i*Z.*(1./K-1./k),r0);

aPPAp = angle ( Parg + pf .* igkK ./ d); %C-C phase for emission in modified long-range approximation

%%

h_300 = figure(300);

plot(fliplr(E),unwrap(fliplr(ePP)),'LineStyle','-','LineWidth',2,'DisplayName','P')
ind = h_300.Children.ColorOrderIndex - 1;
hold on
h_300.Children.ColorOrderIndex = ind;
plot(fliplr(E),unwrap(fliplr(ePPA)),'LineStyle','--','LineWidth',2,'DisplayName','P+A')
h_300.Children.ColorOrderIndex = ind;
plot(fliplr(E),unwrap(fliplr(ePPAp)),'LineStyle','-.','LineWidth',2,'DisplayName','P+A^{\prime}')

ind = h_300.Children.ColorOrderIndex;
plot(fliplr(E),unwrap(fliplr(aPP)),'LineStyle','-','LineWidth',2,'DisplayName','P')
h_300.Children.ColorOrderIndex = ind;
plot(fliplr(E),unwrap(fliplr(aPPA)),'LineStyle','--','LineWidth',2,'DisplayName','P+A')
h_300.Children.ColorOrderIndex = ind;
plot(fliplr(E),unwrap(fliplr(aPPAp)),'LineStyle','-.','LineWidth',2,'DisplayName','P+A^{\prime}')


%hold off
xlim([0,30])
xlabel('Energy (eV)')
ylabel('Phase (rad)')
title('Continuum-Continuum Transition Phase')
h_300.Children.FontSize = 14;
h_300.Children.LineWidth = 2;
legend("ePP", "ePPA", "ePPAp", "aPP", "aPPA", "aPPAp")

CCP = unwrap(fliplr(ePP)) - unwrap(fliplr(aPP));
CCPA = unwrap(fliplr(ePPA)) - unwrap(fliplr(aPPA));
CCPAp = unwrap(fliplr(ePPAp)) - unwrap(fliplr(aPPAp));

h_301 = figure(301);
ax1 = subplot(2,1,1);
plot(fliplr(E), CCP.*(T_L*1000/2/(2*pi)), 'LineStyle','-','LineWidth',2,'DisplayName','P')
ind1 = ax1.ColorOrderIndex - 1;
hold on
ax1.ColorOrderIndex = ind1;
plot(fliplr(E), CCPA.*(T_L*1000/2/(2*pi)), 'LineStyle','--','LineWidth',2,'DisplayName','P+A')
ax1.ColorOrderIndex = ind1;
plot(fliplr(E), CCPAp.*(T_L*1000/2/(2*pi)), 'LineStyle','-.','LineWidth',2,'DisplayName','P+A^{\prime}')
%hold off
xlim([0,30])
xlabel('Energy (eV)')
ylabel('Delay (as)')
title('Continuum-Continuum Delays')


ax2 = subplot(2,1,2);
plot(fliplr(E(1:end-1)), diff(CCP)./diff(E).*(1351/(2*pi)),'LineStyle','-','LineWidth',2)
ind2 = ax2.ColorOrderIndex - 1;
hold on
ax2.ColorOrderIndex = ind2;
plot(fliplr(E(1:end-1)), diff(CCPA)./diff(E).*(1351/(2*pi)),'LineStyle','--','LineWidth',2)
ax2.ColorOrderIndex = ind2;
plot(fliplr(E(1:end-1)), diff(CCPAp)./diff(E).*(1351/(2*pi)),'LineStyle','-.','LineWidth',2)
%hold off
xlabel('Energy (eV)')
ylabel('Slope of Delay (as/eV)')
title('Slope of Continuum-Continuum Delays')

xlim([0,10])

h_301.Children(1).FontSize = 14;
h_301.Children(1).LineWidth = 2;
h_301.Children(2).FontSize = 14;
h_301.Children(2).LineWidth = 2;

%% Method Comparison
%Look at the difference in the models for correcting small energy
%differences. 

deltaE = 0.5;

Ed = E + deltaE;

CCPd = interp1(fliplr(E),CCP,fliplr(Ed));
CCPAd = interp1(fliplr(E),CCPA,fliplr(Ed));
CCPApd = interp1(fliplr(E),CCPAp,fliplr(Ed));

dCCP = CCPd - CCP;
dCCPA = CCPAd - CCPA;
dCCPAp = CCPApd - CCPAp;

h_302 = figure(302);

plot(fliplr(E),dCCP.*(T_L*1000/2/(2*pi)),'LineStyle','-','LineWidth',2,'DisplayName','P')
ind = h_302.Children.ColorOrderIndex - 1;
h_302.Children.ColorOrderIndex = ind;
hold on
plot(fliplr(E),dCCPA.*(T_L*1000/2/(2*pi)),'LineStyle','--','LineWidth',2,'DisplayName','P+A')
h_302.Children.ColorOrderIndex = ind;
plot(fliplr(E),dCCPAp.*(T_L*1000/2/(2*pi)),'LineStyle','-.','LineWidth',2,'DisplayName','P+A^{\prime}')
hold off
xlim([2,10])
xlabel('Energy (eV)')
ylabel('\tau_{CC}(E+\DeltaE) - \tau_{CC}(E) (as)')
title(['Differential Delay for \Delta',sprintf('E = %1.2f eV', deltaE)])
h_302.Children.FontSize = 14;
h_302.Children.LineWidth = 2;