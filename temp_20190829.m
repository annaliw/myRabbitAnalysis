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

pf = ((2.*K).^(1j*Z./K)) ./ ((2.*k).^(1j*Z./k)); %Pre-factor
d = ((K-k).^(1j*Z.*(1./K-1./k))); %Denominator

Parg = pf .* igamma( 1 + l + 1j*Z .* (1./K - 1./k), 0 ) ./ d ; %Arguemnt of angle
%igamma is upper incomplete gamma function
ePP = angle( Parg ); %C-C phase for emission in asymptotic approx. 

alpha = angle(1 + 1j*Z/2 * (1./K.^2 + 1./k.^2) .* (K-k)./(1 + 1j*Z*(1./K - 1./k))); 
ePPA = ePP + alpha; 


%Momentum of electron before C-C transition (absorption)
K = sqrt( k.^2 - 2*w );

pf = ((2.*K).^(1j*Z./K)) ./ ((2.*k).^(1j*Z./k)); %prefactor
d = ((K-k).^(1j*Z.*(1./K-1./k))); %denominator

Parg = pf .* igamma( 1 + l + 1j*Z .* (1./K - 1./k), 0 ) ./ d ;
aPP = angle( Parg ); %C-C phase for absorption in asymptoic approx. 

alpha = angle(1 + 1j*Z/2 * (1./K.^2 + 1./k.^2) .* (K-k)./(1 + 1j*Z*(1./K - 1./k))); 
aPPA = aPP + alpha; 


