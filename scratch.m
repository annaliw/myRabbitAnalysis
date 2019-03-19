%% plot Coulomb scattering phase

hw = 1240/wavelength; 

% k = sqrt(2*((11:1:21)*hw-IP(1))/27.211); 
k = logspace(-1, 0, 1000); 
E = k.^2/2; 
l = 0; 
sigma = unwrap(angle(igamma(l + 1 - 1j./k, 0))); 
delay = 24.2 * diff(sigma)./diff(E); 

figure; 
plot(E, sigma); 
figure; 
plot(27.221*E(2:end), delay, 'o-'); 

%% subtract Coulomb scattering phase from H2 data

k = [12*hw - IP, 14*hw - IP, 16*hw - IP]; 
sigma = angle(igamma(l + 1 - 1j./k, 0)); 
measured_phases = [sb12_param(:,2)', sb14_param(:,2)', sb16_param(:,2)']; 
sub_phases = [sb12_param(:,2)'-sb12_param(1,2), sb14_param(:,2)'-sb14_param(1,2), sb16_param(:)'-sb16_param(1,2)]; 
cc_phase = sub_phases - sigma; 

figure; plot(k, sigma, 'o'); 
figure; plot(k, cc_phase, 'o'); 