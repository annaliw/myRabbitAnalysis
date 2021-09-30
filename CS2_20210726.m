IP = [10.07 12.689 13.854 14.068 14.474 16.183  16.7 16.9]; % X3/2, X1/2, A
IP_label = ["X", "A", "Pi", "sat", "B", "C", "D1", "D2"]; 
HH = 3:2:21; 

XUV = ones([1 numel(Ebins)]); 
SB  = ones([1 numel(Ebins)]); 
for ii=1:numel(HH)
    for jj=1:numel(IP)
        XUV = XUV + 2*exp(-(Ebins-(HH(ii)*1240/810-IP(jj))).^2/(2*0.06.^2));
        SB  = SB  + 0.5*exp(-(Ebins-((HH(ii)+1)*1240/810-IP(jj))).^2/(2*0.06.^2));
    end
end

figure; hold on; 
subplot(2,1,1); hold on; 
ax1 = gca; 
plot(Ebins, XUV, 'b-'); 
plot(Ebins, SB, 'r-'); 
xlabel('photoelectron energy (eV)'); 
ylabel('yield')
AddHarmonicAxis(ax1, IP, IP_label, 810, 1); 
subplot(2,1,2); hold on; 
plot(Ebins, XUV+SB, 'k-'); 
xlabel('photoelectron energy (eV)'); 
ylabel('yield'); 
AddHarmonicAxis(ax1, IP, IP_label, 810, 1); 

