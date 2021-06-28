%% check what RABBITT spectrum might look like

IP = [9.538 9.538+0.626 12.13 15]; % X3/2, X1/2, A
HH = 3:2:21; 

% sigma = 0.06; 
XUV_X1 = ones([1 numel(Ebins)]); 
XUV_X2 = ones([1 numel(Ebins)]); 
XUV_A  = ones([1 numel(Ebins)]); 
XUV_B  = ones([1 numel(Ebins)]); 
XUV_X1e = ones([1 numel(Ebins)]); 
XUV_X2e = ones([1 numel(Ebins)]); 
SB_X1 = ones([1 numel(Ebins)]); 
SB_X2 = ones([1 numel(Ebins)]); 
SB_A  = ones([1 numel(Ebins)]); 
SB_B  = ones([1 numel(Ebins)]); 
SB_X1e = ones([1 numel(Ebins)]); 
SB_X2e = ones([1 numel(Ebins)]); 
for ii=1:numel(HH)
    XUV_X1 = XUV_X1 + 2*exp(-(Ebins-(HH(ii)*1240/810-9.538)).^2/(2*0.06.^2)); 
    SB_X1  = SB_X1 + 2*0.2*exp(-(Ebins-((HH(ii)-1)*1240/810-9.538)).^2/(2*0.06.^2)); 
    
    XUV_X2 = XUV_X2 + 2*exp(-(Ebins-(HH(ii)*1240/810-10.168)).^2/(2*0.06.^2)); 
    SB_X2  = SB_X2 + 2*0.2*exp(-(Ebins-((HH(ii)-1)*1240/810-10.168)).^2/(2*0.06.^2)); 
    
    XUV_A  = XUV_A + 0.5*exp(-(Ebins-(HH(ii)*1240/810-12.13)).^2/(2*0.43.^2)); 
    SB_A   = SB_A + 0.5*0.2*exp(-(Ebins-((HH(ii)-1)*1240/810-12.13)).^2/(2*0.43.^2)); 
    
    XUV_B  = XUV_B + 0.6*exp(-(Ebins-(HH(ii)*1240/810-15)).^2/(2*0.9.^2)); 
    SB_B   = SB_A + 0.6*0.2*exp(-(Ebins-((HH(ii)-1)*1240/810-15)).^2/(2*0.9.^2)); 
    
    XUV_X1e = XUV_X1e + 2*exp(-(Ebins-(HH(ii)*1240/810-15.73+1)).^2/(2*0.08.^2)); 
    SB_X1e  = SB_X1e + 2*0.2*exp(-(Ebins-((HH(ii)-1)*1240/810-15.73+1)).^2/(2*0.08.^2)); 
    
    XUV_X2e = XUV_X2e + 2*exp(-(Ebins-(HH(ii)*1240/810-16.37+1)).^2/(2*0.08.^2)); 
    SB_X2e  = SB_X2e + 2*0.2*exp(-(Ebins-((HH(ii)-1)*1240/810-16.37+1)).^2/(2*0.08.^2)); 
end

figure; hold on;  
subplot(2,1,1); hold on; 
ax = gca;
plot(Ebins, (XUV_X1 + XUV_X2 + XUV_A + XUV_B+SB_X1 + SB_X2 + SB_A + SB_B) ./ ...
    sum(XUV_X1 + XUV_X2 + XUV_A + XUV_B+SB_X1 + SB_X2 + SB_A + SB_B), 'k-'); 
plot(Ebins, (XUV_X1 + XUV_X2 + XUV_A + XUV_B) ./ sum(XUV_X1 + XUV_X2 + XUV_A + XUV_B), 'b-'); 
plot(Ebins, (SB_X1 + SB_X2 + SB_A + SB_B) ./ sum(SB_X1 + SB_X2 + SB_A + SB_B), 'r-'); 
axh = AddHarmonicAxis(ax, IP, ["X3/2", "X1/2", "A", "B"], 810, 1); 

subplot(2,1,2); hold on; 
ax = gca;
plot(Ebins, (XUV_X1 + XUV_X2 + XUV_A + XUV_B+SB_X1 + SB_X2 + SB_A + SB_B + XUV_X1e + XUV_X2e + SB_X1e + SB_X2e) ./ ...
    sum(XUV_X1 + XUV_X2 + XUV_A + XUV_B+SB_X1 + SB_X2 + SB_A + SB_B + XUV_X1e + XUV_X2e + SB_X1e + SB_X2e), 'k-'); 
plot(Ebins, (XUV_X1 + XUV_X2 + XUV_A + XUV_B + XUV_X1e + XUV_X2e) ./ sum(XUV_X1 + XUV_X2 + XUV_A + XUV_B + XUV_X1e + XUV_X2e), 'b-'); 
plot(Ebins, (SB_X1 + SB_X2 + SB_A + SB_B + SB_X1e + SB_X2e) ./ sum(SB_X1 + SB_X2 + SB_A + SB_B + SB_X1e + SB_X2e), 'r-'); 
axh = AddHarmonicAxis(ax, IP, ["X3/2", "X1/2", "A", "B"], 810, 1); 

%%