ESpectra_norm = E_SpectraArray ./ repmat(sum(E_SpectraArray,1), size(E_SpectraArray,1), 1, 1);
XUV_norm = XUV_only ./ sum(XUV_only);
tmp = sum(abs(ESpectra_norm - repmat(XUV_norm, 1, 223, 37)),3);

% make color map
numc = 50; 
cm = ones([numc 3]); 
cm(1:(numc/2),1:2) = repmat((1:1:(numc/2))' ./ (numc/2), 1, 2); 
cm((numc/2+1):end,2:3) = flipud(cm(1:(numc/2),1:2)); 

figure; hold on; imagesc(E, stageTimes*1e15, (sum(abs(tmp'),3)-mean(sum(abs(tmp'),3),1)));
colormap(cm); 
cb = colorbar; 
cb.Label.String = 'XUV/IR - XUV only'; 
cb.Ticks = []; 
caxis([-0.01 0.01]); 
xlabel('electron kinetic energy (eV)'); 
ylabel('XUV/IR delay (fs)'); 
goodplot(18); 