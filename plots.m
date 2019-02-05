% % Load saved figures
% a=hgload('NO_plot1.fig');
% b=hgload('NO_plot2.fig');
% c=hgload('NO_plot3.fig');
% Prepare subplots
figure; hold on; 

ax1=subplot(3,1,1);
x1 = E; 
y1 = mean(abs(E_SpectraArray_early),2); 
plot(x1, y1, 'k'); 
set(ax1, 'XTick', []); 
set(ax1, 'position', [0.1, 0.65, 0.7, 0.3] )


ax2=subplot(3,1,2);
x2 = E; 
y2 = twoOmega_abs; 
plot(x2, y2, 'b'); 
set(ax2, 'XTick', [], 'position', [0.1, 0.35, 0.7, 0.3] ); 

ax3=subplot(3,1,3);
x3 = E; 
y3 = stageTimes*10^(15); 
imagesc(x3, y3, (abs(Espectra_TD)-mean(abs(Espectra_TD),2)).');
ylim(ax3, [0 10]); 
ylabel('XUV and IR pulse delay (fs)'); 
xlabel('Energy (eV)'); 
set(ax3, 'position', [0.1, 0.05, 0.7, 0.3] )




