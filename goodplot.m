function goodplot(fontsize)
% function which produces a nice-looking plot
% and sets up the page for nice printing
% axis square;
box on; 
set(gca,'LineWidth',2);
set(findall(gcf, 'Type', 'Line'), 'MarkerSize', 10, 'LineWidth', 2); 
set(findall(gcf, 'Type', 'ErrorBar'), 'MarkerSize', 10, 'LineWidth', 2); 
set(gca,'FontSize',fontsize);
set(gca,'FontWeight','Bold');
set(get(gca,'xlabel'),'FontSize', fontsize, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fontsize, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fontsize, 'FontWeight', 'Bold');
set(findall(gcf,'Type','Text'),'FontSize',fontsize, 'FontWeight', 'Bold'); 
% box off; 

set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [8 8]);
set(gcf,'PaperPosition',[0.5 0.5 7 7]);
set(gcf,'PaperPositionMode','Manual');
