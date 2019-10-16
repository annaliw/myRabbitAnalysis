function goodplot()
% function which produces a nice-looking plot
% and sets up the page for nice printing
set(get(gca,'xlabel'),'FontSize', 14, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', 14, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', 16, 'FontWeight', 'Bold');
box off; axis square;
set(gca,'LineWidth',2);
set(get(gca, 'Children'), 'MarkerSize', 8); 
set(get(gca, 'Children'), 'LineWidth', 1.5); 
% set(gca,'FontSize',16);
% set(gca,'FontWeight','Bold');

set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [8 8]);
set(gcf,'PaperPosition',[0.5 0.5 7 7]);
set(gcf,'PaperPositionMode','Manual');
