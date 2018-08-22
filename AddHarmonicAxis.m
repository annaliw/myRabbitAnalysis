function axl = AddHarmonicAxis(fh,Ip)

ax1 = fh.Children; %Axis Handle
POS = ax1.Position; %Position of Axis
del = 0.07*POS(4);

n=9:1:19; %Harmonics to plot
lambda = 810; % Fundamental Wavelength
% lambda = 2*810/1.778; 

%Harmonic Labels
for i=1:numel(n)
    nlbl{i} = sprintf('%2.0f',n(i));
end

for i=1:numel(Ip)
    
    axl(i) = axes('Position',POS);
    axl(i).XColor = ax1.ColorOrder(i,:);
    axl(i).Color = 'none';
    axl(i).XAxisLocation = 'top';
    axl(i).XTick = n*(1240/lambda)-Ip(i);
    axl(i).XTickLabel = nlbl;
    axl(i).YTick = 0;
    axl(i).YTickLabel = '';
    axl(i).XGrid = 'on';
    %xlim(axl(i), [0,20]);
    
    POS(4) = POS(4)-del;
end
ax1.Position = POS;

linkaxes([ax1,axl],'x')
end