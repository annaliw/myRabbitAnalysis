function axl = AddHarmonicAxis(ax1, Ip, IP_label, wavelength, yvis)

newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             ];
         
colororder(newcolors); 

POS = ax1.Position; %Position of Axis
del = 0.05*POS(4);

% n=9:1:19; %Harmonics to plot
lambda = wavelength; % Fundamental Wavelength
% lambda = 2*810/1.778; 

% ntrue = 1:1:21;
ntrue = 2:2:20; 
%Harmonic Labels
for i=1:numel(ntrue)
    nlbl{i} = sprintf('%2.0f',ntrue(i));
end

for i=1:numel(Ip)
    
    axl(i) = axes('Position',POS);
    axl(i).FontSize = 14; 
%     axl(i).XColor = ax1.ColorOrder(i,:);
    axl(i).Color = 'none';
    axl(i).XAxisLocation = 'top';
%     axl(i).Position(2) = POS(2)+POS(4); 
    axl(i).XTick = ntrue*(1240/lambda)-Ip(i);
    ylim(axl(i), [0 1]); axl(i).YTick = 1; 
    axl(i).XGrid = 'on';
    axl(i).GridAlpha = 1; 
    %xlim(axl(i), [0,20]); 
    if yvis == 0
        axl(i).YTickLabel = '';
%         ylabel = ''; 
        axl(i).XTickLabel = ''; 
    elseif yvis == 1
        axl(i).XTickLabel = nlbl;
        axl(i).YTickLabel = IP_label(i);
%         ylabel = IP_label(i); 
    end
    
    if yvis==1
        POS(4) = POS(4)-del;  
    end
end
ax1.Position = POS;

linkaxes([ax1,axl],'x'); 
end