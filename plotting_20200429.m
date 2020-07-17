signal = squeeze(sum(twoOmega_signal,2)); 
gauss = @(x,A,mu,sig) A .* exp( -(x-mu).^2 ./ (2.*sig.^2) ); 

SB=12; 
region = [1.6 2.8]; % sideband 12
ar_slope = Ar_SB12_slope(1); 
[~,start] = min(abs(E-region(1))); 
[~,stop]  = min(abs(E-region(2))); 
[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), ...
    abs(signal(start:stop)), signal(start:stop), 1, 0); 

% fh = figure('Position', [10 10 400*2.5 700]); 
figure; 
hold on; 
ylabel('phase'); 
for ii=1:size(paramout,1)
% for ii=2
   [~,center] = min(abs(E-paramout_gauss(ii,2))); 
   window = 4; 
   x = E((center-window):(center+window)); 
   y = signal((center-window):(center+window)).'; 
   t = RABBITT_phase((center-window):(center+window)); 
   
   plot(x-mean(x)-ii*0.3, t-mean(t) + 0.2, 'b-'); 
   
   peak = gauss(paramout_gauss(ii,2),paramout_gauss(ii,1),paramout_gauss(ii,2),paramout_gauss(ii,3)); 
   for jj=1:numel(x)
       if abs(y(jj)) >= 0.75 * peak 
           tmp = angle(y(jj))-ar_slope*x(jj) - mean(angle(y)-ar_slope*x) + 0.2; 
           p1 = line(x(jj)-mean(x)-ii*0.3, tmp, ...
               'Marker', 's', 'LineStyle', '-', 'Color', 'b', 'MarkerFaceColor', 'b'); 
%             p1 = plot(ax1, x(jj)-mean(x)-ii*0.3, angle(y(jj))-mean(paramout(:,1)), ...
%                 '-x', 'MarkerEdgeColor', 'b');
       end
   end   
end

ax1 = gca; 
ax1.XTick = []; 
% xlabel('(relative) photoelectron energy (eV)'); 
goodplot(30); 

SB=14; 
region = [4.75 5.9]; % sideband 14
ar_slope = Ar_SB14_slope(1); 
[~,start] = min(abs(E-region(1))); 
[~,stop]  = min(abs(E-region(2))); 
[paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, E(start:stop), ...
    abs(signal(start:stop)), signal(start:stop), 1, 0); 

for ii=1:size(paramout,1)
% for ii=2
   [~,center] = min(abs(E-paramout_gauss(ii,2))); 
   window = 4; 
   x = E((center-window):(center+window)); 
   y = signal((center-window):(center+window)).'; 
   t = RABBITT_phase((center-window):(center+window)); 
   
   plot(x-mean(x)-ii*0.3, t-mean(t), 'r-'); 
   
   peak = gauss(paramout_gauss(ii,2),paramout_gauss(ii,1),paramout_gauss(ii,2),paramout_gauss(ii,3)); 
   for jj=1:numel(x)
       if abs(y(jj)) >= 0.85 * peak 
           tmp = angle(y(jj))-ar_slope*x(jj) - mean(angle(y)-ar_slope*x); 
           p2 = line(x(jj)-mean(x)-ii*0.3, tmp, ...
               'Parent', ax1, 'Marker', 'o', 'LineStyle', '-', 'Color', 'r', 'MarkerFaceColor', 'r'); 
%             p2 = plot(ax1, x(jj)-mean(x)-ii*0.3, angle(y(jj))-mean(paramout(:,1)), ...
%                 '-o', 'MarkerEdgeColor', 'r', 'IconDisplayStyle', 'off'); 
       end
   end   
end
legend([p1 p2], 'SB12', 'SB14'); 
goodplot(30); 

ax2 = axes('Position', ax1.Position, ...
            'YAxisLocation', 'right', ...
            'Color', 'none'); 
ax2.YGrid = 'off'; 
ax2.XGrid = 'on'; 
% line(E(start:stop)-mean(E(start:stop)), abs(signal(start:stop)), 'Parent', ax2, 'Color', 'none'); 
line(fliplr(-(1:numel(paramout_gauss(:,2)))*0.3), fliplr(-(1:numel(paramout_gauss(:,2)))*0.3), ...
    'Parent', ax2, 'Color', 'none'); 
 
% xticks(flipud(paramout_gauss(:,2)))
xticks(fliplr(-(1:numel(paramout_gauss(:,2)))*0.3)); 
xticklabels(flipud(split(num2str(1:size(paramout,1)))))
xlabel('\nu'); 
ax2.YTick = []; 

% linkaxes([ax1, ax2], 'x')
goodplot(30)

%%
function Yout = Spectrum(E, gaussian, p)
        Yout = 0; 
        Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        if size(p,2) == 1
            Phase = @(x,b,mu) exp(1j .* b); 
            % sum the 2w signal
            for n = 1:size(p,1)
                Amp = gaussian(n,1); 
    %             Amp = 1; 
                E0 = gaussian(n,2); 
                wid = gaussian(n,3); 

                b = mod(p(n,1),2*pi); 

                Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,b,E0);
            end
        elseif size(p,2) == 2
            Phase = @(x,b,c,mu) exp(1j .* (b + c.*(x-mu)) ); 
            % sum the 2w signal
%             clist = [0.17, 0.64, 0.64]; 
            for n = 1:size(p,1)
                Amp = gaussian(n,1); 
    %             Amp = 1; 
                E0 = gaussian(n,2); 
                wid = gaussian(n,3); 

                b = mod(p(n,1),2*pi);
                c = p(n,2); 
%                 c = clist(n); 


                Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,b,c,E0);
            end
        else
            error('invalid guess input'); 
        end
        
        yout_mat = [abs(Yout); angle(Yout)]; 

    end
