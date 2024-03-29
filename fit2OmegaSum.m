function [paramout, fval] = fit2OmegaSum(xin, yin, gaussian, guess, plotting)

    global IP; 
    global IP_label; 
%     global wavelength; 
    wavelength = 810; 
    

    % FMINSEARCH FIT
    yin_abs = yin(1,:); 
    yin_phi = unwrap(yin(2,:)); 
%     chi2 = @(x) sum((real(yin)-real(Spectrum(xin, gaussian, x))).^2)... 
%         + sum((imag(yin)-imag(Spectrum(xin, gaussian, x))).^2); 
    chi2 = @(x) sum((yin_abs.*exp(1j*yin_phi) - Spectrum(xin, gaussian, x))...
        .*conj(yin_abs.*exp(1j*yin_phi) - Spectrum(xin, gaussian, x))); 
    [paramout, fval] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));

%     % LSQCURVEFIT FIT
%     yin_abs = yin(1,:); 
%     yin_phi = unwrap(mod(yin(2,:), 2*pi)); 
%     fun = @(x,xdata) Spectrum(xdata, gaussian, x); 
%     [paramout, fval] = lsqcurvefit(fun, guess, xin, yin); 


%     figure; hold on; 
%     yyaxis left
%     scatter(xin, yin_abs, 'o'); 
%     plot(xin, abs(Spectrum(xin, gaussian, paramout)));
%     ylabel('2\omega Amplitude'); 
%     yyaxis right
%     scatter(xin, yin_phi, '+'); 
%     plot(xin, mod(angle(Spectrum(xin, gaussian, paramout)),2*pi)); 
%     ylabel('2\omega Phase'); 
%     ylim([0 2*pi]); 
%     xlabel('Photoelectron Energy'); 
%     hold off; 
    outnumpts = 100; 
    outdx = (xin(end)-xin(1))/outnumpts; 
    xout = xin(1):outdx:xin(end); 
    yout = Spectrum(xout, gaussian, paramout); 
    yout_abs = abs(yout); 
    yout_phi = mod(angle(yout), 2*pi); 
%     yout_abs = yout(1,:); 
%     yout_phi = mod(yout(2,:),2*pi); 

    
    % text and color settings
    text_size = 12; 
    line_weight = 1.5; 
    abs_color = [0 0 0.8]; 
    phi_color = [0 0.7 1]; 
    lgnd_pad = 1.4; 
    lgnd_pos = 'best';
 
    if plotting == 1
        fh = figure('Position', [10 600 560 400]); 
        hold on; 
        % amplitude axis
        ax1 = gca; 
        ax1_pos = ax1.Position; 
        ax1.XColor = 'k'; 
        ax1.XLabel.String = 'Photoelectron Energy (eV)'; 
        ax1.YLabel.String = 'Amplitude (arbitrary units)';
        ax1.FontSize = text_size; 
        ax1.YColor = abs_color; 
        ax1.LineWidth = line_weight*0.5; 
        hold(ax1, 'on'); 
        s1 = scatter(ax1, xin, yin_abs, 'o'); 
        s1.MarkerEdgeColor = abs_color; 
        s1.LineWidth = line_weight; 
        s1.DisplayName = 'amplitude data'; 
        l1 = line(xout, yout_abs, 'Parent', ax1, ...
            'Color', abs_color, 'LineStyle', '-', 'LineWidth', line_weight, ...
            'DisplayName', 'total amplitude fit'); 
        goodplot(24)
        for ii=1:size(gaussian,1) 
            line(xout, abs(Spectrum(xout, gaussian(ii,:), paramout(ii,:))), 'Parent', ax1, ...
                'Color', abs_color, 'LineStyle', '--', 'LineWidth', line_weight, ...
                'HandleVisibility', 'off'); 
        end

        % phase axis
        ax2 = axes('Position', ax1_pos, ...
            'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
            'Color', 'none', 'LineWidth', line_weight*0.5); 
        ax2.XColor = 'none'; 
        ax2.YColor = phi_color; 
        ax2.YLabel.String = 'Phase (radians)'; 
        ax2.FontSize = text_size; 
        hold(ax2, 'on'); 
        s2 = scatter(ax2, xin, yin_phi, '+'); 
        s2.MarkerEdgeColor = phi_color; 
        s2.LineWidth = line_weight; 
        s2.DisplayName = 'phase data'; 
        l2 = line(xout, unwrap(yout_phi), 'Parent', ax2, ...
            'Color', phi_color, 'LineStyle', '-', 'LineWidth', line_weight, ...
            'DisplayName', 'total phase fit'); 
        goodplot(24)
    %     for i=1:1:(length(paramout(:,1)))
    %         if length(fix) > 1
    %             tmp = mydist_fixpeaks(xout, fix(i), paramout(i,:), slope); 
    %         else
    %             if peakflag==1
    %                 tmp = mydist_fixpeaks(xout, fix(i), paramout(i,:), slope); 
    %             else
    %                 tmp = mydist_fixwidth(xout, fix, paramout(i,:), slope); 
    %             end
    %         end
    % %         tmp_phi = angle(tmp(:,1).*exp(1j*tmp(:,2))); 
    %         tmp_phi = mod(tmp(:,2), 2*pi); 
    % %         tmp_phi = tmp(:,2); 
    %         line(xout, tmp_phi, 'Parent', ax2, ...
    %             'Color', phi_color, 'LineStyle', '--', 'LineWidth', line_weight, ...
    %             'HandleVisibility', 'off'); 
    %     end

        linkaxes([ax1,ax2],'x')
    %     xlabel('photoelectron energy (eV)'); 


        % harmonic labels
        axl = AddHarmonicAxis(ax2, IP, IP_label, wavelength, 1);
        for i = 1:numel(IP)
%             axl(i).XLabel.Position = [ -1.9    0.99];
            axl(i).FontSize = text_size*0.75; 
            axl(i).YAxisLocation = 'right'; 
%             axl(i).XLabel.String = IP_label(i); 
        end
        ax1.XLim = [xout(1), xout(end)]; 
%         ax1.YLim = [0, max(yout_abs)*lgnd_pad]; 
        lgnd = legend([s1, l1, s2, l2], 'Location', lgnd_pos); 
        lgnd.FontSize = text_size*0.6; 

        hold off; 
    end

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

%     function Yout = Spectrum(E, gaussian, p)
%         Yout = 0; 
%         Gauss = @(x,A,mu,sig,alpha) A * alpha/2 * exp(alpha/2 * (2*mu - alpha*sig^2 - 2*x)) ...
%             .* (1 - erf( (mu + alpha*sig^2 -x)/sqrt(2)/sig )); 
%         if size(p,2) == 2
%             Phase = @(x,b,mu) exp(1j .* b); 
%             % sum the 2w signal
%             for n = 1:size(p,1)
%                 Amp = gaussian(n,1); 
%     %             Amp = 1; 
%                 E0 = gaussian(n,2); 
%                 wid = gaussian(n,3); 
%                 alpha = gaussian(n,4); 
% 
%                 b = mod(p(n,1),2*pi); 
% 
%                 Yout = Yout + Gauss(2*E0-E,Amp,E0,wid,alpha).*Phase(E,b,E0);
%             end
%         elseif size(p,2) == 3
%             Phase = @(x,b,c,mu) exp(1j .* (b + c.*(x-mu)) ); 
%             % sum the 2w signal
%             for n = 1:size(p,1)
%                 Amp = gaussian(n,1); 
%     %             Amp = 1; 
%                 E0 = gaussian(n,2); 
%                 wid = gaussian(n,3); 
%                 alpha = gaussian(n,4); 
% 
%                 b = mod(p(n,1),2*pi);
%                 c = p(n,2); 
% 
%                 Yout = Yout + Gauss(2*E0-E,Amp,E0,wid,alpha).*Phase(E,b,c,E0);
%             end
%         else
%             error('invalid guess input'); 
%         end
%         
%         yout_mat = [abs(Yout); angle(Yout)]; 
% 
%     end


end