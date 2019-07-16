function [paramout, fval] = fitGaussianAmplitudes(xin, yin, gaussian, guess, plotting)

%     % FMINSEARCH FIT
%     chi2 = @(x) sum( (yin - Spectrum(xin, x) ).^2 ); 
%     [paramout, fval] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));
%     
%     paramout(4,2) = guess(4,2); 

    % LSQCURVEFIT FIT
    fun = @(x,xdata) Spectrum(xdata, gaussian, x); 
    [paramout, fval] = lsqcurvefit(fun, guess, xin, yin); 
%     paramout(4,2) = guess(4,2); 

    if plotting == 1
        figure; hold on; 
        scatter(xin, yin); 
        plot(xin, Spectrum(xin, gaussian, paramout));
%         plot(xin, Spectrum(xin, gaussian, guess)); 
        legend('data', 'fit', 'guess'); 
        hold off; 
    end


    function Yout = Spectrum(E, gaussian, p)
        Yout = 0; 
        Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        % sum the gaussians
        p = abs(p); 
        for n = 1:length(p)
            Amp = p(n);
            wid = gaussian(n,3);
            E0 = gaussian(n,2);

            Yout = Yout + Gauss(E,Amp,E0,wid);
        end

    end


end