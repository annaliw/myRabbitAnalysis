function [paramout, fval] = fitGaussianSum(xin, yin, guess, plotting)

%     % FMINSEARCH FIT
%     chi2 = @(x) sum( (yin - Spectrum(xin, x) ).^2 ); 
%     [paramout, fval] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));
%     
%     paramout(4,2) = guess(4,2); 

    % LSQCURVEFIT FIT
    fun = @(x,xdata) Spectrum(xdata, x); 
    [paramout, fval] = lsqcurvefit(fun, guess, xin, yin); 
%     paramout(4,2) = guess(4,2); 

    if plotting == 1
        figure; hold on; 
        scatter(xin, yin); 
        plot(xin, Spectrum(xin, paramout));
        plot(xin, Spectrum(xin, guess)); 
        legend('data', 'fit', 'guess'); 
        hold off; 
    end


    function Yout = Spectrum(E, p)
        Yout = 0; 
        Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        % sum the gaussians
        p = abs(p); 
        for n = 1:size(p,1)
            Amp = p(n,1);
            wid = p(n,3);
            E0 = p(n,2);

            Yout = Yout + Gauss(E,Amp,E0,wid);
        end

    end


end