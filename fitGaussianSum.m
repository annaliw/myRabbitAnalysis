function [paramout, fval] = fitGaussianSum(xin, yin, guess)

    chi2 = @(x) sum( (yin - Spectrum(xin, x) ).^2 ); 
    [paramout, fval] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));


    figure; hold on; 
    scatter(xin, yin); 
    plot(xin, Spectrum(xin, paramout));
    plot(xin, Spectrum(xin, guess)); 
    hold off; 


    function Yout = Spectrum(E, p)
        Yout = 0; 
        Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        % sum the gaussians
        for n = 1:size(p,1)
            Amp = p(n,1);
            wid = p(n,3);
            E0 = p(n,2);

            Yout = Yout + Gauss(E,Amp,E0,wid);
        end

    end


end