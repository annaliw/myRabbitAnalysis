function [paramout, fval] = fit2OmegaSum(xin, yin, gaussian, guess)

    chi2 = @(x) sum( (yin - Spectrum(xin, gaussian, x) ).^2 ); 
    [paramout, fval] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));


    figure; hold on; 
    scatter(xin, yin); 
    plot(xin, Spectrum(xin, gaussian, paramout));
    plot(xin, Spectrum(xin, gaussian, guess)); 
    hold off; 


    function Yout = Spectrum(E, gaussian, p)
        Yout = 0; 
        Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        Phase = @(x,a,b) exp( a + 1j .* (b ) ); 
        % sum the 2w signal
        for n = 1:size(p,1)
            Amp = gaussian(n,1); 
            E0 = gaussian(n,2); 
            wid = gaussian(n,3); 
            
            a = p(n,1); 
            b = p(n,2); 

            Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,a,b);
        end

    end


end