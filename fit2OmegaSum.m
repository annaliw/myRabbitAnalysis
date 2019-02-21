function [paramout, fval] = fit2OmegaSum(xin, yin, gaussian, guess)

    figure; hold on; 
    yyaxis left
    scatter(xin, abs(yin)); 
    plot(xin, abs(Spectrum(xin, gaussian, guess)));
    yyaxis right
    scatter(xin, angle(yin)); 
    plot(xin, angle(Spectrum(xin, gaussian, guess))); 
    hold off; 
    
    chi2 = @(x) sum((real(yin)-real(Spectrum(xin, gaussian, x))).^2)... 
        + sum((imag(yin)-imag(Spectrum(xin, gaussian, x))).^2); 
    [paramout, fval] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));


    figure; hold on; 
    yyaxis left
    scatter(xin, abs(yin)); 
    plot(xin, abs(Spectrum(xin, gaussian, paramout)));
    yyaxis right
    scatter(xin, angle(yin)); 
    plot(xin, angle(Spectrum(xin, gaussian, paramout))); 
    hold off; 


    function Yout = Spectrum(E, gaussian, p)
        Yout = 0; 
        Gauss = @(x,A,mu,sig) A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        Phase = @(x,a,b) exp( a + 1j .* (b ) ); 
        % sum the 2w signal
        for n = 1:size(p,1)
            Amp = gaussian(n,1); 
%             Amp = 1; 
            E0 = gaussian(n,2); 
            wid = gaussian(n,3); 
            
            a = p(n,1); 
            b = p(n,2); 

            Yout = Yout + Gauss(E,Amp,E0,wid).*Phase(E,a,b);
        end

    end


end