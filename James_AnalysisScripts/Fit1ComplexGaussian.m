function [X,F] = Fit1ComplexGaussian (Ein,Yin,guess)
% [X,F] = Fit1ComplexGaussian (Ein,Yin,guess)
% guess = [amp, phase, slope, width, center];

Yin = Yin ./ max( abs( Yin ) ); %Normalize

RYin = real(Yin);
IYin = imag(Yin);

%Fitness function
chi2 = @(x) sum( (RYin - real( Spectrum(Ein,x) ) ).^2 ) ...
          + sum( (IYin - imag( Spectrum(Ein,x) ) ).^2 );
    
[X,F] = fminsearch(chi2,guess,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));

figure(990)
ax = plotyy(Ein,abs(Yin),Ein,angle(Yin));
hold (ax(1), 'on')
plot(ax(1), Ein, abs(Spectrum(Ein,X)) );
hold (ax(2), 'on')
plot(ax(2), Ein, angle(Spectrum(Ein,X)) );
hold (ax(1), 'off')
hold (ax(2), 'off')
title(sprintf('X^2 = %2.3f', F))

figure(991)
plot(Ein,RYin)
hold on
plot(Ein,IYin)
plot(Ein,real(Spectrum(Ein,X)),'--')
plot(Ein,imag(Spectrum(Ein,X)),'--')
hold off


%%
    function [Yout] = Spectrum(E,p)
        
        Gauss = @(x,A,mu,sig) A.^2 .* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        Phase = @(x,mu,M,b) exp( 1j .* ( M .* (x-mu) + b ) ); 
        
        Yout = Gauss(E,p(1),p(5),p(4)) .* Phase(E,p(5),p(3),p(2));
    
    end

end