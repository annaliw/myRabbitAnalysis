function [X,F] = FitnessWholePizza (Ein,Yin,config)

%Check Configs
if (~isfield(config,'E0')) % Peak positions:
    E0in = [1.767, 2.000, 2.211, 2.465, 2.722];
    config.E0 = E0in;
else
    E0in = config.E0;
end

if (~isfield(config,'width')) % Peak width:
    peak_width = 0.065;
    config.width = peak_width;
else
    peak_width = config.width;
end
    
if (~isfield(config,'guess')) % Peak width:
    guess = [1, -1.5, 1]; %[Amp, phase, slope]
else
    guess = config.guess;
end

%Pass Initial Guess
param0 = zeros(numel(E0in),numel(guess)+1);
for N = 1: numel(E0in)
    %param0(N,:) = guess;
    param0(N,:) = [guess, E0in(N)];
end

Yin = Yin ./ max( abs( Yin ) ); %Normalize

RYin = real(Yin);
IYin = imag(Yin);

%Fitness function
chi2 = @(x) sum( (RYin - real( Spectrum(Ein,config,x) ) ).^2 ) ...
          + sum( (IYin - imag( Spectrum(Ein,config,x) ) ).^2 );
    
[X,F] = fminsearch(chi2,param0,optimset('MaxFunEvals', 1E6, 'MaxIter',5E5));

figure(996)
ax = plotyy(Ein,abs(Yin),Ein,angle(Yin));
hold (ax(1), 'on')
plot(ax(1), Ein, abs(Spectrum(Ein,config,X)) );
hold (ax(2), 'on')
plot(ax(2), Ein, angle(Spectrum(Ein,config,X)) );
hold (ax(1), 'off')
hold (ax(2), 'off')
title(sprintf('X^2 = %2.3f', F))

figure(997)
plot(Ein,RYin)
hold on
plot(Ein,IYin)
plot(Ein,real(Spectrum(Ein,config,X)),'--')
plot(Ein,imag(Spectrum(Ein,config,X)),'--')
hold off

%p_out = reshape(X,numel(E0in),4);
%p_out = reshape(X(1:end-2),5,4);
%p_off = X(end-1) + 1j * X(end);

%figure(998)
%plot(p_out(:,4)+angle(p_off),'o')
%plot(p_out(:,2),'o')




%%
    function [Yout] = Spectrum(E,config,p)
        %Vector of x-values
        %E0 vector of centroids
        %p vector of parameters, p has size (E0,4) + 1, the last entry is
        %the offset
        
        Gauss = @(x,A,mu,sig) A.^2 .* exp( -(x-mu).^2 ./ (2.*sig.^2) );
        Phase = @(x,mu,M,b) exp( 1j .* ( M .* (x-mu) + b ) ); 
        
        Yout = zeros(size(E));
        
        E0 = config.E0;
        p = reshape(p, [numel(E0),4]);
        
        %E0 = config.E0;
        for n = 1:numel(E0)
            Amp = p(n,1);
            wid = config.width;
            slope = p(n,3);
            b = p(n,2);
            E0 = p(n,4);
            %E0 = E0(n);
            Yout = Yout + Gauss(E,Amp,E0,wid) .* Phase(E,E0,slope,b);
            
        end
        
        %Yout = Yout + offset; 
                
        
    end

end

