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