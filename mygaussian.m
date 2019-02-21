function yout = mygaussian(E, p)
    Gauss = @(x,A,mu,sig) A.^2 .* exp( -(x-mu).^2 ./ (2.*sig.^2) );
    % sum the gaussians
    for n = 1:numel(E0)
        Amp = p(n,1);
        wid = p(n,3);
        E0 = p(n,2);

        Yout = Yout + Gauss(E,Amp,E0,wid);
    end

end