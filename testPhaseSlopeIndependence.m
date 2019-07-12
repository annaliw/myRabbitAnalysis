% make list of peak spacings
delta_x = 0:0.01:0.3; 
% phases
phase1 = 0.5; 
phase2 = 0.55; 
% slopes 
slope1 = 1.2; 
slope2 = 1.2; 

%% make distribution without slopes
p1 = [phase1]; 
p2 = [phase2]; 
phase = [p1; p2]; 
E = -1:0.01:1; 

Yout_array = zeros(numel(delta_x), numel(E)); 
param_array = zeros(2, 4, numel(delta_x)); 

for ii=1:numel(delta_x)
    IP = [0, delta_x]; 
    g1 = [1, 0, 0.07]; 
    g2 = [1, delta_x(ii), 0.07]; 
    gaussian = [g1; g2]; 
    
    Yout = Spectrum(E, gaussian, phase).'; 
    Yout_array(ii,:) = Yout; 
    
    % run fit
    Yout_data = Yout + 0.1*Yout.*rand(size(Yout)); 
    [paramout, fval] = complexfit_section_bootstrap(0, E, Yout_data, gaussian, [0.6; 0.6], 0); 
    param_array(:,:,ii) = paramout; 
end

figure; hold on; 
yyaxis left; 
plot(delta_x, squeeze(param_array(1,4,:)));
plot(delta_x, squeeze(param_array(2,4,:)));
ylabel('phase'); 
xlabel('peak separation'); 
hold off; 

%%

% make distribution with slopes
p1 = [phase1, slope1]; 
p2 = [phase2, slope2]; 
phase = [p1; p2]; 

Yout_array = zeros(numel(delta_x), numel(E)); 
param_array = zeros(2, 5, numel(delta_x)); 

for ii=1:numel(delta_x)
    IP = [0, delta_x]; 
    g1 = [1, 0, 0.07]; 
    g2 = [1, delta_x(ii), 0.07]; 
    gaussian = [g1; g2]; 
    
    Yout = Spectrum(E, gaussian, phase).'; 
    Yout_array(ii,:) = Yout; 
    
    % run fit
    Yout_data = Yout + 0.1*Yout.*rand(size(Yout)); 
    [paramout, fval] = complexfit_section_bootstrap(0, E, Yout_data, gaussian, [0.6, 1; 0.6, 1], 0); 
    param_array(:,:,ii) = paramout; 
end

figure; hold on; 
yyaxis left; 
plot(delta_x, squeeze(param_array(1,4,:)));
plot(delta_x, squeeze(param_array(2,4,:)));
ylabel('phase'); 
yyaxis right; 
plot(delta_x, squeeze(param_array(1,5,:)));
plot(delta_x, squeeze(param_array(2,5,:)));
ylabel('phase slope'); 
xlabel('peak separation'); 
hold off; 





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