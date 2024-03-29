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

%% make distribution with slopes
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

%% loop over different kinds of peaks w/ slopes

% make list of peak spacings
xlist = 0.05:0.025:1; 

% heights
delta_amp = 0.1:0.1:1; 
% widths
widths = 0.02:0.02:0.2; 
% phases
delta_phase = 0.05:0.1:3.05; 
% slopes 
% slopes = 0:0.1:2; 
slopes=0; 

% store bvals
bval_array = zeros([numel(delta_phase) numel(slopes)]); 
trace_array = zeros([size(xlist) numel(delta_phase) numel(slopes)]); 

for nn=1:numel(delta_phase)
    for mm=1:numel(slopes)
        delta_x = xlist*widths(5)*5; 
        p1 = [0, slopes(mm)]; 
        p2 = [delta_phase(nn), slopes(mm)]; 
        phase = [p1; p2]; 

        Yout_array = zeros(numel(delta_x), numel(E)); 
        param_array = zeros(2, 5, numel(delta_x)); 

        for ii=1:numel(delta_x)
            IP = [0, delta_x]; 
            g1 = [1, 0, widths(5)]; 
            g2 = [delta_amp(end), delta_x(ii), widths(5)]; 
            gaussian = [g1; g2]; 

            Yout = Spectrum(E, gaussian, phase).'; 
            Yout_array(ii,:) = Yout; 

            % run fit
            Yout_data = Yout + 0.1*Yout.*rand(size(Yout)); 
            [paramout, fval] = complexfit_section_bootstrap(0, E, Yout_data, gaussian, [0.6, 1; 0.6, 1], 0); 
            param_array(:,:,ii) = paramout; 
        end

        x = (delta_x/widths(5))'; 
        y = squeeze(param_array(2,4,:));
        
        trace_array(:,:,nn,mm) = y; 

        [tmp1, tmp2] = min(abs(y-0.99*y(end))); 
        bval = x(tmp2); 
        bval_array(nn,mm) = bval; 
    end
    
%     figure; hold on; 
%     yyaxis left; 
%     plot(delta_x, squeeze(param_array(1,4,:)));
%     plot(delta_x, squeeze(param_array(2,4,:)));
%     ylabel('phase'); 
%     yyaxis right; 
%     plot(delta_x, squeeze(param_array(1,5,:)));
%     plot(delta_x, squeeze(param_array(2,5,:)));
%     ylabel('phase slope'); 
%     xlabel('peak separation'); 
%     hold off; 

end
check_if_done = 'done!'

% figure; hold on; 
% yyaxis left; 
% plot(delta_x, squeeze(param_array(1,4,:)));
% plot(delta_x, squeeze(param_array(2,4,:)));
% ylabel('phase'); 
% yyaxis right; 
% plot(delta_x, squeeze(param_array(1,5,:)));
% plot(delta_x, squeeze(param_array(2,5,:)));
% ylabel('phase slope'); 
% xlabel('peak separation'); 
% hold off; 

%%
x = delta_x'; 
y = squeeze(param_array(2,4,:));

f = fit(x,y,'exp2'); 
fparam = coeffvalues(f); 
bval = -max(abs(fparam(2)), abs(fparam(4))); 

%%
for nn=1:numel(delta_phase)
    for mm=1:numel(slopes)
        x = (delta_x/peak_width)'; 
        y = trace_array(:,:,nn,mm);

        [tmp1, tmp2] = min(abs(y-0.99*y(end))); 
        bval = x(tmp2); 
        bval_array(nn,mm) = bval; 
    end
end

%% functions
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