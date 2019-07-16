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
% slopes=0; 

peak_width = 0.01; 

% store bvals
bval_array = zeros([numel(delta_phase) numel(slopes)]); 
trace_array = zeros([size(xlist) numel(delta_phase) numel(slopes)]); 

for nn=1:numel(delta_phase)
    for mm=1:numel(slopes)
        delta_x = xlist*peak_width*3; 
        xlist = (-delta_x*2):0.05:(delta_x*3); 
        p1 = [0, slopes(mm)]; 
        p2 = [delta_phase(nn), slopes(mm)]; 
        phase = [p1; p2]; 

        Yout_array = zeros(numel(delta_x), numel(xlist)); 
        param_array = zeros(2, 5, numel(delta_x)); 

        for ii=1:numel(delta_x)
            IP = [0, delta_x]; 
            g1 = [1, 0, peak_width]; 
            g2 = [delta_amp(end), delta_x(ii), peak_width]; 
            gaussian = [g1; g2]; 
            
            Yout = Spectrum(xlist, gaussian, phase).'; 
            Yout_array(ii,:) = Yout; 

            % run fit
            Yout_data = Yout + 0.1*Yout.*rand(size(Yout)); 
            [paramout, fval] = complexfit_section_bootstrap(0, xlist, Yout_data, gaussian, [0.6, 1; 0.6, 1], 0); 
            param_array(:,:,ii) = paramout; 
        end

        x = (delta_x/peak_width)'; 
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
% check_if_done = 'done!'

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

% save data
save('peak_sim_identical.mat', delta_x, peak_width, delta_phase, slopes, trace_array, bval_array)

