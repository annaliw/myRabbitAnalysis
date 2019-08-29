% % make list of peak spacings
% delta_x_spacing = 0.05:0.025:1; 
% 
% % heights
% delta_amp = 0.1:0.05:1; 
% % widths
% widths = 0.01*delta_amp; 
% % phases
% delta_phase = 0.05:0.1:3.05; 
% % slopes 
% slopes = 0:0.1:2; 
% % slopes=0; 

delta_x_spacing = 0.05:0.5:1; 
delta_amp = 0.1:0.5:1; 
widths = 0.1*delta_amp; 
delta_phase = 0.05:0.5:3.05; 
slopes = 0:0.1:2; 

peak_width = 0.01; 

% store bvals
bval_array = zeros([numel(delta_phase) numel(widths) numel(delta_amp)]); 
trace_array = zeros([size(delta_x_spacing) numel(delta_phase) numel(widths) numel(delta_amp)]); 

for ll=1:numel(delta_phase)
for nn=1:numel(widths) 
    for mm=1:numel(delta_amp)
        delta_x = (delta_x_spacing)*peak_width*3;
        xlist = (-3:0.05:6)*peak_width; 
        p1 = [0, slopes(10)]; 
        p2 = [delta_phase(ll), slopes(10)];
	% p1 = 0; 
	% p2 = delta_phase(nn);  
        phase = [p1; p2]; 

        Yout_array = zeros(numel(delta_x), numel(xlist)); 
        param_array = zeros(2, 5, numel(delta_x)); 

        for ii=1:numel(delta_x)
            IP = [0, delta_x]; 
            g1 = [1, 0, peak_width]; 
            g2 = [delta_amp(mm), delta_x(ii), widths(nn)]; 
            gaussian = [g1; g2]; 

            Yout = Spectrum(xlist, gaussian, phase).'; 
            Yout_array(ii,:) = Yout; 

            % run fit
            Yout_data = Yout + 0.1*(real(Yout).*rand(size(Yout)) + 1j*imag(Yout).*rand(size(Yout))); 
            [paramout, fval] = complexfit_section_bootstrap(0, xlist, Yout_data, gaussian, [0.6, 1; 0.6, 1], 0); 
            param_array(:,:,ii) = paramout; 
        end

        x = (delta_x/peak_width)'; 
        y = squeeze(param_array(2,4,:));
        
        trace_array(:,:,ll,nn,mm) = y; 

        [tmp1, tmp2] = min(abs(y-0.99*y(end))); 
        bval = x(tmp2); 
        bval_array(ll,nn,mm) = bval; 
    end
    

end
end

% save data
% save('/farmshare/user_data/annaliw/peak_sim_deltaAmpWidthPhase.mat', 'delta_x', 'delta_phase', 'widths', 'delta_amp', 'trace_array', 'bval_array')
save('test.mat', 'delta_x', 'delta_phase', 'widths', 'delta_amp', 'trace_array', 'bval_array'); 
