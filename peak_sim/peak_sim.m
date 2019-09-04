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

peak_spacing = 0.05:0.1:3.05; 
amps = 0.1:0.5:1; 
peak_width = 0.01; 
widths = peak_width*amps; 
phases = 0.05:0.5:3.05; 
slopes = 0:0.1:2; 


% decide whether to run this variable
phase_flag = 1; 
slope_flag = 1; 
amp_flag = 0; 
width_flag = 0; 

if phase_flag == 1 
    phase_loop = phases; 
else
    phase_loop = 0; 
end
    
if slope_flag == 1
    slope_loop = slopes; 
else
    slope_loop = 1; 
end

if amp_flag == 1
    amp_loop = amps; 
else
    amp_loop = 1; 
end

if width_flag == 1
    width_loop = width; 
else
    width_loop = 1*peak_width; 
end

delta_x = (peak_spacing)*peak_width*3;
xlist = (-3:0.05:6)*max(delta_x); 
% store bvals
% bval_array = zeros([numel(phase_loop) numel(slope_loop) numel(amp_loop) numel(width_loop)]); 
phaseError_array = zeros([numel(peak_spacing) numel(phase_loop) numel(slope_loop) numel(amp_loop) numel(width_loop)]); 
sim_spectrum_array = zeros([numel(xlist) numel(peak_spacing) numel(phase_loop) numel(slope_loop) numel(amp_loop) numel(width_loop)]); 
fit_spectrum_array = sim_spectrum_array; 

for phase_ind=1:numel(phase_loop)
for slope_ind=1:numel(slope_loop)
for amp_ind=1:numel(amp_loop)
for width_ind=1:numel(width_loop)
    p1 = [0, slope_loop(slope_ind)]; 
    p2 = [phase_loop(phase_ind), slope_loop(slope_ind)];
    phase = [p1; p2]; 

%     Yout_array = zeros(numel(delta_x), numel(xlist)); 
    param_array = zeros(2, 5, numel(delta_x)); 

    for ii=1:numel(delta_x)
        IP = [0, delta_x]; 
        g1 = [1, 0, peak_width]; 
        g2 = [amp_loop(amp_ind), delta_x(ii), width_loop(width_ind)]; 
        gaussian = [g1; g2]; 

        Yout = Spectrum(xlist, gaussian, phase).'; 
%         Yout_array(ii,:) = Yout; 

        % run fit
        Yout_data = Yout + 0.1*(real(Yout).*rand(size(Yout)) + 1j*imag(Yout).*rand(size(Yout))); 
        sim_spectrum_array(:,ii,phase_ind,slope_ind,amp_ind,width_ind) = Yout_data; % save simulated spectrum
        [paramout, fval] = complexfit_section_bootstrap(0, xlist, Yout_data, gaussian, [0.6, 1; 0.6, 1], 0); 
        param_array(:,:,ii) = paramout; 
        fit_spectrum_array(:,ii,phase_ind,slope_ind,amp_ind,width_ind) = Spectrum(xlist, paramout(:,1:3), paramout(:,4:5)).'; 
    end

    x = (delta_x/peak_width)'; 
    y = squeeze(abs(phase_loop(phase_ind)-param_array(2,4,:))./phase_loop(phase_ind))'; 

    phaseError_array(:,phase_ind,slope_ind,amp_ind,width_ind) = y; 

%     [tmp1, tmp2] = min(abs(y-0.99*y(end))); 
%     bval = x(tmp2); 
%     bval_array(phase_ind,slope_ind,amp_ind,width_ind) = bval; 
end
end
end
end

% save data
% save('/farmshare/user_data/annaliw/peak_sim_deltaAmpWidthPhase.mat', 'delta_x', 'delta_phase', 'widths', 'delta_amp', 'trace_array', 'bval_array')
save('test.mat', 'delta_x', 'phase_loop', 'slope_loop', 'amp_loop', 'width_loop', 'phaseError_array'); 
