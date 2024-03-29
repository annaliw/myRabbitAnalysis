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

if phase_flag == 0
    phases = 0; 
end
    
if slope_flag == 0
    slopes = 1; 
end

if amp_flag == 0
    amps = 1; 
end

if width_flag == 0
    widths = 1*peak_width; 
end

delta_x = (peak_spacing)*peak_width*3;
xlist = (-3:0.05:6)*max(delta_x); 
% store bvals
% bval_array = zeros([numel(phases) numel(slopes) numel(amps) numel(widths)]); 
phaseError_array = zeros([numel(peak_spacing) numel(phases) numel(slopes) numel(amps) numel(widths)]); 
sim_spectrum_array = zeros([numel(xlist) numel(peak_spacing) numel(phases) numel(slopes) numel(amps) numel(widths)]); 
fit_spectrum_array = sim_spectrum_array; 

for phase_ind=1:numel(phases)
for slope_ind=1:numel(slopes)
for amp_ind=1:numel(amps)
for width_ind=1:numel(widths)
    
    phase = [0, slopes(slope_ind); phases(phase_ind), slopes(slope_ind)]; 

%     Yout_array = zeros(numel(delta_x), numel(xlist)); 
    param_array = zeros(2, 5, numel(delta_x)); 

    for ii=1:numel(delta_x)
%         IP = [0, delta_x]; 
        gaussian = [1, 0, peak_width; amps(amp_ind), delta_x(ii), widths(width_ind)]; 

        Yout = Spectrum(xlist, gaussian, phase).'; 

        % run fit
        sim_spectrum_array(:,ii,phase_ind,slope_ind,amp_ind,width_ind) = Yout + 0.1*(real(Yout).*rand(size(Yout)) + 1j*imag(Yout).*rand(size(Yout))); 
        [paramout, fval] = complexfit_section_bootstrap(0, xlist, sim_spectrum_array(:,ii,phase_ind,slope_ind,amp_ind,width_ind), gaussian, [0.6, 1; 0.6, 1], 0); 
%         param_array(:,:,ii) = paramout; 
        fit_spectrum_array(:,ii,phase_ind,slope_ind,amp_ind,width_ind) = Spectrum(xlist, paramout(:,1:3), paramout(:,4:5)).'; 
    end

    phaseError_array(:,phase_ind,slope_ind,amp_ind,width_ind) =  squeeze(abs(phases(phase_ind)-param_array(2,4,:))./phases(phase_ind))'; 

end
end
end
end

% save data
% save('/farmshare/user_data/annaliw/peak_sim_deltaAmpWidthPhase.mat', 'delta_x', 'delta_phase', 'widths', 'delta_amp', 'trace_array', 'bval_array')
save('test.mat', 'delta_x', 'phases', 'slopes', 'amps', 'widths', 'phaseError_array'); 
