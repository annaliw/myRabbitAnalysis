% run CO2_cross_sections.m first section for H2 XUV only data and save
% output as 'H11_param', etc. 

%% calculate branching ratios from data
numstates = 4; 
vstate_params = cat(3, H11_param2, H13_param, H15_param); 
volume_list = zeros([size(vstate_params,3), numstates]); 
br_list = volume_list; 
for jj=1:size(vstate_params,3)
%     for ii=1:size(vstate_params,1)
    for ii=1:numstates
        volume_list(jj,ii) = vstate_params(ii,1,jj)*vstate_params(ii,3,jj); 
    end
    br_list(jj,:) = volume_list(jj,:) ./ sum(volume_list(jj,:)); 
end

%% plot branching ratios

figure; hold on; 
for jj=1:size(vstate_params,3)
    plot((1:numstates)-1, br_list(jj,:), 'o-'); 
end
legend('H11', 'H13', 'H15'); 
xlabel('v-state'); 
ylabel('branching ratio'); 
goodplot(); 

%% deconvolve

range = min([numel(H11_xdata), numel(H13_xdata), numel(H15_xdata)]); 
xdata_array = cat(1, H11_xdata(1:range)+11*1240/810, ...
                     H13_xdata(1:range)+13*1240/810, ...
                     H15_xdata(1:range)+15*1240/810); 
ycalc_array = cat(1, H11_ycalc(1:range), H13_ycalc(1:range), H15_ycalc(1:range)); 

ygate_array = zeros(size(xdata_array)); 
yreal_array = zeros(size(xdata_array)); 
yrmdr_array = zeros(size(xdata_array)); 
for ii=1:size(xdata_array,1)
    ygate_array(ii,:) = Gauss(xdata_array(ii,:), ...
        APT_param(ii,1), APT_param(ii,2), APT_param(ii,3)); 
    [yreal_array(ii,:), yrmdr_array(ii,:)] = deconv(ycalc_array(ii,:), ygate_array(ii,:)); 
end

% figure; hold on; 
% plot(xdata_array(1,:), ycalc_array(1,:), 'DisplayName', 'data'); 
% plot(xdata_array(1,:), yrmdr_array(1,:), 'DisplayName', 'data'); 
% plot(xdata_array(2,:), ycalc_array(2,:), 'DisplayName', 'data'); 
% plot(xdata_array(2,:), yrmdr_array(2,:), 'DisplayName', 'deconv'); 
% plot(xdata_array(3,:), ycalc_array(3,:), 'DisplayName', 'deconv'); 
% plot(xdata_array(3,:), yrmdr_array(3,:), 'DisplayName', 'deconv'); 
% legend; 

% fit deconvolutions
deconv_param = zeros([3 size(H11_param,1) size(H11_param,2)]); 
for ii=1:size(deconv_param,1)
    [paramout, paramout_gauss, fval] = complexfit_section_full(wavelength, ...
        xdata_array(ii,:)-(2*ii+9)*1240/810, abs(yrmdr_array(ii,:))', yrmdr_array(ii,:)', 1, 1); 
    deconv_param(ii,:,:) = paramout_gauss; 
end

% calculate branching ratios from deconvolution
numstates = 6; 
vstate_params = deconv_param; 
volume_list = zeros([size(vstate_params,3), numstates]); 
br_list = volume_list; 
for jj=1:size(vstate_params,3)
%     for ii=1:size(vstate_params,1)
    for ii=1:numstates
        volume_list(jj,ii) = vstate_params(ii,1,jj)*vstate_params(ii,3,jj); 
    end
    br_list(jj,:) = volume_list(jj,:) ./ sum(volume_list(jj,:)); 
end


%%
function yout = Gauss(x,A,mu,sig)
    yout = A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
end




