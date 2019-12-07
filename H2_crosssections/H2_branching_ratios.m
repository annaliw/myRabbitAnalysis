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

%% import known spectrum
tmp = csvread('/Users/annaliw/code/myRabbitAnalysis/H2_crosssections/H2_PIsigma_1.csv',0,0); 
tmp(:,2) = tmp(:,2)/squeeze(sum(tmp(:,2))); 
H2_PIsigma = tmp; 
tmp = csvread('/Users/annaliw/code/myRabbitAnalysis/H2_crosssections/H2_PIsigma_2.csv',0,0); 
tmp(:,2) = tmp(:,2)/squeeze(sum(tmp(:,2))); 
H2_PIsigma = [H2_PIsigma; tmp]; 
tmp = csvread('/Users/annaliw/code/myRabbitAnalysis/H2_crosssections/H2_PIsigma_3.csv',0,0); 
tmp(:,2) = tmp(:,2)/squeeze(sum(tmp(:,2))); 
H2_PIsigma = [H2_PIsigma; tmp]; 

H2_PIsigma(:,1) = 1240 ./ (H2_PIsigma(:,1)/10); 

%%
function yout = Gauss(x,A,mu,sig)
    yout = A.* exp( -(x-mu).^2 ./ (2.*sig.^2) );
end




