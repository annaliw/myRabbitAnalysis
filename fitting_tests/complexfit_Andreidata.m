open('plotforAnna.fig'); 
h = gcf; 

%%
h = findobj(gca,'Type','line'); 
x=get(h,'Xdata'); 
y=get(h,'Ydata'); 

x1 = cell2mat(x(1)); 
x2 = cell2mat(x(2)); 
x3 = cell2mat(x(34)); 

y1 = cell2mat(y(1)); 
y2 = cell2mat(y(2)); 
y3 = cell2mat(y(34)); 

figure; hold on; 
plot(x1, y1, 'DisplayName', 'x1'); 
plot(x2, y2, 'DisplayName', 'x2'); 
yyaxis('right')
plot(x3, y3, 'DisplayName', 'x3'); 
legend; 

%% 331-389 triple peak 
stpt = 331; 
edpt = 389; 
test_x = x1(stpt:edpt); 
test_yamp = y1(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 

% four known peaks
guess = [2e-03 x2(360) 0.2 -2 ; 4e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2 ; 1.5e-03 8.1 0.15 -3]; 
ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1]; 
lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1];  
% double center peaks
% guess = [2e-03 x2(357) 0.2 -2 ; 2e-03 x2(364) 0.2 -2.2]; 
% ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1]; 
% lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 -0.1];  
% three peaks
% guess = [2e-03 x2(360) 0.2 -2 ; 4e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2]; 
% lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
% ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
% two peaks
% guess = [4e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2]; 
% lb = guess - guess.*[0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
% ub = guess + guess.*[0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1]; 

%% 210-257 triple peak
stpt = 210; 
edpt = 257; 
test_x = x1(stpt:edpt); 
test_yamp = y1(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 

% three peaks
guess = [3e-03 4.878 0.2 -2.5 ; 1.5e-03 5.18 0.2 2.5 ; 2e-03 5.36 0.2 -2.9]; 
lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];
ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1]; 

%% 400-493 triple peak
stpt = 400; 
edpt = 493; 
test_x = x1(stpt:edpt); 
test_yamp = y1(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 

% two peaks
guess = [2.5e-03 9.144 0.2 -2.9 ; 1.5e-03 9.5 0.2 2 ; 1.5e-03 9.656 0.1 1.5 ; 0.6e-03 9.767 0.2 1.4 ; 1e-03 10.61 0.2 0.7]; 
lb = guess - guess.*[1 0.01 0.1 -2 ; 1 0.01 0.1 2 ; 1 0.01 0.1 2 ; 1 0.01 0.1 2 ; 1 0.01 0.1 2];
ub = guess + guess.*[5 0.01 0.1 -2 ; 5 0.01 0.1 2 ; 5 0.01 0.1 2 ; 5 0.01 0.1 2 ; 5 0.01 0.1 2]; 

%% fmincon
xin = test_x; 
yin = test_ycom; 

% [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess, [], [], [], [], lb, ub); 
[paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess); 

x_out = linspace(xin(1),xin(end),length(xin)*100);
y_out = mydist(x_out, paramout); 

% figure; hold on; 
% yyaxis('left')
% scatter(xin, abs(yin), 'o'); 
% plot(x_out, abs(y_out), '--'); 
% plot(x_out, abs(mydist(x_out, guess)), '-'); 
% yyaxis('right')
% scatter(xin, angle(yin), '+'); 
% plot(x_out, angle(y_out)); 
% plot(x_out, angle(mydist(x_out, guess)), '-'); 
% hold off; 

fh = figure;
tmp = plot(x_out, abs(y_out)); 
IP = [13.778, 17.706, 18.077, 19.394]; 
axl = AddHarmonicAxis(fh,IP,810);

axl(1).XLabel.String = 'X';
axl(2).XLabel.String = 'A';
axl(3).XLabel.String = 'B';
axl(4).XLabel.String = 'C';

hold on; 
yyaxis('left')
scatter(xin, abs(yin), 'o'); 
plot(x_out, abs(y_out), '-'); 
plot(x_out, abs(mydist(x_out, paramout(1,:))), '--'); 
plot(x_out, abs(mydist(x_out, paramout(2,:))), '--'); 
plot(x_out, abs(mydist(x_out, paramout(3,:))), '--'); 
plot(x_out, abs(mydist(x_out, paramout(4,:))), '--');  
% plot(x_out, abs(mydist(x_out, paramout(5,:))), '--');  
yyaxis('right')
scatter(xin, angle(yin), '+'); 
plot(x_out, angle(y_out), '-'); 
plot(x_out, angle(mydist(x_out, paramout(1,:))), '--'); 
plot(x_out, angle(mydist(x_out, paramout(2,:))), '--'); 
plot(x_out, angle(mydist(x_out, paramout(3,:))), '--'); 
plot(x_out, angle(mydist(x_out, paramout(4,:))), '--'); 
% plot(x_out, angle(mydist(x_out, paramout(5,:))), '--'); 

delete(tmp); 
xlim([x_out(1), x_out(end)])
hold off; 


%% initial guess tests
xin = test_x; 
yin = test_ycom; 

start_guess = [3e-03 4.878 0.2 -2.5 ; 1.5e-03 5.18 0.2 2.5 ; 2e-03 5.36 0.2 -2.9]; 

p = 2; 
dist = 0.5; 
nsteps = 20; 
guess(1,p) = start_guess(1,p) - start_guess(1,p)*(dist/2)*nsteps; 
lslist = 1:1:nsteps; 

for i=1:1:nsteps
    guess(1,p) = guess(1,p) + start_guess(1,p)*(dist/2); 
    [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess); 
    lslist(i) = fval; 
end

plist = (1:1:nsteps) .* start_guess(1,p)*(dist/2)*nsteps; 
figure; hold on; 
scatter(plist, lslist); 
hold off; 

%% full spectrum line matching 
fh = figure; 
plot(x2, abs(y2)); 
IP = [13.778, 17.706, 18.077, 19.394]; 
axl = AddHarmonicAxis(fh,IP,810);

axl(1).XLabel.String = 'X';
axl(2).XLabel.String = 'A';
axl(3).XLabel.String = 'B';
axl(4).XLabel.String = 'C';

hold on; 
plot(x1, abs(y1), 'r'); 

hold off; 









