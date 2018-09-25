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
edpt = 368; 
test_x = x2(stpt:edpt); 
test_yamp = y2(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 


%% fmincon
xin = test_x; 
yin = test_ycom; 
guess = [2e-03 x2(360) 0.2 2.08 ; 4e-03 x2(345) 0.3 0.09]; 
ub = guess + guess.*[0.1 0.01 0.1 0.1 ; 0.1 0.01 0.1 0.1]; 
lb = guess - guess.*[0.1 0.01 0.1 0.1 ; 0.1 0.01 0.1 0.1];  

[paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess, [], [], [], [], lb, ub); 

x_out = linspace(xin(1),xin(end),length(xin)*100);
y_out = mydist(x_out, paramout); 

figure; hold on; 
yyaxis('left')
scatter(xin, abs(yin), 'o'); 
plot(x_out, abs(y_out)); 
yyaxis('right')
scatter(xin, angle(yin), '+'); 
plot(x_out, angle(y_out)); 
hold off; 

%% test gaussian log...
gauss_x = -50:1:50; 
gauss_1 = 5*exp(-(gauss_x-20).^2/10^2); 
gauss_2 = 6*exp(-(gauss_x).^2/15^2); 
gauss_3 = 5*exp(-(gauss_x+10).^2/5^2); 

gauss_tot = gauss_1 + gauss_2 + gauss_3; 

figure; hold on; 
plot(gauss_x, gauss_tot); 
yyaxis('right'); 
plot(gauss_x, log(gauss_tot)); 







