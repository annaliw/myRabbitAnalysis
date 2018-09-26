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
% edpt = 368; 
edpt = 389; 
test_x = x2(stpt:edpt); 
test_yamp = y2(stpt:edpt); 
test_ypha = y3(stpt:edpt); 
test_ycom = test_yamp .* exp(-1j*test_ypha); 


%% fmincon
xin = test_x; 
yin = test_ycom; 
guess = [2.5e-03 x2(360) 0.2 -2 ; 2e-03 x2(345) 0.3 0.1 ; 2e-03 x2(378) 0.2 -2.2]; 
ub = guess + guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1]; 
lb = guess - guess.*[0.5 0.01 0.1 -0.1 ; 0.5 0.01 0.1 0.1 ; 0.5 0.01 0.1 -0.1];  
% guess = [1.5e-03 x2(345) 0.3 -0.09]; 
% ub = guess + guess.*[0.2 0.01 0.1 -0.1]/2; 
% lb = guess - guess.*[0.2 0.01 0.1 -0.1]/2;  

[paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess, [], [], [], [], lb, ub); 
% [paramout, fval] = fmincon(@(x) myfit(xin, yin, x), guess); 

x_out = linspace(xin(1),xin(end),length(xin)*100);
y_out = mydist(x_out, paramout); 

figure; hold on; 
yyaxis('left')
scatter(xin, abs(yin), 'o'); 
plot(x_out, abs(y_out), '--'); 
plot(x_out, abs(mydist(x_out, guess)), '-'); 
yyaxis('right')
scatter(xin, angle(yin), '+'); 
plot(x_out, angle(y_out)); 
plot(x_out, angle(mydist(x_out, guess)), '-'); 
hold off; 


%% fmincon parameter robustness
xin = test_x; 
yin = test_ycom; 
% starting_guess = [1e-03 xin(1) 0.2 -1 ; 1e-03 xin(ceil(length(xin)/2)) 0.2 -1]; 
starting_guess = [1e-03 xin(1) 0.2 -1 ; 2e-03 xin((ceil(length(xin)/3))) 0.2 -1 ; 2e-03 xin((2*ceil(length(xin)/3))) 0.2 -1]; 
guess = starting_guess; 
[paramout, fval] = fmincon(@(x) myfit(xin, yin, x), starting_guess); 

x_out = linspace(xin(1),xin(end),length(xin)*100);
y_out = mydist(x_out, paramout); 

figure; hold on; 
yyaxis('left')
scatter(xin, abs(yin), 'o'); 
plot(x_out, abs(y_out), '--'); 
plot(x_out, abs(mydist(x_out, starting_guess)), '-'); 
yyaxis('right')
scatter(xin, angle(yin), '+'); 
plot(x_out, angle(y_out)); 
plot(x_out, angle(mydist(x_out, starting_guess)), '-'); 
hold off; 
%% 

figure; hold on; 
yyaxis('left')
scatter(xin, abs(yin), 'o'); 
plot(x_out, abs(y_out), '-'); 
plot(x_out, abs(mydist(x_out, paramout(1,:))), '--'); 
plot(x_out, abs(mydist(x_out, paramout(2,:))), '--'); 
plot(x_out, abs(mydist(x_out, paramout(3,:))), '--'); 
yyaxis('right')
scatter(xin, angle(yin), '+'); 
plot(x_out, angle(y_out), '-'); 
plot(x_out, angle(mydist(x_out, paramout(1,:))), '--'); 
plot(x_out, angle(mydist(x_out, paramout(2,:))), '--'); 
plot(x_out, angle(mydist(x_out, paramout(3,:))), '--'); 
hold off;





