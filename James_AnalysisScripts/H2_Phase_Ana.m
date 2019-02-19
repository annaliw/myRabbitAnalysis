%% Load Data
load('/Users/annaliw/code/myRabbitAnalysis/HydrogenSpectra.mat')
% Change this to point to location of Hydrogen Data
X = E;
Y2 = abs(restackedSpectra);
Y3 = angle(restackedSpectra);

load('/Users/annaliw/code/myRabbitAnalysis/ArgonSpectra.mat')
%Change this to point to location of Ar Ref Data
XA = E;
YA = restackedSpectra;

% addpath('/Users/jcryan/Dropbox/AttoLab/RABBITT/')
% %This should point to directory where fitting functions live.

%% Fit Ar spectrum
E0 = [2.61 5.69, 8.73 11.75];
p0 = [-1.1 -0.45 0.40 1.5]; 

for n = 1:numel(E0)
    
%    n=2;
    ind = find(XA > (E0(n)-0.7) & XA < (E0(n)+0.7) );
    
    [A,fval] = Fit1ComplexGaussian(XA(ind),YA(ind), [1,p0(n),0.7,0.2,E0(n)]);
    
    ArSB(n) = A(2);
    ArSSB(n) = A(3);
end

figure; 
plotyy(12:2:18, ArSB, 12:2:18, ArSSB)
%% 12th SB
ind = find(X>1.5 & X <2.856);

Y_in = Y2(ind).* exp(1j.*Y3(ind));

config.E0 = [1.767, 1.99, 2.211, 2.465, 2.722];
config.width = 0.07;
config.guess = [1,-1.5,1];

[A,fval] = FitnessWholePizza(X(ind),Y_in, config); 
fval
p = reshape(A,5,4);

figure (12);
subplot(2,1,1)
plot(5:-1:1,p(:,2)','o')
hold on
plot(3,ArSB(1),'o')
subplot(2,1,2)
plot(5:-1:1,p(:,3)');
hold on
plot(3,ArSSB(1),'o')

SB12 = p(:,2)';
A12 = A;
%% 14th SB
ind = find(X>4.7 & X <5.97);

Y_in = Y2(ind).* exp(1j.*Y3(ind));

config.E0 = [4.85, 5.078, 5.31, 5.55, 5.78];
config.width = 0.08;
config.guess = [1,-0.8,0.7];

[A,fval] = FitnessWholePizza(X(ind),Y_in, config); 
fval
p = reshape(A,5,4);

figure (14);
subplot(2,1,1)
plot(5:-1:1,p(:,2)','o')
hold on
plot(3,ArSB(2),'o')
subplot(2,1,2)
plot(5:-1:1,p(:,3)');
hold on
plot(3,ArSSB(2),'o')

SB14 = p(:,2)';
A14 = A;
%% 16th SB
%Not working so great.... This is where I stopped.
ind = find(X>7.72 & X<8.92);

Y_in = Y2(ind).* exp(1j.*Y3(ind));

config.E0 = [7.9300    8.1580    8.3900    8.6300    8.8600];
config.width = 0.15;
config.guess = [1,-0.8,0.1];

[A,fval] = FitnessWholePizza(X(ind),Y_in, config); 
fval
p = reshape(A,5,3);

figure (16);
subplot(2,1,1)
plot(5:-1:1,p(:,2)','o')
hold on
plot(3,ArSB(3),'o')
subplot(2,1,2)
plot(5:-1:1,p(:,3)');
hold on
plot(3,ArSSB(3),'o')

SB16 = p(:,2)';

%%
ind = find(X>10.72 & X<12.23);

Y_in = Y2(ind).* exp(1j.*Y3(ind));

config.E0 = [7.9300    8.1580    8.3900    8.6300    8.8600] + 3.08;
config.width = 0.23;
config.guess = [1,-0.8,0.1];

[A,fval] = FitnessWholePizza(X(ind),Y_in, config); 
fval
p = reshape(A,5,3);

figure (16);
plotyy(5:-1:1,p(:,2)',5:-1:1,p(:,3)')

SB18 = p(:,2)';

%%
figure (899)
plot(5:-1:1,SB12 - ArSB(1))
hold on
plot(5:-1:1,SB14 - ArSB(2))
plot(5:-1:1,SB16 - ArSB(3))
plot(5:-1:1,SB18 - ArSB(4))
hold off

%% Compare 12th and 14th SB phase

figure(20)
subplot(2,1,1)
ind = find(X>1.5 & X <2.856);
plot(X(ind)-mean(X(ind)),Y2(ind))
subplot(2,1,2)
plot(X(ind)-mean(X(ind)),Y3(ind)-mean(Y3(ind)))

ind = find(X>4.7 & X <5.97);
subplot(2,1,1)
hold on
plot(X(ind)-mean(X(ind))+0.066,Y2(ind))
subplot(2,1,2)
hold on
plot(X(ind)-mean(X(ind))+0.066,Y3(ind)-mean(Y3(ind)))

ind = find(X>7.72 & X<8.92);
subplot(2,1,1)
plot(X(ind)-mean(X(ind))+0.049,Y2(ind))
hold off
subplot(2,1,2)
plot(X(ind)-mean(X(ind))+0.049,Y3(ind)-mean(Y3(ind)))
hold off
