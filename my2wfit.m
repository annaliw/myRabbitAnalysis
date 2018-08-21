function [paramout,fval,x_out,y_out]=my2wfit(xin,yin,Guess)
%Fits data to two overlapping peaks of fft output
%(must select region of data with only two overlapping peaks)
%[paramout,fval,x_out,y_out]=GausExpFit(xin,yin,ein,[A,x0,sigma,tau,offset])
%paramout=[Ag,Ae,x0,sigma,tau,offset]
% 
% if nargin<1
%     xin = get(gco,'XData');
%     yin = get(gco,'Ydata');
%     ein = sqrt(yin);
%     x0g = [sum(xin.*yin)./sum(yin) sum(xin.*yin)./sum(yin)+5];
%     sg  = zeros(1,2) .* sqrt(sum((xin-x0g).^2.*yin))/sum(yin);
%         
%     Guess=[max(yin)./min(yin)/2 max(yin)./min(yin)/2; 
%             x0g; 
%             sg; 
%             2*sg; 
%             min(yin) min(yin)]; % 2x5 array of input parameters
% 
% elseif nargin<4
%     if length(ein)<length(yin)
%         Guess=ein;
%         ein=sqrt(yin);
%         % I don't get this part
%     else
%         x0g=[sum(xin.*yin)./sum(yin) sum(xin.*yin)./sum(yin)+5];
%         sg=sqrt(sum((xin-x0g).^2.*yin))/sum(yin);
%         
%         Guess=[max(yin)./min(yin)/2,max(yin)./min(yin)/2,x0g,sg,2*sg,min(yin)];
%     end
% end
    
%%
options = optimset('MaxFunEvals',100000, 'MaxIter', 50001, 'TolFun',1E-6, 'TolX',1E-6);
% Guess=[Amp,x0,s,2*s,min(yin)];
% A=mean(abs(yin));
% yin=yin./A;
% ein=ein./A;
% Guess(end)=Guess(end)./A;

[paramout, fval] = fminsearch(@(x) myfit(xin,yin,x),Guess,options); 

x_out=linspace(xin(1),xin(end),length(xin)*100);

A1 = paramout(1); 
A2 = paramout(2); 
x1 = paramout(3); 
x2 = paramout(4); 
s1 = paramout(5);
s2 = paramout(6); 
b1 = paramout(7); 
b2 = paramout(8); 
o  = paramout(9); 

y1_r = abs(A1) * cos(b1) .* exp(-(x_out-x1).^2/(2*s1));  
y2_r = abs(A2) * cos(b2) .* exp(-(x_out-x2).^2/(2*s2));  
yf_r = y1_r + y2_r + o; 

y1_i = abs(A1) * sin(b1) .* exp(-(x_out-x1).^2/(2*s1));  
y2_i = abs(A2) * sin(b2) .* exp(-(x_out-x2).^2/(2*s2));  
yf_i = y1_i + y2_i + o; 

y_out = yf_r + 1i*yf_i; 

if nargin<1
    hold on
    plot(x_out,y_out,'k')
end

end


function chi2=myfit(xin,yin,para)
% parameters are an array of the coefficients for each peak
    
%     A1 = para(1);
%     A2 = para(2); 
%     E1 = para(3);
%     E2 = para(4); 
%     s1 = para(5);
%     s2 = para(6); 
%     a1 = para(7);
%     a2 = para(8); 
%     b1 = para(9);
%     b2 = para(10); 
%     o = para(11); 
%     
%     y1 = abs(A1) .* exp((xin-E1).^2/(2.*s1)) .* exp(1i * (a1.*(xin-E1) + b1));
%     y2 = abs(A2) .* exp((xin-E2).^2/(2.*s2)) .* exp(1i * (a2.*(xin-E2) + b2));
%     yf = y1 + y2 + o; 
% 
%     diff = yf - yin; 
%     
%     chi2=sum(real(diff).^2 + imag(diff).^2); 

% simple gaussian test 
A1 = para(1); 
A2 = para(2); 
x1 = para(3); 
x2 = para(4); 
s1 = para(5);
s2 = para(6); 
b1 = para(7); 
b2 = para(8);  
o  = para(9); 

y1_r = abs(A1) * cos(b1) .* exp(-(xin-x1).^2/(2*s1));  
y2_r = abs(A2) * cos(b2) .* exp(-(xin-x2).^2/(2*s2));  
yf_r = y1_r + y2_r + o; 
chi2_r = (real(yin) - yf_r).^2/real(yin).^2; 

y1_i = abs(A1) * sin(b1) .* exp(-(xin-x1).^2/(2*s1));  
y2_i = abs(A2) * sin(b2) .* exp(-(xin-x2).^2/(2*s2));  
yf_i = y1_i + y2_i; 
chi2_i = (imag(yin) - yf_i).^2/imag(yin).^2; 

chi2 = chi2_r + chi2_i; 


end