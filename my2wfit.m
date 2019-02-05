function [paramout,fval,x_out,y_out]=my2wfit(xin,yin,peaks,Guess)
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
    
options = optimset('MaxFunEvals',100000, 'MaxIter', 50001, 'TolFun',1E-6, 'TolX',1E-6);
% Guess=[Amp,x0,s,2*s,min(yin)];
% A=mean(abs(yin));
% yin=yin./A;
% ein=ein./A;
% Guess(end)=Guess(end)./A;

[paramout, fval] = fminsearch(@(x) myfit(xin,yin,peaks,x),Guess,options); 

x_out=linspace(xin(1),xin(end),length(xin)*100);

A1 = paramout(1); 
A2 = paramout(2);  
s1 = paramout(3);
s2 = paramout(4); 
% a1 = paramout(5); 
% a2 = paramout(6); 
a1 = 0; 
a2 = 0; 
b1 = paramout(7); 
b2 = paramout(8); 
o  = paramout(9); 

x1 = peaks(1); 
x2 = peaks(2); 

y1 = abs(A1) * exp(1i*(a1*(x_out-x1)+b1)) .* exp(-(x_out-x1).^2/(2*s1)); 
y2 = abs(A2) * exp(1i*(a2*(x_out-x2)+b2)) .* exp(-(x_out-x2).^2/(2*s2)); 

y_out = y1 + y2 + o; 

if nargin<1
    hold on
    plot(x_out,y_out,'k')
end

end


function chi2=myfit(xin,yin,peaks,para)
% parameters
A1 = para(1); 
A2 = para(2); 
s1 = para(3);
s2 = para(4); 
% a1 = para(5); 
% a2 = para(6); 
a1 = 0; 
a2 = 0; 
b1 = para(7); 
b2 = para(8); 
o  = para(9); 

x1 = peaks(1); 
x2 = peaks(2); 

y1 = abs(A1) * exp(1i*(a1*(xin-x1)+b1)) .* exp(-(xin-x1).^2/(2*s1)); 
y2 = abs(A2) * exp(1i*(a2*(xin-x2)+b2)) .* exp(-(xin-x2).^2/(2*s2)); 

chi2_r = (real(yin) - (real(y1+y2) + o)).^2/real(yin).^2; 
chi2_i = (imag(yin) - (imag(y1+y2))).^2/imag(yin).^2; 

chi2 = chi2_r + chi2_i; 


end