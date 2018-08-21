function [Counts, Energy]=Convert_Eng_V2(tof,tcounts,param,Espec)
%[Counts, Energy]=Convert_Eng(tof,tcounts,[t_prompt,A],[Emin,Emax,nbins])
%tof and t_prompt must be in same units.  
%A is the conversion from time units to energy units 
%E = A(1) + A(2)/(tof-t_prompt)^2 + A(3)/(tof-t_prompt)^4 + ...
%tcounts is the histogram of time bins. 
%if A is not given, a previous value is used. 

%Get Prompt time and ToF Conversion
if numel(param)>2
    t_prompt = param(1);
    A = param(2:end);
elseif numel(param) == 2
    t_prompt = param(1);
    A = [0,param(2)];
elseif numel(param) == 1
    t_prompt = param(1);
    % A = [0,6.9386E5];
    A = [0.7980, 4.6454e+06, -2.5990e+10];
else
    t_prompt = 0;
    A = [0.7980, 4.6454e+06, -2.5990e+10];
end

%Convert time to energy
t0 = t_prompt;%-1.993E-9; %t0=prompt time + light travel time to detector
t = (tof-t0); %Adjust for prompt time and convert to ns
x = 1./t.^2;

Eng = 0;
for i = 0:numel(A)-1
    Eng = Eng + A(i+1) .* x.^i;
end
Eng(tof<=t_prompt)=0;

%Make Linear Energy axis
Emin=Espec(1); Emax=Espec(2); nbins=Espec(3);
Estep = (Emax-Emin) / nbins;
Energy = Emin:Estep:(Emax-Estep);

%Overlap Matrix
OM = zeros( numel(Energy), numel(Eng) );

for ind1 = 1:numel(Energy)-1
    for ind2 = 1:numel(Eng)-1
         
        dE2 = Eng(ind2) - Eng(ind2+1);
        dE = min([Energy(ind1+1), Eng(ind2)]) - max([Energy(ind1), Eng(ind2+1)]);
        
        if dE > 0
            OM( ind1, ind2 ) = dE ./ dE2 ;
        end               
    end
end

%Convert ToF to Energy.
if (numel(tof) == size(tcounts,1))
    Counts = zeros( size(tcounts,2), numel(Energy) );
    for ind = 1:size(tcounts,2)
        %Counts(ind, :) = OM * [tcounts(2:end,ind);0];
        Counts(ind, :) = OM * circshift(tcounts(:,ind),-1);
    end
elseif (numel(tof) == size(tcounts,2))
    Counts = zeros( size(tcounts,1), numel(Energy) );
    for ind = 1:size(tcounts,1)
        %Counts(ind, :) = OM * [tcounts(ind,2:end)';0];
        Counts(ind, :) = OM * circshift(tcounts(ind,:)',-1);
    end
else
    Counts = 0;
end

Energy = Energy + Estep/2;

end

