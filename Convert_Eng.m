function [Counts, Energy]=Convert_Eng(tof,tcounts,param,Espec)
%[Counts, Energy]=Convert_Eng(tof,tcounts,[t_prompt,A],[Emin,Emax,nbins])
%tof and t_prompt must be in same units.  A is the conversion from time
%units to energy units E=A/(tof-t_prompt)^2
%tcounts is the histogram of time bins. 
%if A is not given, a previous value is used. 

%Get Prompt time and ToF Conversion
if length(param)>1
    t_prompt=param(1);
    A=param(2);
else
    t_prompt=param(1);
%     A=6.9386E5;
    A = 9016000; 
    %A=( 9.10939E-31 * (0.728)^2 ) / (2 * (1E-9)^2 * 1.602177E-19);
end

%Get Ebin info
Emin=Espec(1); Emax=Espec(2); nbins=Espec(3);

%Convert time to energy
t0=t_prompt;%-1.993E-9;%t0=prompt time + light travel time to detector
t=(tof-t0); %Adjust for prompt time and convert to ns

Eng=A./(t.^2);
Eng(tof<=t_prompt)=0;

Counts=zeros(1,length(Eng));
Energy=zeros(1,length(Eng));

    Estep=(Emax-Emin)/nbins; %Energy stepsize
    Counts=zeros(1,nbins);
    Energy=zeros(1,nbins);
    delta=t(2)-t(1);
    
    ind=find(Eng);
    Eng=Eng(ind);
    t=t(ind);
    tcounts=tcounts(ind);
    
    
    for i=1:nbins
    
        E1=Emin+Estep*(i-1); %Low Energy of Bin
        E2=Emin+Estep*i; %High Energy of Bin
        
        if E2<=min(Eng) %if High energy (short time) side of bin outside tof range
            Counts(i)=0;
            
        else %High energy side of bin is inside tof Range
            at1=interp1(Eng,t,E2); %Actual tof corresponding to E2
            index=find(t<at1); %find index closest to E2
            ind_t1=max(index)+1;
            t1=t(ind_t1); %time according to previous bin number
            
            if E1>min(Eng) %if low Energy side of bin is inside tof range
                at2 = interp1(Eng,t,E1); %Actual tof corresponding to E1
                index=find(t<at2); %find index closest to E2
                ind_t2=index(end);
                t2=t(ind_t2);
            
            else %if low Energy side of bin is outside tof range, use edge of tof range
                at2=t(end);
                ind_t2=length(t);
                t2=t(end);
            end

            if (ind_t2>=ind_t1)
                if (ind_t2>ind_t1)
                    for j=ind_t1+1:ind_t2
                        Counts(i)=Counts(i)+tcounts(j); %Add all timebins inside energy bin
                        
                    end
                end
                %split edges of time bin
                Counts(i)=Counts(i)+tcounts(ind_t2+1)*abs(at2-t2)/delta;
                Counts(i)=Counts(i)+tcounts(ind_t1)*abs(at1-t1)/delta;
           
            elseif (ind_t2 < ind_t1)
                Counts(i)=Counts(i)+tcounts(ind_t2+1)*abs(at2-at1)/delta;
            end
        end
        Energy(i)=(E1+E2)/2;
    end
end
