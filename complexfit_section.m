function [paramout_fixwidth, fval, exitflag, output] = complexfit_section(wavelength, E, data, fitRegion, varargin)
    % ONLY FIXED WIDTH VERSION OF FITTING, fminsearch on lsq
    
    global IP
%     global IP_label
    plotting = 0; 
    fixguess = 0; 
    
    if (~isempty(varargin))
        
        for c=1:length(varargin)
            switch varargin{c}
                case {'FixWidth'}
                    fixpeaks = 0;
                    fixwidth = 1; 
                
                case {'FixPeaks'}
                    fixpeaks = 1; 
                    fixwidth = 0; 

                case {'FitSlope'}
                    slopeflag = 1; 
                    if length(fitRegion) == 2
                        slope_guess = 1; 
                    elseif length(fitRegion) == 3
                        slope_guess = fitRegion(3); 
                    elseif length(fitRegion) == 4
                        slope_guess = fitRegion(4); 
                    else
                        error('Bad fitRegion size.')
                    end
    
%                 case {'FixGuess'}
%                     guess = 
                case {'Plot'}
                    plotting = 1;
                    
                case {'FixGuess'}
                    fixguess = 1; 
                    
            % (continued)
            
            otherwise         
                error(['Invalid optional argument, ', ...
                    varargin{c}]);
            end % switch
        end % for
    else 
        % default to fixed width, no slope, no plotting
        fixpeaks = 0; 
        fixwidth = 1; 
        slopeflag = 0; 
        if length(fitRegion) == 2
            width = 0.1; 
        elseif length(fitRegion) == 3
            width = fitRegion(3); 
        end
        plotting = 0; 
    end % if
    
%     if ~exist('peakflag', 'var')
%         peakflag = 0; 
%     end
%     if ~exist('IP_label', 'var')
%         plotting = 0; 
%     else
%         plotting = 1; 
%     end
%     if length(fitRegion)==2
%         slopeflag = 0; 
%         slope_guess = 0; 
%     elseif length(fitRegion)==3
%         slopeflag = 0;
%         slope_guess = 0; 
%         width = fitRegion(3); 
%     elseif length(fitRegion)==4
%         slopeflag = 1; 
%         slope_guess = fitRegion(4); 
%         width = fitRegion(3); 
%     else
%         error('invalid fitRegion size'); 
%     end
    
    
    % full data set-up
    x = E; 
    y_abs = abs(data).'; 
    y_phi = mod(angle(data), 2*pi).'; 
    n=9:1:19; 
    peaks = zeros(length(IP), length(n)); 
    for i=1:1:length(IP)
        peaks(i,:) = n*(1240/wavelength)-IP(i); 
    end
    peaks = peaks(:).'; 
    
    % fit section set-up
    start = find(abs(x-fitRegion(1))<0.05, 1); 
    stop = find(abs(x-fitRegion(2))<0.05, 1); 
    xin = x(start:stop); 
    yin_abs = y_abs(start:stop); 
    yin_phi = unwrap(y_phi(start:stop)); 
    yin = [yin_abs ; yin_phi].'; 

    % find peaks and their indices
    peaks_guess = peaks; 
    remove = [find(peaks_guess < xin(1)) find(peaks_guess > xin(end))]; 
    peaks_guess(remove) = []; 
    peak_ind = 1:1:length(peaks_guess); 
    for i=1:1:length(peaks_guess)
        peak_ind(i) = find(abs(xin - peaks_guess(i)) < 0.03, 1);  
    end
    amp_guess = yin_abs(peak_ind); % form amplitude guess
    phi_guess = yin_phi(peak_ind); % form phase guess
    slope_guess = zeros(size(peak_ind)) + slope_guess; % form slope guess
    guess_fixwidth = [amp_guess; peaks_guess; phi_guess; slope_guess].'; % full guess matrix
    if slopeflag==0
        guess = guess_fixwidth(:,1:3); 
    else
        guess = guess_fixwidth; 
    end
    if fixguess == 1 % reassign amp, phase, phase slope guesses
        guess_flag = 1; 
        if slopeflag==0
            guess = guess_fixwidth(:,1:3);
            guess(:,1) = fixguess(1); 
%             guess(:,2) = fixguess(2); 
            guess(:,3) = fixguess(3); 
        else
            guess = guess_fixwidth; 
            guess(:,1) = fixguess(1); 
%             guess(:,2) = fixguess(2); 
            guess(:,3) = fixguess(3); 
            guess(:,4) = fixguess(4); 
        end
    end
    
%     if length(fitRegion)==2
%         if fixpeaks==1
%             width = guess(:,2).'; % width will now hold fixed peaks value
%             guess(:,2) = ones(size(guess(:,2)))*0.1; % turn peak center guesses into width guesses
%         else
%             width = fitRegion(3); 
%         end
%     else
%         width = fitRegion(3); 
%     end
    if fixpeaks == 1
        width = guess(:,2).'; % width will now hold fixed peaks value
        guess(:,2) = ones(size(guess(:,2)))*0.1; % turn peak center guesses into width guesses
    else
        width = fitRegion(3); 
    end
    
    
    % do fit
    fun = @(guess, data) complex_lsq(data, width, guess, slope, peakflag); 
    data = [xin; yin_abs.*exp(1j*yin_phi)]; 
    options = optimset('MaxFunEvals', 200000, 'MaxIter', 100000); 
    [paramout_fixwidth, fval, exitflag, output] = fminsearch(@(x) complex_lsq(data, width, x, slopeflag, 0), guess, options); 
    paramout_fixwidth(:,3) = mod(paramout_fixwidth(:,3), 2*pi); 
    
    if plotting == 1
        plotfun_fit(9:1:19, wavelength, xin, yin, width, paramout_fixwidth, slopeflag, fixpeaks); 
    end
    

end