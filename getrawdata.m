function [HistTot_array, XUV_only, stageTimes, freqAxis] = getrawdata(folderString, delayscan_start, wavelength)

    % data location
    %     folderString = '/Users/annaliw/code/KrCO2_scan/'; 
    saveString = 'extracted_data'; 

    % get data files
    dataList = dir(fullfile(folderString, '*.mat'));
    numberSubScans = length(dataList);

    % load in all the files 
    numFiles = numberSubScans; 
    for ii=1:1:numFiles
        creationCheck = exist('HistTot'); 
        if creationCheck == 0 
            makeArrays = 1; 
        else
            makeArrays = 0; 
        end
        currentSubScan = dataList(ii); 
        load(char(string(currentSubScan.folder) + '/' + string(currentSubScan.name))); 
        %SigSum is sum of anode signal
        %HistTot is the set of histograms
        %stage_positions are the stage positions for each scan
        if makeArrays == 1  
            SigSum_array = zeros([size(SigSum), numFiles]); 
            HistTot_array = zeros([size(HistTot), numFiles]); 
            stagePositions_array = zeros([size(stage_positions), numFiles]); 
        end
        % saving each subscan data in an array
        SigSum_array(:,:,ii) = SigSum; 
        HistTot_array(:,:,ii) = HistTot; 
        stagePositions_array(:,:,ii) = stage_positions; 
    end
    
    % convert to pump/probe delay time
    stage_times = (stage_positions*2*1e-3)/(2.9979e8); % tmp variable
    stageCenter = stage_times(round(length(stage_times)/2)); % find center
    stageTimes_array = (stagePositions_array*2*1e-3)/(2.9979e8) - stageCenter; % full array
    
    % clip out data where XUV and IR are not overlapped
    XUV_only = HistTot_array(:,end,:); 
    HistTot_array = HistTot_array(:,delayscan_start:(end-1),:); 
    stageTimes = squeeze(stageTimes_array(:,delayscan_start:(end-1),1)); 
    
    % make fourier space axis
    df = 1/(stageTimes(end)-stageTimes(1)); 
    freqAxis = (-(numel(stageTimes)-1)/2:1:(numel(stageTimes)-1)/2)*df*(wavelength*10^(-9))/(3*10^8); 
    

end