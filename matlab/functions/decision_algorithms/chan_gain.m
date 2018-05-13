function [ nextL ] = chan_gain( currentChanGain, Params)
% Using the BER curves previously calculated, we will select the highest
% rate (smallest L) of which we can fulfil the SINR requirement. 


    % Load the SINR thresholds previously calculated 
    % These indicate the prefilt SINR required to maintain BER > 1e-3
    persistent sinrThresholds codeLengths txGainRef;
    if isempty(sinrThresholds)
        
        fId = fopen('chan_gain_thresholds.csv');
       	firstLine = fgetl(fId);
        txGainRef = str2num(firstLine);
        fclose(fId);
        
        csvMat = csvread('chan_gain_thresholds.csv',1);
        codeLengths = csvMat(:,1);
        sinrThresholds = csvMat(:,2);
    end
    
    % Adjust to the different TX gain
    %txGainRef = 95; % The thresholds are calculated for 75dB
    gainDiff = txGainRef - Params.snrUser;
    sinrAdjusted = sinrThresholds + gainDiff;
    

    % Add a margin in dB
    sinrMargin = sinrAdjusted + Params.targetMargin;
    
    % Find the shortest L supported by current SINR value
    % This will be the first positive number in diffVec
    diffVec = currentChanGain - sinrMargin;
    sinrInd = find((diffVec>0),1);
    
    % Use the longest code for lowest SINR
    if isempty(sinrInd)
        sinrInd = length(sinrThresholds);
    end
    nextL = codeLengths(sinrInd);
    
    
    
end

