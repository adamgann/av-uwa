% full_loop_ct.m
%
% Simulate the adaptive algorithm and code length updates over a 
% time-varying channel. 
%
% Adam Gannon, SUNY Buffalo, 2018


clear variables;
clc;
close all;

warning('off','MATLAB:nearlySingularMatrix')                                % Supress the nearly signular matrix warning, we know.
addpath('functions')                                                        % This script uses functions. Add them to the PATH


%% Set Parameters

Params = struct;

% Debugging 
Params.checkConvergence = false;                                           % Plot norm(g) and |wmaxsinr - wav|^2 plots    
Params.plotPsd = false;
Params.plotScatter = false;
Params.plotChanEst = false;

% Channel Settings
Params.dynamicChannels = false;
Params.estimateChan = false;
Params.useFutureChan = false;

% AV Settings
Params.avSelectionMethod = 'BlindJdiv';
Params.requirePlateau = false;
Params.skipAvCalc = false;
Params.maxSnapshots = true;                                                % If TRUE, use all packet bits to calculate autocorrelation matrix. 

% Code Length Update Settings 
Params.decisionMethod = 'ChanGain';
Params.genColoredNoise = true;

% Adaptiation
Params.enableAdaptation = true;
Params.targetMargin = 1.5;
  

% Include the propagation delay calculated by the distance between Tx/Rx
Params.calcFeedbackDelay = true;


% Signal Model
Params.M=200;                                                               % Number of taps to use for channel simulation
Params.Mrec = 200;                                                           % Number of taps to use for chan estimation and filtering
Params.K=0;                                                               	% Number of interferers
Params.L = 24;
Params.snrInterf = 90;                                                     % Power level of interferers

Params.snrUser = 93.5;                                                         % Power level of user of interest

Params.av_maxiters = 100;

Params.fstart = 95e3;                                                        
Params.fstop = 145e3; 
 
Params.alpha = 0.35;
Params.sps = 4;

% Number of packets per simulation
nPkts = 2000;



%% Calculated parameters 

% Primary SINR
amplInterf = sqrt(10^(Params.snrInterf/10));
amplUser = sqrt(10^(Params.snrUser/10));

% Calc bandwidth and Tchip
B = Params.fstop-Params.fstart;
Params.Tc = (1+Params.alpha)/B;


% Init storage variables for previous SINR 
previousPrefilt = 0;
previousPostfilt = 0;


% Coherence time will dictate the packet length. We will set the packets
% such that their duration is slightly less than the coherence time.
Params.targetPacketDuration = 200e-3;


% Set up the temp variable for delayed L updates
nextL = Params.L;


% Pulse-shape filtering
if (exist('rcosdesign','file'))
    gT = rcosdesign(Params.alpha,6,Params.sps,'sqrt');
else
    gT = rcosine(1,Params.sps,'sqrt',Params.alpha);
end
    

% Load channel sim data
load acoustic_channel_simulator/channels/chan_wide_0.mat
hmatUser = hmat;
hmatUser = hmatUser(1:Params.M,:);      % Take the first M taps
Params.dt = dt;

% Load interference matrix
if (Params.K==1)
    load channel_sim/channels/chan_full_loop_0.mat
    hmatInterf = hmat(1:Params.M,:);
end



% Calculate the propagation delay for the feedback packets
% T_prop gives unidrectional delay, we want twice that. 
% At a fixed rate we don't need to wait for feedback packets before sending
% the next transmission, so rttTime is zero. 
rttTime = 0;
if ((Params.calcFeedbackDelay) && (Params.enableAdaptation))
    rttTime = 2*T_prop;
end
dtPerRoundTrip = round(rttTime/dt)
rttTime




% Elapsed time required before a code length update is processed. This is
% to take into account the prop time for the feedback packet.
timeReqVec = ([1:nPkts].*Params.targetPacketDuration)+rttTime; 




%% Simulation


% Storage

sinrAv= zeros(1,nPkts);
sinrPrefilt= zeros(1,nPkts);
berAv= zeros(1,nPkts);

packetDurStore = zeros(1,nPkts);
packetLenStore = zeros(1,nPkts);
chosenLStore = zeros(1,nPkts);
fullGainStore = zeros(1,nPkts);

% Index for keeping track of which channel vector we're using. 
% Start with the first column of the channel matrix 
iChanVec = 1; 

% Time elapsed will be advanced by packetDuration, but first we need to
% take into account the initial propagation time of the first packet
timeElapsed = (rttTime/2);


for iPkt = 1:nPkts

    %% Calculate packet duration parameters
    
    symbolDuration = Params.L*Params.Tc;          
    Params.N = floor(Params.targetPacketDuration / symbolDuration);                          % Adjust packet duration to coherence time
    Params.numSnapshots = Params.N;                                            %Use all packet bits
    packetDuration = (Params.N*Params.L*Params.Tc);                                   % Actual packet duration
    
    packetDurStore(iPkt) = packetDuration;
    packetLenStore(iPkt) = Params.N;
    
    % Number of chunks
    symbolsPerDt = floor(dt/symbolDuration);
    numDt = ceil(Params.N/symbolsPerDt);  
    
    %% Call the run_one_trial function
    
    Params.numSim = 1;
    [berOutVec,sinrOutVec] = run_one_trial(Params,hmatUser(:,iChanVec:end));
    hsVec = hmatUser(:,iChanVec);
        
    % Store the SNR and BER metrics
    berAv(iPkt) = berOutVec('AV');
    sinrAv(iPkt) = sinrOutVec('AV');
    sinrPrefilt(iPkt) = sinrOutVec('PREFILT');
        
    
    % Advance the channel index by the duration of the packet
    iChanVec = iChanVec + numDt;

    
    %% Calculate channel gain 
    
    % The channel gain of this current packet. 
    fullGainStore(iPkt) = 10*log10(hsVec'*hsVec);


    if (Params.useFutureChan)
        % Get the known future channel gain to get an upper bound on how good
        % adaptation can be. First calculate how many dt in the
        % future I should retrieve. 
        indAdvance = ceil(ceil(rttTime/dt)/numDt)*numDt;
        knownGainNext = 10*log10(hmatUser(:,iChanVec+indAdvance)'*hmatUser(:,iChanVec+indAdvance));
        chanGain = knownGainNext;
    else
        chanGain = fullGainStore(iPkt);
    end
    
    %% Decide Next L
    

    if (Params.enableAdaptation)
        [chosenLStore(iPkt)] = decide_next_code(Params,chanGain);
    else
        chosenLStore(iPkt) = Params.L;   % L doesn't change
    end
    
    
    %% Update Parameters for next run of simulation
    
    % Advance time elpased counter.
    % Current time is Tpacket + Tprop
    % We should see coherenceTime ~ packetDuration
    timeElapsed = timeElapsed + Params.targetPacketDuration;
    
    % Get the index of chosenLStore we should use to update Params.L
    % When rttTime = 0, iUpdate = iPkt
    % Delay the update of L to simulate prop delay. 
    diffReqVec = timeElapsed - timeReqVec;
    if any(diffReqVec >= 0)    %Make sure at least one req is met first
       diffReqVec(diffReqVec < 0) = inf; 
       [~,iUpdate] = min(diffReqVec);
       
       Params.L = chosenLStore(iUpdate);
    end
    
    
    disp(sprintf('Pkt: %d/%d\tCurrent L: %d\tSINR (pre/post): (%f/%f)\n',iPkt,nPkts,Params.L,sinrPrefilt(iPkt),sinrAv(iPkt)))
end

%% Calculate rate and data transfer

rateStore = (packetLenStore)./(1000*(packetDurStore)); %rate in kbps
rawDataTransfer = sum(packetDurStore.*rateStore);

%% Calc channel gain


lenHmat = size(hmatUser,2);
Gtilde = zeros(lenHmat,1);
for iChan = 1:lenHmat
   Gtilde(iChan) = hmatUser(:,iChan)'*hmatUser(:,iChan); 
end
Gtilde_db = 10*log10(Gtilde);


sinrTvec = cumsum(packetDurStore);
chanTvec = (0:length(Gtilde)-1)*dt; %Each value is dt apart, by chan sim



%% Plot figures

% Plot rate vs time
figure;
[hAx,hLine1,hLine2]=plotyy(sinrTvec,sinrPrefilt,sinrTvec,rateStore)
ylabel(hAx(1),'Prefilt SINR (dB)') % left y-axis 
ylabel(hAx(2),'Rate (kbps)') % right y-axis

% BER vs Time
figure;
[hAx,hLine1,hLine2]=plotyy(sinrTvec,sinrAv,sinrTvec,berAv)
ylabel(hAx(1),'Postfilt SINR (dB)') % left y-axis 
ylabel(hAx(2),'Bit Error Rate') % right y-axis

% Postfilt vs Channel Gain
figure;
subplot(2,1,1);
plot(sinrTvec,fullGainStore)
title('Channel Gain')
%xlim([0 sinrTvec(end)])
subplot(2,1,2);
plot(sinrTvec,sinrAv);
title('SINR AV')


%% Save the workspace

varSaveList = {'berAv','sinrAv','sinrPrefilt','sinrTvec','packetDurStore','rateStore','chanTvec','Gtilde_db','sinrTvec','fullGainStore'};


if (Params.enableAdaptation)
    save('workspace_adaptiveL.mat',varSaveList{:})
else
    fileName = sprintf('workspace_fixedL%d.mat',Params.L);
    save(fileName,varSaveList{:})
end



