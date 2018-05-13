% sinr_vs_changain.m
%
% Calculate the output SINR of the AV filter as a function of channel gain.
% Store the thresholds required to reach an output SINR of 10dB. 
%
% Adam Gannon, SUNY Buffalo, 2018

clear variables;
clc;
close all;

warning('off','MATLAB:nearlySingularMatrix')                                % Supress the nearly signular matrix warning, we know.

addpath('functions')                                                        % This script uses functions. Add them to the PATH


%% parameters

Params = struct;

% Debugging 
Params.checkConvergence = false;                                           % Plot norm(g) and |wmaxsinr - wav|^2 plots     
Params.plotScatter = false;                                                % Plot output of filter in scatterplot
Params.plotPsd = false;
Params.skipAvCalc = false;                                                  % Don't calculate AV filter.
Params.plotChanEst = false;

% Settings
%Params.use_channel_sim_data = true;                                        % Use the simulator for the user of interest's channels
Params.genColoredNoise = true;                                            % Use colored noise (vs AWGN)
Params.avSelectionMethod = 'BlindJdiv';   %Dont use output sinr
Params.requirePlateau = false;
Params.maxSnapshots = true;                                                % If TRUE, use all packet bits to calculate autocorrelation matrix. 

%Channel Settings
Params.dynamicChannels = false;                                             % If TRUE, channel will change within the packet.
Params.estimateChan = false;


Params.numSim=100;                                                         %Number of channel realizations
Params.M=200;                                                               	% Single-path
Params.Mrec = Params.M;

Params.snrUser = 93.35;
Params.snrInterf=0;                                                         	% SINR of the interferers
Params.K=0;                                                               	 % Number of interferers

Params.av_maxiters = 100;



Params.fstart = 95e3;                                                        
Params.fstop = 145e3; 
Params.chan_files = ['acoustic_channel_simulator/channels/chan_wide_0.mat'; ...
                     'acoustic_channel_simulator/channels/chan_wide_1.mat'];



% Calc bandwidth and Tchip
Params.alpha = 0.35;
Params.sps = 4;
B = Params.fstop-Params.fstart
Params.Tc = (1+Params.alpha)/B


% Assumptions of channel
Params.targetPacketDuration = 200e-3;
ber_thresh = 1e-3;                                                          % Desired bit-error rate threshold

% Sweep parameters
codeVec = [8,12,16,24,32,64];
%code_vec = [8,16,64];
codeVec = [8,9,10,11,12,16,24,32,64];
codeVec = [8,9,10,11,12,14,16,20,24,28,32,48,64];

% Storage Vector
berMfStore = zeros(length(codeVec),1);
berMvdrStore = zeros(length(codeVec),1);
berAvStore = zeros(length(codeVec),1);

sinrMfStore = zeros(length(codeVec),1);
sinrMvdrStore = zeros(length(codeVec),1);
sinrAvStore = zeros(length(codeVec),1);

chanFitStore = zeros(3,length(codeVec));


%% Run simulation 

% Load the interfer channel matrix
if (Params.K)
    load(Params.chan_files(2,:))
    hmatInterf = hmat(1:Params.M,:);
end

% Load channel matrices and store to variables.
load(Params.chan_files(1,:))
hmatUser = hmat(1:Params.M,:);
Params.dt = dt;



for iCode=1:length(codeVec)
    
    % Set L and calculate num snapshots achievable within coherence
    current_L = codeVec(iCode)
    Params.L = current_L;
    
    disp('Current L:'); Params.L
    
    if (Params.K)
        [berVec,sinrOutVec] = run_one_trial(Params,hmatUser,hmatInterf);
    else
        [berVec,sinrOutVec,chanFit] = run_one_trial(Params,hmatUser);
    end

   

    % Update storage vectors
    berMfStore(iCode) = berVec('MF');
    berMvdrStore(iCode) = berVec('MaxSINR');
    berAvStore(iCode) = berVec('AV');
    
    sinrMfStore(iCode) = sinrOutVec('MF');
    sinrMvdrStore(iCode) = sinrOutVec('MaxSINR');
    sinrAvStore(iCode) = sinrOutVec('AV');

    chanFitStore(:,iCode) = chanFit;


end

%% Get the channel gain
hmat2=hmatUser.';
Gtilde=zeros(length(hmat2),1);
for countt= 1:size(hmat2,1)
    Gtilde(countt)= ((hmat2(countt, :)*hmat2(countt, :)'));
end
gtildeDb = 10*log10(Gtilde);
gainMin = min(gtildeDb);
gainMax = max(gtildeDb);

%% Plot Data
colorMat = get_color_spec();

format_vec = ['-*';'-o';'-v';'-^';'-<';'->';'-s';'-x';'-d'];
nMarkers = 10;
%inputVec = [gainMin:0.5:gainMax];
inputVec = linspace(gainMin,gainMax,200);

nSkip = ceil(length(inputVec)/nMarkers);

% What output SINR do we need to get our desired BER?
berThresh = 10; 
threshLine = berThresh.*ones(size(inputVec));

codeVecInd = [1:2:11]
codeVec(codeVecInd)


sinrThresh = zeros(1,size(chanFitStore,2));
lineHandles = zeros(1,length(codeVecInd));
figure; hold on
for jj = 1:length(codeVecInd)
    iCode = codeVecInd(jj);
    fitVec = chanFitStore(1,iCode)*inputVec+chanFitStore(2,iCode);
    diffVec = abs(fitVec-berThresh);
    [~,ind] = min(diffVec);
    sinrThresh(iCode)=inputVec(ind);
    
    lineHandles(jj) = plot(inputVec(1:nSkip:end),fitVec(1:nSkip:end),format_vec(jj,:));
    
end
plot(inputVec,threshLine,'k--');
hold off; box on

% Swap colors to be consistent
set(lineHandles(end-1),'color',colorMat(6,:));
set(lineHandles(end),'color',colorMat(5,:));

for hh=1:length(lineHandles)
    set(lineHandles(hh),'LineWidth',1.25);
    set(lineHandles(hh),'MarkerSize',9);
end

xH = xlabel('Channel Gain (dB)');
yH = ylabel('Pre-Detection SINR (dB)');
lH = legend('L=8','L=10','L=12','L=16','L=24','L=32','SINR Target','location','NorthWest');

set(lH,'Interpreter','latex','Fontsize',12);
set(xH,'Interpreter','latex','Fontsize',14);
set(yH,'Interpreter','latex','Fontsize',14);

xlim([-65 -50])
ylim([2 20])

%% Write the results to a text file

% Write the TX gain that was used
fileName = 'functions/decision_algorithms/chan_gain_thresholds.csv'
fId = fopen(fileName,'w');
fprintf(fId,'%f\n',Params.snrUser);
fclose(fId);

%Write the code-SINR pairs as a CSV
csvMat = [codeVec.', sinrThresh.']
dlmwrite(fileName,csvMat,'-append')


