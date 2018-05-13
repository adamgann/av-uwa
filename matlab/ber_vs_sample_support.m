% ber_vs_sample_support.m
% 
% Plot the performance of the three different filtering algorithms as a
% function of sample support. 
%
% Adam Gannon, SUNY Buffalo, 2018



%clear all;
clear variables;
clc;
close all;

% Supress the nearly signular matrix warning, we know.
warning('off','MATLAB:nearlySingularMatrix')    

% Add functions/ to our current path
addpath('functions')

%% Settings

Params = struct;

% Debugging 
Params.checkConvergence = false;                                           % Plot norm(g) and |wmaxsinr - wav|^2 plots     
Params.plotScatter = false;                                                % Plot output of filter in scatterplot
Params.plotPsd = false;
Params.skipAvCalc = false;                                                  % Don't calculate AV filter.
Params.plotChanEst = false;

% Settings
Params.genColoredNoise = false;                                           % Use colored noise (vs AWGN)
Params.avSelectionMethod = 'SupervisedJdiv';  
Params.requirePlateau = false;
Params.maxSnapshots = false;                                                % If TRUE, use all packet bits to calculate autocorrelation matrix. 

%Channel Settings
Params.dynamicChannels = false;                                             % If TRUE, channel will change within the packet.
Params.estimateChan = false;


%% Parameters

Params.numSim = 10                                                           %Number of channel realizations                                                                      
                                                                            
                                                                            
Params.L=8;                                                                 % Code Length
Params.M=200;                                                               % Single-path
Params.Mrec = Params.M;

% SINR of the user and interferer(s)
Params.snrUser=80;                                                          
Params.snrInterf=85;                                                        

% Number of interferers to consider
Params.K=0;                                                            

% Target duration of the packet in seconds
Params.targetPacketDuration = 200e-3;


    
Params.av_maxiters = 200;


Params.fstart = 95e3;                                                        
Params.fstop = 145e3; 
Params.chan_files = ['acoustic_channel_simulator/channels/chan_wide_0.mat'; ...
                     'acoustic_channel_simulator/channels/chan_wide_1.mat'];


% Calc bandwidth and Tchip
Params.alpha = 0.35;
Params.sps = 4;
B = Params.fstop-Params.fstart
Params.Tc = (1+Params.alpha)/B


snapshotsVec = [1,4,8,10,20,40,50,100,200,500,1000,2e3,5e3,1e4,2e4,5e4];
Params.N=max(snapshotsVec); 


% Load the interfer channel matrix
if (Params.K)
    load(Params.chan_files(2,:))
    hmatInterf = hmat(1:Params.M,:);
end

% Load channel matrices and store to variables.
load(Params.chan_files(1,:))
hmatUser = hmat(1:Params.M,:);
Params.dt = dt;


%% Run Trials

% Storage Vectors
berMfStore = zeros(size(snapshotsVec));
berMvdrStore = zeros(size(snapshotsVec));
berAvStore = zeros(size(snapshotsVec));

sinrMfStore = zeros(size(snapshotsVec));
sinrMvdrStore = zeros(size(snapshotsVec));
sinrAvStore = zeros(size(snapshotsVec));

for jj=1:length(snapshotsVec)

    sprintf('Num Snapshots: %d',snapshotsVec(jj))
    Params.numSnapshots = snapshotsVec(jj);  
    
    if (Params.K)
        [ber_vec,sinr_vec] = run_one_trial(Params,hmatUser,hmatInterf);
    else
        [ber_vec,sinr_vec] = run_one_trial(Params,hmatUser);
    end
    
    % Update storage vectors
    berMfStore(jj) = ber_vec('MF');
    berMvdrStore(jj) = ber_vec('MaxSINR');
    berAvStore(jj) = ber_vec('AV');
    
    sinrMfStore(jj) = sinr_vec('MF');
    sinrMvdrStore(jj) = sinr_vec('MaxSINR');
    sinrAvStore(jj) = sinr_vec('AV');
   
    
end


%% 


% Plot SNR
figure;
semilogx(snapshotsVec,sinrMfStore,'LineWidth',2)
hold on
semilogx(snapshotsVec,sinrMvdrStore,'g-s','MarkerFaceColor','None','MarkerSize',9)
semilogx(snapshotsVec,sinrAvStore,'r-*','MarkerFaceColor','None','MarkerSize',9)
hold off

xH = xlabel('Data Record Size As Multiples of (L+M-1) Samples')
yH = ylabel('Pre-Detection SINR (dB)')
lH = legend('MF','MVDR','AV','location','NorthWest')

xlim([min(snapshotsVec) max(snapshotsVec)])
ylim([-10 30])

set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(lH,'Interpreter','latex','Fontsize',12)
set(xH,'Interpreter','latex','Fontsize',14)
set(yH,'Interpreter','latex','Fontsize',14)

