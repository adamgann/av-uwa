% ber_vs_sinr_av.m
%
% Sweep TX gain and calculate BER for different code lengths using both AV
% and MVDR filters. 
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


Params.numSim=1000;                                                         %Number of channel realizations
Params.M=200;                                                               	% Single-path
Params.Mrec = Params.M;

Params.snrInterf=80;                                                         	% SINR of the interferers
Params.K=0;                                                               	% Number of interferers

Params.av_maxiters = 100;


% Frequency and channel settings 
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
codeVec = [8,16,32];

sinrInVec = [85,86,87,88,89,90:0.5:100];

% Storage Vectors
sinrOutMf = zeros(length(codeVec),length(sinrInVec));
sinrOutMvdr = zeros(length(codeVec),length(sinrInVec));
sinrOutAv = zeros(length(codeVec),length(sinrInVec));
sinrOutPrefilt = zeros(length(codeVec),length(sinrInVec));


berOutMf = zeros(length(codeVec),length(sinrInVec));
berOutMvdr = zeros(length(codeVec),length(sinrInVec));
berOutAv = zeros(length(codeVec),length(sinrInVec));


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



for code_ind=1:length(codeVec)
    
    % Set L and calculate num snapshots achievable within coherence
    current_L = codeVec(code_ind)
    Params.L = current_L;

    

    for sinr_ind=1:length(sinrInVec)
        
        disp('Current SINR: ');sinrInVec(sinr_ind)
        disp('Current L:'); Params.L
        Params.snrUser = sinrInVec(sinr_ind);

        if (Params.K)
            [berVec,sinrOutVec] = run_one_trial(Params,hmatUser,hmatInterf);
        else
            [berVec,sinrOutVec] = run_one_trial(Params,hmatUser);
        end



        % Update storage vectors
        berOutMf(code_ind,sinr_ind) = berVec('MF');
        berOutMvdr(code_ind,sinr_ind) = berVec('MaxSINR');
        berOutAv(code_ind,sinr_ind) = berVec('AV');

        sinrOutMf(code_ind,sinr_ind) = sinrOutVec('MF');
        sinrOutMvdr(code_ind,sinr_ind) = sinrOutVec('MaxSINR');
        sinrOutAv(code_ind,sinr_ind) = sinrOutVec('AV');
        sinrOutPrefilt(code_ind,sinr_ind) = sinrOutVec('PREFILT');
    

    end

end

%% Plot Results

colorMat = get_color_spec;
format_vec = ['b-*';'g-o';'k-v';'r-^';'c->';'y-<']
markerInd = 1:2:length(sinrInVec);

figure;
lineH(1) = semilogy(sinrInVec,berOutAv(1,:),'-*','color',colorMat(1,:));
hold on
lineH(2) = semilogy(sinrInVec,berOutAv(2,:),'-^','color',colorMat(4,:));
lineH(3) = semilogy(sinrInVec,berOutAv(3,:),'->','color',colorMat(5,:));
lineH(4) = semilogy(sinrInVec,berOutMvdr(1,:),':*','color',colorMat(1,:));
lineH(5) = semilogy(sinrInVec,berOutMvdr(2,:),':^','color',colorMat(4,:));
lineH(6) = semilogy(sinrInVec,berOutMvdr(3,:),':>','color',colorMat(5,:));

% Format line properties
for hh=1:length(lineH)
   set(lineH(hh),'MarkerIndices',markerInd);
   set(lineH(hh),'LineWidth',1.0);
   set(lineH(hh),'MarkerSize',9);
end

% Plot a reference line 
refline = 1e-4.*ones(size(sinrInVec));
plot(sinrInVec,refline,'k--','LineWidth',1)

% And label
xlabel('Transmit Symbol Energy (dB re $\mu$Pa)','Interpreter','latex','FontSize',14)
ylabel('Bit Error Rate','Interpreter','latex','FontSize',14)
legend('AV, L=8','AV, L=16','AV, L=32','MVDR, L=8', 'MVDR, L=16','MVDR, L=32','location','SouthWest')
ylim([5e-5 1e-0])

return 
%% Other useful plots 

%close all

format_vec = ['b-*';'g-o';'k-v';'r-^';'c->';'y-<']


% Plot Prefilt as a function of TX power
figure;
plot(sinrInVec,sinrOutPrefilt(1,:),'b-*')
hold on
for ii=2:size(sinrOutPrefilt,1)
    semilogy(sinrInVec,sinrOutPrefilt(ii,:),format_vec(ii,:))
end
legend(string(codeVec),'location','SouthEast')
xlabel('TX Gain (dB)')
ylabel('Prefiltering SINR (dB)')


% Plot AV SINR as a function of TX power
figure;
plot(sinrInVec,sinrOutAv(1,:),'b-*')
hold on
for ii=2:size(sinrOutAv,1)
    semilogy(sinrInVec,sinrOutAv(ii,:),format_vec(ii,:))
end
legend(string(codeVec),'location','SouthEast')
xlabel('TX Gain (dB)')
ylabel('Postfiltering (AV) SINR (dB)')






% Plot AV BER as a function of input TX power 
figure;
semilogy(sinrInVec,berOutAv(1,:),'b-*')
hold on
for ii=2:size(berOutAv,1)
    semilogy(sinrInVec,berOutAv(ii,:),format_vec(ii,:))
end
refLine = berThresh.*ones(size(sinrInVec));
semilogy(sinrInVec,refLine,'k:','LineWidth',2)
legend(string(codeVec),'location','SouthWest')
xlabel('TX Gain (dB)')
ylabel('Bit Error Rate (AV)')
ylim([5e-5 1e-1])


% Plot AV BER as a function of postfilt SINR
figure;
semilogy(sinrOutAv(1,:),berOutAv(1,:),'b-*')
hold on
for ii=2:size(berOutAv,1)
    semilogy(sinrOutAv(ii,:),berOutAv(ii,:),format_vec(ii,:))
end
legend(string(codeVec),'location','SouthEast')
xlabel('Postfilt SINR (dB)')
ylabel('Bit Error Rate (AV)')




% Plot AV BER as a function of Prefilt SINR
figure;
semilogy(sinrOutPrefilt(1,:),berOutAv(1,:),'b-*')
hold on
for ii=2:size(berOutAv,1)
    semilogy(sinrOutPrefilt(ii,:),berOutAv(ii,:),format_vec(ii,:))
end
legend(string(codeVec),'location','SouthEast')
xlabel('Prefilt SINR (dB)')
ylabel('Bit Error Rate (AV)')



% Plot MF BER as a function of input TX power 
figure;
semilogy(sinr_in_vec,ber_out_mf(1,:),'b-*')
hold on
for ii=2:size(ber_out_mf,1)
    semilogy(sinr_in_vec,ber_out_mf(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('TX Gain (dB)')
ylabel('Bit Error Rate (MF)')


% Plot MF BER as a function of postfilt SINR 
figure;
semilogy(sinr_out_mf(1,:),ber_out_mf(1,:),'b-*')
hold on
for ii=2:size(ber_out_mf,1)
    semilogy(sinr_out_mf(ii,:),ber_out_mf(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('Postfiltering SINR (MF) (dB)')
ylabel('Bit Error Rate (MF)')


% Plot MaxSINR BER as a function of TX power
figure;
semilogy(sinr_in_vec,ber_out_maxsinr(1,:),'b-*')
hold on
for ii=2:size(ber_out_av,1)
    semilogy(sinr_in_vec,ber_out_maxsinr(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('TX Gain (dB)')
ylabel('Bit Error Rate (MMSE)')

% Plot MF SINR as a function of TX Gain
figure;
semilogy(sinr_in_vec,sinr_out_mf(1,:),'b-*')
hold on
for ii=2:size(ber_out_av,1)
    semilogy(sinr_in_vec,sinr_out_mf(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('TX Gain (dB)')
ylabel('Post-Filtering SINR (MF)')

% Plot MMSE SINR as a function of TX Gain
figure;
semilogy(sinr_in_vec,sinr_out_maxsinr(1,:),'b-*')
hold on
for ii=2:size(ber_out_av,1)
    semilogy(sinr_in_vec,sinr_out_maxsinr(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('TX Gain (dB)')
ylabel('Post-Filtering SINR (MMSE)')


% Plot MMSE BER as a function of Postfilt SINR
figure;
semilogy(sinr_out_maxsinr(1,:),ber_out_maxsinr(1,:),'b-*')
hold on
for ii=2:size(ber_out_av,1)
    semilogy(sinr_out_maxsinr(ii,:),ber_out_maxsinr(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('Post-Filtering SINR - MMSE (db)')
ylabel('Bit Error Rate (MMSE)')

% Plot MF BER as a function of Postfilt SINR
figure;
semilogy(sinr_out_mf(1,:),ber_out_mf(1,:),'b-*')
hold on
for ii=2:size(ber_out_av,1)
    semilogy(sinr_out_mf(ii,:),ber_out_mf(ii,:),format_vec(ii,:))
end
legend(string(code_vec),'location','SouthEast')
xlabel('Post-Filtering SINR - MF (db)')
ylabel('Bit Error Rate (MF)')
