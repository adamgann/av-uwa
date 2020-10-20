function [ berVec, sinrVec, chanFit ] = run_one_trial( Params , hmatUser, hmatInterf, interfOn )
% Get one datapoint of the simulation.
%
% Adam Gannon, SUNY Buffalo, 2018

%% Set up input / persistent parameters

persistent masterCodeSeq gT;

if isempty(masterCodeSeq)
    % In this sim, the signatures are static Gold sequences, only the length is
    % varied. Generate from comm toolbox, if exists. Otherwise load .mat 
    if (exist('comm.GoldSequence'))
        hGold = comm.GoldSequence;
        hGold.SamplesPerFrame = 64;
        masterCodeSeq = step(hGold);
        masterCodeSeq(masterCodeSeq == 0) = -1;
    else
      load masterCodeSeq.mat
    end
end

if isempty(gT)
    if (exist('rcosdesign'))
        g_T = rcosdesign(Params.alpha,6,Params.sps,'sqrt');
    else
        g_T = rcosine(1,Params.sps,'sqrt',Params.alpha);
    end
end



%% Grab parameters from struct

numSim = Params.numSim;
K = Params.K;

% Unpack params for ease of use
M = Params.M;
L_M = Params.L+M-1;

% Calculated parameters
amplUser=sqrt(10^(Params.snrUser/10));
amplInterf=sqrt(10^(Params.snrInterf/10));

% Storage Vectors
berMf= zeros(1,numSim);
berMaxsinr = zeros(1,numSim);
berAv = zeros(1,numSim);

sinrMf = zeros(1,numSim);
sinrMvdr = zeros(1,numSim);
sinrAv = zeros(1,numSim);
sinrPrefilt = zeros(1,numSim);

chanGainStore = zeros(1,numSim);


% Iterator for channels.
iChanVec = 1;



iRun=1; %Iteration index
while iRun<(numSim+1)

    %% Calculate parameters for this iteration
    if (numSim>1)
        iRun
    end

    % Set duration and snapshots dynamically based on the target packet
    % duration and the symbol duration
    symbolDuration = Params.L*Params.Tc;
    if (Params.maxSnapshots)

        Params.N = floor(Params.targetPacketDuration / symbolDuration);                          % Adjust packet duration to coherence time
        packetDuration = (Params.N*Params.L*Params.Tc);                                   % Actual packet duration
        Params.numSnapshots = Params.N;                                            %Use all packet bits
    end

    % Number of chunks
    symbolsPerDt = floor(Params.dt/symbolDuration);
    numDt = ceil(Params.N/symbolsPerDt);

    % Primary (Interferers)
    sigInterf=(1/sqrt(Params.L)).*sign(randn(Params.L,Params.K));
    bitsInterf=sign(randn(Params.K,Params.N));

    % Secondary (User of Interest)
    sigUser = masterCodeSeq(1:Params.L)./sqrt(Params.L);
    bitsUser=sign(randn(1,Params.N));  % Data bits

    %% Create the TX signal vector

    if (Params.dynamicChannels)

        % Get the TX signal over a chunk of the channel matrix.
        hs_mat = hmatUser(:,iChanVec:iChanVec+numDt-1);
        sigVec = get_signal_dynamic_chan(sigUser, bitsUser, amplUser, hs_mat, Params, numDt, symbolsPerDt);
        iChanVec = iChanVec + numDt;

    else

        % Get the TX signal over only one channel iteration
        hs_mat = hmatUser(:,iChanVec);
        sigVec = get_signal_static_chan(sigUser, bitsUser, amplUser, hs_mat, Params);
        iChanVec = iChanVec + numDt;

    end


    %% Channel Simulation


    % Generate interference vector
    interfVec = zeros(size(sigVec));
    for kk=1:K
        if (interfOn)
            interfVec = interfVec + get_signal_static_chan(sigInterf(:,kk),bitsInterf(kk,:), amplInterf, hmatInterf(:,1), Params);
        end
    end

    % Generate noise
    if (Params.genColoredNoise)
        Noise=get_colored_noise(Params.L,M,Params.N+10,Params);
        noiseVec = reshape(Noise,[],1);
        noiseVec = noiseVec(1:length(sigVec));
    else
        Noise=(randn(L_M,Params.N)+1i*randn(L_M,Params.N))/sqrt(2);
        noiseVec = (randn(size(sigVec))+1i*randn(size(sigVec)))/sqrt(2);
    end

    % Add noise and interference to signal vector
    disturbanceVec = interfVec + noiseVec;
    rxVec = sigVec + disturbanceVec;

    % Pulse-matched filtering
    symFilt = conv(sigVec, g_T);
    interfFilt = conv(interfVec, g_T);
    rxFilt = conv(rxVec, g_T);


    % Deterministic Symbol Timing and Downsample
    sampInd = ((length(g_T))); %Going through filter twice, delays by half number of taps per filter. %FIXME: Which samp ind??
    symSampled = symFilt(sampInd:Params.sps:end);
    interfSampled = interfFilt(sampInd:Params.sps:end);
    rxSampled = rxFilt(sampInd:Params.sps:end);


    % Create multipath matrices
    sigMat = vec_to_multipath_mat(symSampled,Params,Params.Mrec);
    interfMat = vec_to_multipath_mat(interfSampled,Params,Params.Mrec);
    rxMat = vec_to_multipath_mat(rxSampled,Params,Params.Mrec);


    % Calculate RN matrix
    ySamples = rxMat(:,1:Params.numSnapshots);                                      % Limit the number of snapshots used to calculate R
    R_Rx=(ySamples*ySamples')/size(ySamples,2);                          % Calculate the sample-limited R



    %% Receive signal and estimate channel


    % Channel-processed signature matrix
    SS_mat = toeplitz([sigUser;zeros(Params.Mrec-1,1)],[sigUser(1),zeros(1,Params.Mrec-1)]);        %Multipath-processed signature matrix

    % Estimate the channel using a fixed number of known bits
    if (Params.estimateChan)
        nTrain = 100; % Keep this constant for all code lengths
        yTrain = rxMat(:,1:nTrain);
        bTrain = bitsUser(1:nTrain).';
        hestVec = (pinv(SS_mat)*yTrain*bTrain)/(nTrain*amplUser);
        hs_vec = hestVec;
    else
        % Use the first channel instance for the entire packet
        hs_vec = hs_mat(1:Params.Mrec,1);
    end

    chanGainStore(iRun) = 10*log10(hs_vec'*hs_vec);

    if (Params.plotChanEst)
        figure;
        hold on;
        plot(real(hs_mat(1:Params.Mrec,1)),'r')
        plot(real(hs_mat(1:Params.Mrec,end)),'k')
        plot(real(hestVec))
        legend('Start','End','Est')
    end


    if (Params.plotPsd)
        Nfft = 1024;
        fs=5e6;

        y_vec = reshape(Noise,numel(Noise),1);  %Reshape into vector
        y_cut = y_vec(1:Nfft*floor(length(y_vec)/Nfft));
        y_split=reshape(y_cut,Nfft,[]);
        Y_split = fft(y_split,Nfft);

        Pk = mean( abs(Y_split).^2, 2) / Nfft;
        Pk_wattsHz = Pk / fs;
        Pk_dbHz = 10*log10(Pk_wattsHz);
        Pk_ss = Pk_dbHz(1:floor(end/2));


        y_vec = reshape(Noise,numel(Noise),1);  %Reshape into vector
        y_cut = y_vec(1:Nfft*floor(length(y_vec)/Nfft));
        y_split=reshape(y_cut,Nfft,[]);
        Y_split = fft(y_split,Nfft);

        Pk = mean( abs(Y_split).^2, 2) / Nfft;
        Pk_wattsHz = Pk / fs;
        Pk_dbHz = 10*log10(Pk_wattsHz);
        Pk_ds = Pk_dbHz(1:floor(end/2));

        deltaf = fs/Nfft;
        fx = (0:Nfft/2-1)*deltaf;


        y_split=reshape(sigMat(:,1:Nfft*floor(Params.N/Nfft)),Nfft,[]);

        Y_split = fft(y_split,Nfft);

        Sk = mean( abs(Y_split).^2, 2) / Nfft;
        Sk_wattsHz = Sk / fs;
        Sk_dbHz = 10*log10(Sk_wattsHz);
        Sk_ss = Sk_dbHz(1:floor(end/2));


        figure;
        plot(fx,Pk_ss,'r','LineWidth', 2)
        hold on
        plot(fx,Pk_ds,'g')
        plot(fx,Sk_ss)
        legend('Vector','CT','UOI')
        title('PSD of Interference Plus Noise')


        return
    end

    %% Calculate Filters



    % Matched filter
    %SS_mat = toeplitz([sigUser;zeros(M-1,1)],[sigUser(1),zeros(1,M-1)]);        %Multipath-processed signature matrix
    w_mf = (SS_mat*hs_vec)/norm(SS_mat*hs_vec).^2;


    % MaxSINR, built with SMI autocorrelation matrix
    w_sinr = (R_Rx\w_mf)/(w_mf'*(R_Rx\w_mf));



    % Calc AV filter
    if (Params.skipAvCalc)
        w_av = zeros(size(w_sinr));
    else
        opt = struct;
        opt.Y_Rx = rxMat;
        opt.Y_noise = Noise;
        opt.Y_S = sigMat;
        opt.B_S = bitsUser;
        %opt.w_sinr_ideal = (R_ideal\w_mf)/(w_mf'*(R_ideal\w_mf));

        try
            w_av = calc_av_filter(R_Rx,w_mf,Params,opt);
        catch ME
            if strcmp(ME.message,'Stopping to check convergence')
                return
            elseif strcmp(ME.message,'No plateau')
                disp('Skipping...')
                continue
            else
                ME.message
                return
            end

        end
    end



    % Now that convergence has been tested, for convergence
    % both filters must be normalized to get proper bit energy
    w_av = w_av./(norm(SS_mat*hs_vec).^2);
    w_sinr = w_sinr./(norm(SS_mat*hs_vec).^2) ;


    %% Bit and SINR Recovery

    % Discard the last few bits of the packet for calculating BER. This is
    % to account for responses of filters and other simulation effects
    nBitsToCut = 5;
    bitsUser = bitsUser(1:end-nBitsToCut);


    % Recover matched filter bits
    rxMf = (w_mf'*rxMat)/amplUser;
    rxMf = rxMf(1:end-nBitsToCut);

    bitRecMf=sign(real(rxMf));
    errvecMf = bitsUser ~= bitRecMf;
    berMf(iRun) = nnz(errvecMf)/length(errvecMf);



    % Recover maxSINR bits
    rxMaxsinr = (w_sinr'*rxMat)/amplUser;
    rxMaxsinr = rxMaxsinr(1:end-nBitsToCut);
    bitRecMaxsinr=sign(real(rxMaxsinr));

    errvecMaxsinr = bitsUser ~= bitRecMaxsinr;
    berMaxsinr(iRun) = nnz(errvecMaxsinr)/length(errvecMaxsinr);


    % Recover AV bits
    rxAv = (w_av'*rxMat)/amplUser;
    rxAv = rxAv(1:end-nBitsToCut);
    bitRecAv=sign(real(rxAv));


    errvecAv = bitsUser ~= bitRecAv;
    berAv(iRun)=nnz(errvecAv)/length(errvecAv);


    evmMf = sqrt(mean(abs(rxMf-bitsUser).^2));
    evmMvdr = sqrt(mean(abs(rxMaxsinr-bitsUser).^2));
    evmAv = sqrt(mean(abs(rxAv-bitsUser).^2));

    sinrMf(iRun) = 10*log10(1/evmMf.^2);
    sinrMvdr(iRun) = 10*log10(1/evmMvdr.^2);
    sinrAv(iRun) = 10*log10(1/evmAv.^2);

    % Prefiltering SINR
    sigPower = sum(abs(sigVec).^2)/length(sigVec);
    distPower = sum(abs(disturbanceVec).^2)/length(disturbanceVec);
    sinrPrefilt(iRun) = 10*log10(sigPower/distPower);




    % Plot scatterplots
    if (Params.plotScatter)

        if exist('scatterplot')
            scatterplot(rxMf(1:100))
            title('MF')

            scatterplot(rxMaxsi(1:100))
            title('MaxSINR')

            scatterplot(rxAv(1:100))
            title('AV')
        else

            figure;scatter(real(rxMf(1:100)),imag(rxMf(1:100)))
            title('MF')

            figure;scatter(real(rxAv(1:100)),imag(rxAv(1:100)))
            title('AV')


        end

    end


iRun = iRun + 1;   % Increase index
end

% Return the average BER and SINR for each filter


sinrMfStore = mean(sinrMf);
sinrMvdrStore = mean(sinrMvdr);
sinrAvStore= mean(sinrAv);
sinrPrefiltStore = mean(sinrPrefilt);


berMfStore = mean(berMf);
berMvdrStore = mean(berMaxsinr);
berAvStore = mean(berAv);

keySet = {'MF','MaxSINR','AV','PREFILT'};
sinrSet = [sinrMfStore,sinrMvdrStore,sinrAvStore,sinrPrefiltStore];
berSet = [berMfStore, berMvdrStore, berAvStore,];

berVec = containers.Map(keySet(1:3),berSet);
sinrVec = containers.Map(keySet,sinrSet);


%% Fit channel gain to output SINR

% Get the polynomial fit
chanFit = polyfit(chanGainStore,sinrAv,1);

% Calculate the 95% confidence interval
if exist('polyfit')
    [chanPoly,S,mu] = polyfit(chanGainStore,sinrAv,1);
    [fitLine,chanDelta] = polyval(chanPoly,chanGainStore,S,mu);
    deltaMean = 2*mean(chanDelta);
else
    deltaMean = 0;
end
chanFit = [chanFit, deltaMean];


% Plot the fit and the bounds
if (0)
   figure;
   scatter(chanGainStore,sinr_av)
   hold on
   plot(chanGainStore,fitLine,'r')
   plot(chanGainStore,(fitLine-chanDelta),'k-')
   plot(chanGainStore,(fitLine+chanDelta),'k-')
end

return
