function [ sigVec ] = get_signal_dynamic_chan( sigUser, bitUser, amplUser, chanMat, Params, numDt, symbolsPerDt )
% Creates a TX signal vector for the case that the channel is changing

    sigVec = [];
    iChan = 1;

    for tt=1:numDt

        % Cut the bit vector into a chunk
        % If not the first chunk, pad the vector by nSymbolsPad symbols
        % Let's make nSymbolsPad half the chunk length, by default
        nSymbolsPad = ceil(symbolsPerDt/2);
        ind_start = (tt-1)*symbolsPerDt+1;
        if (tt ~= 1)
            ind_start = ind_start -nSymbolsPad;
        end
        ind_stop = min(tt*symbolsPerDt,Params.N);  
        ind_vec = ind_start:ind_stop;
        B_chunk = bitUser(ind_vec);




        % Get the channel corresponding to this time index 
        hs_vec = chanMat(:,iChan);

        % Generate the CT signal for this time instance
        % Reset the filter prototype for the first chunk
        % Pad the last chunk with zeros to flush filter
        resetFilt = (tt==1);
        padFilt = (tt==numDt);
        sigChunk = get_signal_static_chan(sigUser, B_chunk, amplUser, hs_vec, Params, resetFilt, padFilt);

        % Cut the padded symbol from the front of the filtered chunk and
        % append
        if (tt ~= 1)
            nCut = nSymbolsPad*(Params.L*Params.sps)+1;
            sigChunk = sigChunk(nCut:end);
        end
        sigVec = [sigVec; sigChunk];

        % Advance the index of the channel matrix
        iChan = iChan + 1;


    end

end

