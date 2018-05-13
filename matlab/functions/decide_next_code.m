function [ nextL ] = decide_next_code( Params, currentChanGain)
% Decision process for selecting next L
% Top level file calls algos in decision_algorithms

    

    addpath('functions/decision_algorithms')  
    
    if strcmp(Params.decisionMethod,'ChanGain')
        [nextL] = chan_gain(currentChanGain,Params);
        
    else
        error('No entry for that decision method.')
    end
    

end

