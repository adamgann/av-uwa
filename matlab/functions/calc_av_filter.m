function [ w_av ] = calc_av_filter( R_Rx, w_mf, params, opt )
% Runs the AV algorithm and chooses the best filter 

    
    term_at_convergence = false;                                         % Exit when norm(g) == 0
    desired_response = 1;                                               % We want wav'*y = 1
    [g_mat,W_av_mat] = av(R_Rx,w_mf,desired_response,term_at_convergence,params.av_maxiters);
    
    
    %% Take a fixed index every time
    if strcmp(params.avSelectionMethod,'Fixed')
       w_av = W_av_mat(:,params.av_ind); 
    
    
    %% Use convergence to ideal MMSE
    elseif strcmp(params.avSelectionMethod,'IdealMMSE')
        
        
        norm_vec = zeros(1,size(g_mat,2));
        for av_ind=1:size(g_mat,2)
            norm_vec(av_ind) = norm(W_av_mat(:,av_ind) - opt.w_sinr_ideal).^2;
        end
        [~,ind] = min(norm_vec);
        w_av = W_av_mat(:,ind);
        
        
        if (params.checkConvergence)
            figure;
            plot(norm_vec)
            xlabel('Algorithm Iteration')
            ylabel('||w_d - w_{mvdr} ||^2')
            title('MSE with Ideal MMSE')
        end
        
    %% Use minimum norm of aux vector
    elseif strcmp(params.avSelectionMethod,'NormG')
    
        g_norm = zeros(1,size(g_mat,2));
        for av_ind=1:size(g_mat,2)
            g_norm(av_ind) = norm(g_mat(:,av_ind));
        end
        [~,ind] = min(g_norm);
        w_av = W_av_mat(:,ind);
        
        if (params.checkConvergence)
            figure;
            plot(g_norm)
            xlabel('Algorithm Iteration')
            ylabel('||g||')
            title('Norm of G Vector')
        end
     
    %% MaxSINR
    elseif strcmp(params.avSelectionMethod,'OutputSINR')
        
        sinr_vec = zeros(1,size(g_mat,2));
        for av_ind=1:size(g_mat,2)
            wav_temp = W_av_mat(:,av_ind);
            sinr_vec(av_ind) = 10*log10(mean(abs(wav_temp'*opt.Y_S).^2)/mean(abs(wav_temp'*opt.Y_noise(:,1:size(opt.Y_S,2))).^2));
        end
        [~,ind] = max(sinr_vec);
        w_av = W_av_mat(:,ind);
        
        if (params.checkConvergence)
            figure;
            plot(sinr_vec)
            xlabel('Algorithm Iteration')
            ylabel('SINR (dB)')
            title('Post-Filtering SINR')
        end
        
        if (params.requirePlateau)
            sinr_diff = sinr_vec(end)-sinr_vec(end-100);
            if (sinr_diff>0.0001)
                error('No plateau')
            end
        end
        
        
    %% Blind J-divergence
    elseif strcmp(params.avSelectionMethod,'BlindJdiv')
        
        jdiv_store = zeros(1,size(g_mat,2));
        for av_ind=1:size(g_mat,2)
           rx_temp = W_av_mat(:,av_ind)'*opt.Y_Rx;
           rx_abs = abs(real(rx_temp));
           jdiv_store(av_ind) = (4*(mean(rx_abs).^2))/var(rx_abs);
        end
        [~,ind]=max(jdiv_store);
        w_av = W_av_mat(:,ind);
        
        if (params.checkConvergence)
            figure;
            plot(jdiv_store)
            xlabel('Algorithm Iteration')
            title('J-Divergence (Blind)')
        end
        
        
    %% Supervised J-divergence
    elseif strcmp(params.avSelectionMethod,'SupervisedJdiv')
        
        % Supervised J-Div
        jdiv_sup = zeros(1,size(g_mat,2));
        for av_ind=1:size(g_mat,2)
            rx_temp = W_av_mat(:,av_ind)'*opt.Y_Rx;
            rx_bit = opt.B_S.*real(rx_temp);
            jdiv_sup(av_ind) = (4*(mean(rx_bit).^2))/var(rx_bit);
        end
        [~,ind]=max(jdiv_sup);
        w_av = W_av_mat(:,ind);
        
        
        if (params.checkConvergence)
            figure;
            plot(jdiv_sup)
            xlabel('Algorithm Iteration')
            title('J-Divergence (Supervised)')
        end 
        
        if (params.requirePlateau)
            j_diff = jdiv_sup(end)-jdiv_sup(end-10);
            if (j_diff>0.001)
                error('No plateau')
            end
        end
        
    %% CV-MOV
    elseif strcmp(params.avSelectionMethod,'MOV')
        
        cmov_vec = zeros(1,params.av_maxiters);
        for iteration_index = 1:params.av_maxiters;    % n=1:params.av_maxiters

            variance_vec = zeros(1,params.num_snapshots);
            for exclude_index = 1:params.num_snapshots;  % j=1:params.num_snapshots

                % Calculate the R/j matrix 
                copy_inds = [1:exclude_index-1 exclude_index+1:params.num_snapshots];
                Y_excl = opt.Y_Rx(:,copy_inds);
                R_excl = Y_excl*Y_excl';

                [~,W_av_mat] = av(R_excl,w_mf,desired_response,term_at_convergence,iteration_index+1);
                wav_temp = W_av_mat(:,iteration_index); %Chose the n-th filter (also the last) 
                variance_val = wav_temp'*(opt.Y_Rx*opt.Y_Rx')*wav_temp;
                variance_vec(exclude_index) = variance_val;
            end
           cmov_vec(iteration_index) = sum(variance_vec);
        end
        [~,cvmov_ind]=min(cmov_vec);
        w_av = W_av_mat(:,cvmov_ind);
        
        if (params.checkConvergence)
            figure;
            plot(abs(cmov_vec))
            xlabel('Algorithm Iteration')
            title('CV-MOV')
        end  
      
    %% Else
        
    else
        error('No entry for that algorithm')
    end
    
    
    % Hold if checking convergence
    if (params.checkConvergence)
        %pause
        error('Stopping to check convergence')
    end

end

