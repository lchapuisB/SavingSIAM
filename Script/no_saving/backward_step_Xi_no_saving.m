function Xi_t = backward_step_Xi_no_saving(Xi_t1, P_t, P_next, Pi_all_next)
% Computes one backward step of Xi for the NO-SAVING counterfactual.
% Uses P_t for the flow utility and P_next for the continuation.

    J = P_t.J;
    Nz = P_t.Nz;
    gam_t = P_t.gamma;
    beta_t = P_t.beta;
    lambda_bar_t = P_t.lambda_bar;
    
    A_t   = P_t.A(:); 
    w_til_t = P_t.w_tilde(:);
    y_til_t = P_t.y_tilde(:);
    sig2_eta_t = P_t.sig2_eta(:);
    z_grid_t = P_t.z_grid(:)'; 

    gam_next = P_next.gamma;
    nu_next  = P_next.nu;
    psi_next = P_next.psi; 

    F_exp = zeros(J, Nz);
    for ell = 1:J
        F_exp(ell, :) = -gam_t * y_til_t(ell) - gam_t * w_til_t(ell) * z_grid_t - A_t(ell) + 0.5 * (gam_t * w_til_t(ell))^2 * sig2_eta_t(ell);
    end

   
    nu_B_ns = -nu_next * psi_next; 

    Xi_t = zeros(J, Nz);
    log_Ez = zeros(J, Nz);

    % 1. Compute log( E_{z'|z,m}[exp(-gam_next*Xi_t1(m, z'))] ) using t+1 transition
    for m = 1:J
        V_m = -gam_next * Xi_t1(m, :);
        V_max = max(V_m);
        E_val_vec = exp(V_m - V_max) * Pi_all_next{m}; 
        log_Ez(m, :) = V_max + log(E_val_vec);
    end

    % 2. Compute the new Xi_t for each origin ell and past shock iz
    for iz = 1:Nz

        arg_exp = nu_B_ns - nu_next * log_Ez(:, iz)'; 
        
      
        max_arg = max(arg_exp, [], 2);
        mean_exp = mean(exp(arg_exp - max_arg), 2);
        log_sum_m = max_arg + log(mean_exp);
        
       
        log_beta_Omega = log(beta_t) + log(lambda_bar_t) - (1/nu_next) * log_sum_m; 
        
     
        flow_val = F_exp(:, iz);
        cont_val = log_beta_Omega; 
        
  
        max_val = max(flow_val, cont_val);
        log_sum_terms = max_val + log(exp(flow_val - max_val) + exp(cont_val - max_val));
        
      
        Xi_t(:, iz) = -(1/gam_t) * log_sum_terms;
    end
end
