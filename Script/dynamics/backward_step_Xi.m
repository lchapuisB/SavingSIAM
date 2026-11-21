function Xi_t = backward_step_Xi(Xi_t1, P_t, P_next, kappa_t1, Pi_all_next)
    factor = 1 / (1 + kappa_t1*(1+P_t.r));
    J = P_t.J;
    Nz = P_t.Nz;

    B1hat_next = compute_B1hat(P_next);

    log_Ez_mat = zeros(J, Nz);
    for m = 1:J
        v    = -P_next.gamma * Xi_t1(m, :);
        vmax = max(v);
        ev   = exp(v - vmax);
        log_Ez_mat(m, :) = vmax + log(ev * Pi_all_next{m});  
    end

    V_det = P_t.nu * ((P_next.A(:)' + B1hat_next(:)') - P_t.psi);

    log_Ez_3d = reshape(log_Ez_mat, 1, J, Nz);
    arg_exp_3d = V_det - P_t.nu * log_Ez_3d;  

    mx = max(arg_exp_3d, [], 2);                            
    ls = mx + log(mean(exp(arg_exp_3d - mx), 2));           
    ls_2d = squeeze(ls);                                     

    z_row = P_t.z_grid(:)';                                 
    rhs = log(P_t.beta*(1+P_t.r)) ...
          - P_t.gamma*kappa_t1*(1+P_t.r) * P_t.w_tilde(:) * z_row ...
          + P_t.A(:) + log(P_t.lambda_bar) - (1/P_t.nu)*ls_2d;

    Xi_t = -(1/P_t.gamma) * factor * rhs;
end