function Xi_cf_local = compute_Xi_cf_local(P_base, P_cf, Xi_base_guess)
% Computes the local counterfactual continuation value for the shift-share decomposition.
% For each state j, it assumes only state j receives the P_cf parameters,
% while all other states remain at P_base parameters. It then solves the
% general equilibrium Xi and extracts the j-th row.

    J = P_base.J;
    Nz = P_base.Nz;
    Xi_cf_local = zeros(J, Nz);
    
    fprintf('  Computing Local Counterfactual Xi state by state...\n');
    for j = 1:J
        P_temp = P_base;
        P_temp.w(j) = P_cf.w(j);
        P_temp.A(j) = P_cf.A(j);
        P_temp.w_tilde(j) = P_cf.w_tilde(j);
        P_temp.y_tilde(j) = P_cf.y_tilde(j);
        P_temp.sig2_eps(j) = P_cf.sig2_eps(j);
        P_temp.sig2_eta(j) = P_cf.sig2_eta(j);
        
        Xi_temp = compute_Xi(P_temp, 1e-7, Xi_base_guess);
        
        Xi_cf_local(j, :) = Xi_temp(j, :);
        
        if mod(j, 10) == 0
            fprintf('    Computed %d/%d states\n', j, J);
        end
    end
end
