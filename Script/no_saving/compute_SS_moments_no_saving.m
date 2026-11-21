function [Flow_jl, Z_mean_jl, out_migration_rate, mean_wealth, mean_earn, wealth_to_earn, Lambda, a_bar] = compute_SS_moments_no_saving(P, Xi)
    %% 1. INITIALIZATION %%
    J = P.J; 
    Nz = P.Nz; 
    delta = P.delta; 
    iz0 = round(Nz / 2); 

    % transition proba
    Pi_trans = zeros(Nz, Nz, J);
    Pi_all = cell(J, 1);
    for m = 1:J
        [~ , PI_temp] = discretize_normal(P.Nz, 0, P.rho, sqrt(P.sig2_eps(m)));
        Pi_trans(:, :, m) = PI_temp;  
        Pi_all{m} = PI_temp;
    end
    
    % proba of moving using the no-saving function
    mu = compute_mu_no_saving(Xi, P, Pi_all);
    
    %% 2. Transition matrix (Survivors and Newborns) %%
    N_states = J * Nz;
    num_elements = J * Nz * J * Nz;
    
    row_idx = zeros(num_elements, 1);
    col_idx = zeros(num_elements, 1);
    vals_T  = zeros(num_elements, 1);
    
    row_idx_d = zeros(num_elements, 1);
    col_idx_d = zeros(num_elements, 1);
    vals_d    = zeros(num_elements, 1);
    
    idx = @(j, iz) (j - 1) * Nz + iz;
    
    counter = 1;
    counter_d = 1;
    
    for j = 1:J
        for iz = 1:Nz
            s = idx(j, iz);
            for l = 1:J
                prob_move = mu(j, l, iz) * (1 - delta);
                prob_move_nb = mu(j, l, iz0) * delta;
                
                for iz_new = 1:Nz
                    s_prime = idx(l, iz_new);
                    
                    row_idx(counter) = s;
                    col_idx(counter) = s_prime;
                    vals_T(counter)  = prob_move * Pi_trans(iz, iz_new, l);
                    counter = counter + 1;
                    
                    row_idx_d(counter_d) = s;
                    col_idx_d(counter_d) = s_prime;
                    vals_d(counter_d) = prob_move_nb * Pi_trans(iz0, iz_new, l);
                    counter_d = counter_d + 1;
                end
            end
        end
    end
    
    T_mat = sparse(row_idx, col_idx, vals_T, N_states, N_states);
    T_death = sparse(row_idx_d, col_idx_d, vals_d, N_states, N_states);
    T_mat_full = T_mat + T_death;
    
    %% 3. Computation of SS population
    
    is_disconnected = max(P.psi(:)) > 1e7;
    
    if is_disconnected
        Lambda = zeros(J, Nz);
        for j = 1:J
            [eig_z, ~] = eigs(Pi_trans(:, :, j)', 1, 'lr');
            stat_z = abs(real(eig_z));
            stat_z = stat_z / sum(stat_z);
            Lambda(j, :) = P.pop_share(j) * stat_z';
        end
    else
        opts_eigs.disp = 0;
        [eig_vec, ~] = eigs(T_mat_full', 1, 'lr', opts_eigs);
        Lambda_vec = abs(real(eig_vec))';
        Lambda_vec = Lambda_vec / sum(Lambda_vec);
        Lambda = reshape(Lambda_vec, [Nz, J])';
    end
    
    %% 4. Model's moments (No Saving adjustments)
    a_bar = zeros(J, Nz);
    mean_wealth = zeros(J, 1);
    mean_earn = zeros(J, 1);
    
    for l = 1:J
        mass_l = sum(Lambda(l, :));
        if mass_l > 0
            mean_earn(l) = sum(P.w(l) * (1 + P.z_grid') .* Lambda(l, :)) / mass_l;
        end
    end
    
    wealth_to_earn = zeros(J, 1);
    
    Flow_jl = zeros(J, J);         
    Z_mean_jl = zeros(J, J);        
    
    for j = 1:J
        for l = 1:J
            flow_vec = Lambda(j, :) .* squeeze(mu(j, l, :))';
            Flow_jl(j, l) = sum(flow_vec);
            if Flow_jl(j, l) > 1e-15
                Z_mean_jl(j, l) = sum(flow_vec .* P.z_grid') / Flow_jl(j, l);
            else
                Z_mean_jl(j, l) = 0;
            end
        end
    end
    
    pop_j = sum(Lambda, 2);
    out_flow_j = sum(Flow_jl, 2) - diag(Flow_jl);
    out_migration_rate = out_flow_j ./ pop_j;
end
