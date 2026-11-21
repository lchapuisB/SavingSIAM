function [Flow_jl, Z_mean_jl, out_migration_rate, mean_wealth, mean_earn, wealth_to_earn, Lambda, a_bar] = compute_SS_moments(P, Xi)
    %% 1. INITIALIZATION %%
    J = P.J; 
    Nz = P.Nz; 
    delta = P.delta; 
    r = P.r;
    kappa = r / (1 + r);
    iz0 = round(Nz / 2); 

    % transition proba
    Pi_trans = zeros(Nz, Nz, J);
    for m = 1:J
        [~ , PI_temp] = discretize_normal(P.Nz, 0, P.rho, sqrt(P.sig2_eps(m)));
        Pi_trans(:, :, m) = PI_temp';  
    end
    
    % proba of moving 
    mu = zeros(J, J, Nz);
    log_Ez = zeros(J, Nz);
    for m = 1:J
        V_m = -P.gamma * Xi(m, :);
        V_max = max(V_m);
        for iz = 1:Nz
            E_val = sum(exp(V_m - V_max) .* Pi_trans(iz, :, m));
            log_Ez(m, iz) = V_max + log(E_val);
        end
    end
    
    B1hat = P.gamma * kappa * P.y_tilde - 0.5 * (kappa * P.gamma * P.w_tilde).^2 .* P.sig2_eta;
    
    for j = 1:J
        for iz = 1:Nz
            log_w = zeros(J, 1);
            for l = 1:J
                B1_jl = P.A(l) - P.psi(j, l) + B1hat(l);
                log_w(l) = P.nu * B1_jl - P.nu * log_Ez(l, iz);
            end
            w_norm = exp(log_w - max(log_w));
            mu(j, :, iz) = w_norm / sum(w_norm);
        end
    end
    
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
    
    %% 3. Computation of SS population and wealth
    
    is_disconnected = max(P.psi(:)) > 1e7;
    
    if is_disconnected
       
        Lambda = zeros(J, Nz);
        for j = 1:J
            [eig_z, ~] = eigs(Pi_trans(:, :, j)', 1, 'lr');
            stat_z = abs(real(eig_z));
            stat_z = stat_z / sum(stat_z);
            
            Lambda(j, :) = P.pop_share(j) * stat_z';
        end
        Lambda_vec = reshape(Lambda', 1, []);
    else
        
        opts_eigs.disp = 0;
        [eig_vec, ~] = eigs(T_mat_full', 1, 'lr', opts_eigs);
        Lambda_vec = abs(real(eig_vec))';
        Lambda_vec = Lambda_vec / sum(Lambda_vec);
        Lambda = reshape(Lambda_vec, [Nz, J])';
    end
    % --- INJECTION RICHESSE NOUVEAU-NÉS ---
    ArrivingNB = zeros(1, N_states);
    pop_j = sum(Lambda, 2); 
    total_deaths_j = delta * pop_j;
    
    for j = 1:J
        for l = 1:J
            mass_nb_to_l = total_deaths_j(j) * mu(j, l, iz0);
            for iz_new = 1:Nz
                s_prime = idx(l, iz_new);
                ArrivingNB(s_prime) = ArrivingNB(s_prime) + mass_nb_to_l * Pi_trans(iz0, iz_new, l);
            end
        end
    end
    
    C_income = zeros(1, N_states);
    for l = 1:J
        for iz_new = 1:Nz
            s_prime = idx(l, iz_new);
            C_income(s_prime) = P.y_tilde(l) + (1+r)*(P.w_tilde(l) * P.z_grid(iz_new) - Xi(l, iz_new));
        end
    end
    
    % Le vecteur B_A intègre l'arrivée d'agents déjà pourvus de a0 dans la destination
    B_A = Lambda_vec .* C_income + ArrivingNB .* P.a0;  
    
    A_vec = B_A / (speye(N_states) - T_mat);
    A_mass = reshape(A_vec, [Nz, J])';
    
    %% 4. Model's moments 
    a_bar = zeros(J, Nz);
    mean_wealth = zeros(J, 1);
    mean_earn = zeros(J, 1);
    
    for l = 1:J
        valid_idx = Lambda(l, :) > 1e-15;
        a_bar(l, valid_idx) = A_mass(l, valid_idx) ./ Lambda(l, valid_idx);
        
        mass_l = sum(Lambda(l, :));
        if mass_l > 0
            mean_wealth(l) = sum(A_mass(l, :)) / mass_l;
            mean_earn(l) = sum(P.w(l) * (1 + P.z_grid') .* Lambda(l, :)) / mass_l;
        end
    end
    
    wealth_to_earn = mean_wealth ./ mean_earn;
    
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