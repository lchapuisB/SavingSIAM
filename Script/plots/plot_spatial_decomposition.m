function plot_spatial_decomposition(P_0, P_1, Xi_0, Xi_1, output_dir, scenario_name)
    % Executes an exact shift-share decomposition on Per-Capita Average Wealth
    % utilizing the Time-Reversed Markov Chain.
    
    J = P_0.J; Nz = P_0.Nz; N = J * Nz;
    I_mat = speye(N);
    
    [T_0, C_0, NB_0, L_0, mu_0, Pi_0] = extract_model_operators(P_0, Xi_0);
    [T_1, C_1, NB_1, L_1, mu_1, Pi_1] = extract_model_operators(P_1, Xi_1);
    
    % Counterfactual Xi: Autarky Continuation Value
    P_aut_0 = P_0;
    P_aut_0.psi = 1e8 * ones(P_0.J, P_0.J);
    P_aut_0.psi(1:P_0.J+1:end) = 0; % Zero on diagonal
    Xi_aut_0 = compute_Xi(P_aut_0, 1e-7, Xi_0);
    [T_aut_0, ~, ~, L_aut_0, ~, ~] = extract_model_operators(P_aut_0, Xi_aut_0);
    
    P_aut_1 = P_1;
    P_aut_1.psi = 1e8 * ones(P_1.J, P_1.J);
    P_aut_1.psi(1:P_1.J+1:end) = 0;
    Xi_aut_1 = compute_Xi(P_aut_1, 1e-7, Xi_1);
    
    [F_Dam_0,  F_Risk_0,  F_Spat_0]  = decompose_flows(P_0, Xi_0, Xi_aut_0, mu_0, Pi_0);
    [F_Dam_1,  F_Risk_1,  F_Spat_1]  = decompose_flows(P_1, Xi_1, Xi_aut_1, mu_1, Pi_1);
    
    % Per-capita savings flows
    c_0 = C_0 + NB_0 ./ L_0;
    
    % Intensive Margin Sub-components (Change in behavioural flows)
    dc_Dam  = F_Dam_1  - F_Dam_0;
    dc_Risk = F_Risk_1 - F_Risk_0;
    dc_Spat = F_Spat_1 - F_Spat_0;
    
    % Time-Reversed Markov Matrices
    T_tilde_0 = spdiags(1./L_0, 0, N, N) * T_0' * spdiags(L_0, 0, N, N);
    T_tilde_1 = spdiags(1./L_1, 0, N, N) * T_1' * spdiags(L_1, 0, N, N);
    T_tilde_aut_0 = spdiags(1./L_aut_0, 0, N, N) * T_aut_0' * spdiags(L_aut_0, 0, N, N);
    
    % Baseline conditional per-capita wealth
    a_0 = (I_mat - T_tilde_0) \ c_0;
    
    % Per-capita savings flows in 2100
    c_1 = C_1 + NB_1 ./ L_1;
    a_1 = (I_mat - T_tilde_1) \ c_1;

    idx_or = (38 - 1) * Nz + (1:Nz);
    % 2. Get Oregon's baseline population shares for each z-state
    pi_0_or = L_0(idx_or) / sum(L_0(idx_or));
    % 3. Calculate the pure local behavioral effect for Damage Dampening
    pure_effect_Dam = sum(dc_Dam(idx_or) .* pi_0_or);
    fprintf('Pure Local Damage Dampening for Oregon: %f\n', pure_effect_Dam);
    
    % 1. Intensive: Behavioural Savings (Evaluated at AUTARKY network)
    % This guarantees they perfectly match the 'no_migration' scenario
    a_Dam  = (I_mat - T_tilde_aut_0) \ dc_Dam;
    a_Risk = (I_mat - T_tilde_aut_0) \ dc_Risk;
    
    % Spatial Option absorbs the network premium so the sum identity holds
    dc_total = dc_Dam + dc_Risk + dc_Spat;
    a_Spat = ((I_mat - T_tilde_0) \ dc_total) - a_Dam - a_Risk;
    
    % 2. Spatial Transition: Change in network history (Evaluated at NEW behavior)
    a_Trans = a_1 - ((I_mat - T_tilde_0) \ c_1);
    
    %% Aggregation to State-Level Local Averages
    idx = @(j, iz) (j - 1) * Nz + iz;
    
    da_Dam_bar = zeros(J, 1); da_Risk_bar = zeros(J, 1); da_Spat_bar = zeros(J, 1);
    da_Trans_bar = zeros(J, 1); da_Comp_bar = zeros(J, 1);
    
    for j = 1:J
        % State population totals
        L0_j = sum(L_0(idx(j, 1:Nz)));
        L1_j = sum(L_1(idx(j, 1:Nz)));
        
        for iz = 1:Nz
            s = idx(j, iz);
            pi_0 = L_0(s) / L0_j;
            pi_1 = L_1(s) / L1_j;
            
            % Aggregated Micro Responses (weighted by BASELINE demographic shares)
            da_Dam_bar(j)   = da_Dam_bar(j)   + pi_0 * a_Dam(s);
            da_Risk_bar(j)  = da_Risk_bar(j)  + pi_0 * a_Risk(s);
            da_Spat_bar(j)  = da_Spat_bar(j)  + pi_0 * a_Spat(s);
            da_Trans_bar(j) = da_Trans_bar(j) + pi_0 * a_Trans(s);
            
            % Local Compositional Sorting (Change in shares * NEW wealth)
            da_Comp_bar(j)  = da_Comp_bar(j)  + (pi_1 - pi_0) * a_1(s);
        end
    end
    
    Total_da = da_Dam_bar + da_Risk_bar + da_Spat_bar + da_Trans_bar + da_Comp_bar;
    
    % Display & Save CSV
    T_res = table(P_0.state_names(:), Total_da, da_Dam_bar, da_Risk_bar, da_Spat_bar, da_Trans_bar, da_Comp_bar, ...
        'VariableNames', {'State', 'Total_da', 'Income_Dampening', 'Precautionary_Risk', 'Spatial_Option', 'Spatial_Transition', 'Compositional'});
    
    csv_path = fullfile(output_dir, sprintf('exact_percapita_decomposition_%s.csv', scenario_name));
    writetable(T_res, csv_path);
    fprintf('  -> Saved per-capita decomposition table to: %s\n', csv_path);
    
    %% --- Plot Stacked Bar Chart ---
    [~, sort_idx] = sort(Total_da, 'descend');
    
    idx_pos = sort_idx(1:min(5, J));
    idx_neg = sort_idx(max(1, J-4):J);
    idx_top = unique([idx_pos; idx_neg], 'stable');
    top_n = length(idx_top);
    
    fig = figure('Visible', 'off', 'Position', [150 150 1000 600]);
    bar_data = [da_Dam_bar(idx_top), da_Risk_bar(idx_top), da_Spat_bar(idx_top), da_Trans_bar(idx_top), da_Comp_bar(idx_top)];
    
    b = bar(bar_data, 'stacked');
    b(1).FaceColor = [0.82, 0.42, 0.35]; % Terracotta 
    b(2).FaceColor = [0.85, 0.70, 0.30]; % Muted Gold 
    b(3).FaceColor = [0.35, 0.68, 0.48]; % Emerald 
    b(4).FaceColor = [0.45, 0.80, 0.60]; % Light Emerald 
    b(5).FaceColor = [0.25, 0.52, 0.72]; % Cobalt Blue 
    
    xticks(1:top_n);
    xticklabels(P_0.state_names(idx_top));
    xtickangle(45);
    
    hold on;
    plot(1:top_n, Total_da(idx_top), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'DisplayName', 'Total \Delta \bar{a}');
    
    ylabel('Per-Capita Average Wealth');
    legend({'Income Dampening', 'Precautionary Risk', 'Spatial Option', 'Spatial Transition', 'Compositional', 'Total'}, ...
           'Location', 'southwest', 'FontSize', 12, 'Box', 'off');
    grid on;
    
    graph_path = fullfile(output_dir, sprintf('exact_percapita_decomposition_bar_%s.png', scenario_name));
    saveas(fig, graph_path);
    fprintf('  -> Saved per-capita decomposition graph to: %s\n', graph_path);
end
%% Helper Function 1: Extract Model Operators
function [T_mat, C_income, NB_vec, L_vec, mu, Pi_trans] = extract_model_operators(P, Xi)
    J = P.J; Nz = P.Nz; delta = P.delta; r = P.r; kappa = r / (1 + r); iz0 = round(Nz / 2);
    
    Pi_trans = zeros(Nz, Nz, J);
    for m = 1:J
        [~ , PI_temp] = discretize_normal(Nz, 0, P.rho, sqrt(P.sig2_eps(m)));
        Pi_trans(:, :, m) = PI_temp;  
    end
    
    mu = zeros(J, J, Nz); log_Ez = zeros(J, Nz);
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
    
    N_states = J * Nz; num_elements = J * Nz * J * Nz;
    row_idx = zeros(num_elements, 1); col_idx = zeros(num_elements, 1); vals_T = zeros(num_elements, 1);
    idx = @(j, iz) (j - 1) * Nz + iz; counter = 1;
    
    for j = 1:J
        for iz = 1:Nz
            s = idx(j, iz);
            for l = 1:J
                prob_move = mu(j, l, iz) * (1 - delta);
                for iz_new = 1:Nz
                    s_prime = idx(l, iz_new);
                    row_idx(counter) = s; col_idx(counter) = s_prime;
                    vals_T(counter)  = prob_move * Pi_trans(iz, iz_new, l);
                    counter = counter + 1;
                end
            end
        end
    end
    T_mat = sparse(row_idx, col_idx, vals_T, N_states, N_states);
    
    if max(P.psi(:)) > 1e7
        Lambda = zeros(J, Nz);
        for j = 1:J
            [eig_z, ~] = eigs(Pi_trans(:, :, j)', 1, 'lr');
            stat_z = abs(real(eig_z)); stat_z = stat_z / sum(stat_z);
            Lambda(j, :) = P.pop_share(j) * stat_z';
        end
        L_vec = reshape(Lambda', 1, [])';
    else
        row_idx_d = zeros(num_elements, 1); col_idx_d = zeros(num_elements, 1); vals_d = zeros(num_elements, 1);
        counter_d = 1;
        for j = 1:J
            for iz = 1:Nz
                s = idx(j, iz);
                for l = 1:J
                    prob_move_nb = mu(j, l, iz0) * delta;
                    for iz_new = 1:Nz
                        s_prime = idx(l, iz_new);
                        row_idx_d(counter_d) = s; col_idx_d(counter_d) = s_prime;
                        vals_d(counter_d) = prob_move_nb * Pi_trans(iz0, iz_new, l);
                        counter_d = counter_d + 1;
                    end
                end
            end
        end
        T_death = sparse(row_idx_d, col_idx_d, vals_d, N_states, N_states);
        opts_eigs.disp = 0;
        [eig_vec, ~] = eigs((T_mat + T_death)', 1, 'lr', opts_eigs);
        L_vec = abs(real(eig_vec)); L_vec = L_vec / sum(L_vec);
    end
    
    ArrivingNB = zeros(N_states, 1);
    Lambda = reshape(L_vec, [Nz, J])';
    pop_j = sum(Lambda, 2); total_deaths_j = delta * pop_j;
    
    for j = 1:J
        for l = 1:J
            mass_nb_to_l = total_deaths_j(j) * mu(j, l, iz0);
            for iz_new = 1:Nz
                s_prime = idx(l, iz_new);
                ArrivingNB(s_prime) = ArrivingNB(s_prime) + mass_nb_to_l * Pi_trans(iz0, iz_new, l);
            end
        end
    end
    NB_vec = ArrivingNB .* P.a0;
    
    C_income = zeros(N_states, 1);
    for l = 1:J
        for iz_new = 1:Nz
            s_prime = idx(l, iz_new);
            C_income(s_prime) = P.y_tilde(l) + (1+r)*(P.w_tilde(l) * P.z_grid(iz_new) - Xi(l, iz_new));
        end
    end
end

%% Helper Function 2: Analytic Decomposition of Savings Flow
function [F_Dam, F_Risk, F_Spat] = decompose_flows(P, Xi, Xi_aut, mu, Pi_trans)
    % Xi     = equilibrium Xi (full model)
    % Xi_aut = autarky Xi (no spatial mobility)
    J = P.J; Nz = P.Nz; N = J * Nz;
    F_Dam = zeros(N, 1); F_Risk = zeros(N, 1); F_Spat = zeros(N, 1);
    idx = @(j, iz) (j - 1) * Nz + iz;

    for j = 1:J
        V_m = -P.gamma * Xi(j, :);
        V_max = max(V_m);
        for iz = 1:Nz
            s = idx(j, iz);

            E_val = sum(exp(V_m - V_max) .* Pi_trans(iz, :, j));
            log_Ez = V_max + log(E_val);
            
            E_Xi = sum(Pi_trans(iz, :, j) .* Xi(j, :));
            E_Xi_aut = sum(Pi_trans(iz, :, j) .* Xi_aut(j, :));

            % 1. Damage Dampening (Income & Mean Reversion)
            % Cleaned of spatial option value by using Xi_aut
            F_Dam(s) = P.y_tilde(j)/(1+P.r) + P.w_tilde(j)*P.z_grid(iz) - E_Xi_aut;

            % 2. Risk Premium (Second Moment)
            % The + E_Xi isolates the variance (Jensen's wedge)
            Risk_Pers = (1/P.gamma) * log_Ez + E_Xi;
            Risk_Trans = (P.r^2 * P.gamma / (2 * (1+P.r)^2)) * (P.w_tilde(j)^2 * P.sig2_eta(j));
            F_Risk(s) = Risk_Pers + Risk_Trans;

            % 3. Spatial Option Value
            % Current option value + expected future option value
            F_Spat(s) = (1 / (P.gamma * P.nu)) * log(mu(j, j, iz)) - (E_Xi - E_Xi_aut);
        end
    end
end