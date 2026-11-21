function results = run_scenario(expectation_type, P,counterfactual)
    set_aesthetics();
% ============================================================
% Forward Projection: 2010 -> 2100 Climate Shock Scenario
% ============================================================
script_dir = fileparts(mfilename('fullpath'));
set(0, 'DefaultFigureVisible', 'off');
addpath(fullfile(script_dir, '..', 'plots'), fullfile(script_dir, '..', 'dynamics'), fullfile(script_dir, '..', 'helpers'), fullfile(script_dir, '..', 'static'), fullfile(script_dir, '..', 'no_saving'));

%% ── STEP 1 : BASELINE (2010) AND TERMINAL (2100) STATES ─────────────────
fprintf('\n=== Step 1: Defining Baseline (2010) and Terminal (2100) States ===\n');

P_2010 = P; 
if strcmp(counterfactual, 'no_migration')
    P_2010.psi = 1e8 * ones(P.J, P.J); 
    for j = 1:P.J
        for l = 1:P.J
            if j == l
                P_2010.psi(j, l) = 0;
            end
        end
    end
    P_2100 = P_2010;
else 
    P_2100 = P_2010; 
end

wage_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_wage_2080.csv');
opts = detectImportOptions(wage_damages_path);
damages_wage = readtable(wage_damages_path, opts);

amen_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_amenity_2080.csv');
opts = detectImportOptions(amen_damages_path);
damages_amen = readtable(amen_damages_path, opts);

housing_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_housing_2080.csv');
opts = detectImportOptions(housing_damages_path);
damages_housing = readtable(housing_damages_path, opts);

shock_housing =damages_housing.shock_housing;
shock_amen = damages_amen.shock_amenity;
shock_wage = damages_wage.shock_wage;

P_2100.w = P_2010.w .* (1 + shock_wage);
P_2100.w_tilde = P_2010.w_tilde .* (1 + shock_wage);
P_2100.q = P_2010.q .* (1 + shock_housing);
P_2100.A = P_2010.A - shock_amen;
P_2100.y_tilde = P_2100.w_tilde - P_2100.q;


P_trans_start = P_2010;
if isfield(P_2010, 'uncertainty')
    if P_2010.uncertainty == 1
        P_trans_start.sig2_eta = P_2010.sig2_eta .* (1 + 50 * abs(shock_wage));
        P_2100.sig2_eta = P_2010.sig2_eta;
    elseif P_2010.uncertainty == 2
        P_trans_start.sig2_eta = P_2010.sig2_eta;
        P_2100.sig2_eta = P_2010.sig2_eta .* (1 + 50 * abs(shock_wage));
    end
end




%% ── STEP 2 : STATIONARY XI (2010 AND 2100) ──────────────────────────────
fprintf('\n=== Step 2: Computing Terminal Steady State Xi ===\n');
fprintf('  Xi SS 2010 (Baseline)...\n');
if strcmp(counterfactual, 'no_saving_struct')
    Xi_SS_2010 = compute_Xi_no_saving(P_2010);
    fprintf('  Xi SS 2100 (Terminal)...\n');
    Xi_SS_2100 = compute_Xi_no_saving(P_2100);
else
    Xi_SS_2010 = compute_Xi(P_2010);
    fprintf('  Xi SS 2100 (Terminal)...\n');
    Xi_SS_2100 = compute_Xi(P_2100);
end

iz_mid   = round(P_2010.Nz / 2);
diff_SS  = Xi_SS_2100(:, iz_mid) - Xi_SS_2010(:, iz_mid);
fprintf('  Xi_2100 - Xi_2010 at z=0: min=%.5f  max=%.5f  mean=%.5f\n', ...
        min(diff_SS), max(diff_SS), mean(diff_SS));


%% ── STEP 3 : STATIONARY DISTRIBUTIONS (2010) ────────────────────────────
fprintf('\n=== Step 3: Baseline Distributions and Initial Wealth ===\n');
if strcmp(counterfactual, 'no_saving_struct')
    [~, ~, ~, mean_w_2010, mean_e_2010, wye_2010, Lambda_SS_2010, a_bar_SS_2010] = ...
        compute_SS_moments_no_saving(P_2010, Xi_SS_2010);

    [~, ~, ~, mean_w_2100, mean_e_2100, wye_2100, Lambda_SS_2100, a_bar_SS_2100] = ...
        compute_SS_moments_no_saving(P_2100, Xi_SS_2100);
        
    delta_w = P_2100.w_tilde(:) - P_2010.w_tilde(:);
    kappa_stat = P_2010.r / (1 + P_2010.r);
    Xi_cf_wage = Xi_SS_2010 + kappa_stat * (delta_w * P_2100.z_grid(:)');
    
    [~, ~, ~, ~, ~, ~, ~, a_bar_cf] = ...
        compute_SS_moments_no_saving(P_2100, Xi_cf_wage);
        
    Xi_cf_local = compute_Xi_cf_local(P_2010, P_2100, Xi_SS_2010);
    [~, ~, ~, ~, ~, ~, ~, a_bar_cf_local] = ...
        compute_SS_moments_no_saving(P_2100, Xi_cf_local);
else
    [~, ~, ~, mean_w_2010, mean_e_2010, wye_2010, Lambda_SS_2010, a_bar_SS_2010] = ...
        compute_SS_moments(P_2010, Xi_SS_2010);

    [~, ~, ~, mean_w_2100, mean_e_2100, wye_2100, Lambda_SS_2100, a_bar_SS_2100] = ...
        compute_SS_moments(P_2100, Xi_SS_2100);
        
    delta_w = P_2100.w_tilde(:) - P_2010.w_tilde(:);
    kappa_stat = P_2010.r / (1 + P_2010.r);
    Xi_cf_wage = Xi_SS_2010 + kappa_stat * (delta_w * P_2100.z_grid(:)');
    
    [~, ~, ~, ~, ~, ~, ~, a_bar_cf] = ...
        compute_SS_moments(P_2100, Xi_cf_wage);
        
    Xi_cf_local = compute_Xi_cf_local(P_2010, P_2100, Xi_SS_2010);
    [~, ~, ~, ~, ~, ~, ~, a_bar_cf_local] = ...
        compute_SS_moments(P_2100, Xi_cf_local);
end

fprintf('  W/Y 2010 = %.3f\n', sum(mean_w_2010.*sum(Lambda_SS_2010,2)) / sum(mean_e_2010.*sum(Lambda_SS_2010,2)));
fprintf('  W/Y 2100 = %.3f (Targeting post-shock equilibrium)\n', sum(mean_w_2100.*sum(Lambda_SS_2100,2)) / sum(mean_e_2100.*sum(Lambda_SS_2100,2)));



%% ── STEP 4 : BUILD TRANSITION PATH (2010 -> 2100) ───────────────────────
fprintf('\n=== Step 4: Building transition path ===\n');

TRANSITION_TYPE = 'linear'; 
T_shock = 90;  
T_sim   = 500; 

paths_temp = build_paths(P_trans_start, P_2100, T_shock, TRANSITION_TYPE);

keys = {'w','w_tilde','y_tilde','sig2_eps','sig2_eta', 'A', 'q'};
paths = struct(); 

for k = 1:numel(keys)
    terminal_state = paths_temp.(keys{k})(:, end);
    steady_state_padding = repmat(terminal_state, 1, T_sim - T_shock);
    paths.(keys{k}) = [paths_temp.(keys{k}), steady_state_padding];
end

fprintf('  T_shock=%d, T_sim=%d, type=%s\n', T_shock, T_sim, TRANSITION_TYPE);

T = T_sim;
fprintf('  T=%d, type=%s\n', T, TRANSITION_TYPE);


%% ── STEP 5 : SOLVE EXPECTATIONS ──────────────────────────
if strcmp(expectation_type, 'myopia')
    fprintf('\n=== Step 5: Myopic expectations (Solving stationary Xi_t at each t) ===\n');
    if strcmp(counterfactual, 'no_saving_struct')
        error('Myopia not yet implemented for no_saving_struct');
    end
    Xi_trans = solve_transition_myopia(P_2010, Xi_SS_2010, paths, T);
else
    fprintf('\n=== Step 5: Backward induction (Solving Xi_t) ===\n');
    if strcmp(counterfactual, 'no_saving_struct')
        Xi_trans = solve_transition_backward_no_saving(P_2100, Xi_SS_2100, paths, T);
    else
        Xi_trans = solve_transition_backward(P_2100, Xi_SS_2100, paths, T);
    end
end



%% ── STEP 6 : FORWARD SIMULATION (POPULATION & WEALTH) ───────────────────  
fprintf('\n=== Step 6: Forward simulation of Lambda ===\n');
if strcmp(counterfactual, 'no_saving_struct')
    Lambda_trans = simulate_Lambda_forward_no_saving(P_2010, Lambda_SS_2010, Xi_trans, paths, T);
else
    Lambda_trans = simulate_Lambda_forward(P_2010, Lambda_SS_2010, Xi_trans, paths, T);
end


%% ── STEP 7 : DYNAMIC MOMENTS ────────────────────────────────────────────
fprintf('\n=== Step 7: Computing dynamic moments ===\n');
if strcmp(counterfactual, 'no_saving_struct')
    [agg_earn, agg_wealth, agg_cons, agg_mobility, mean_earn_t, mean_wealth_t, WY_t, mean_ce_t, ~] = ...
        compute_moments_dyn_no_saving(Lambda_trans, Xi_trans, paths, P_2010);
else
    [agg_earn, agg_wealth, agg_cons, agg_mobility, mean_earn_t, mean_wealth_t, WY_t, mean_ce_t, ~] = ...
        compute_moments_dyn(Lambda_trans, Xi_trans, a_bar_SS_2010, paths, P_2010);
end


%% ── STEP 8 : FIGURES ────────────────────────────────────────────────────
fprintf('\n=== Step 8: Plotting Transition ===\n');

if P.uncertainty == 1
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_uncertainty_damage_path';
elseif P.uncertainty == 2
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_uncertainty_variance_shocks';
elseif strcmp(expectation_type, 'myopia')      
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_myopia';
elseif strcmp(counterfactual, 'no_migration')
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_no_migration';
elseif strcmp(counterfactual, 'big_migration')
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_big_migration';
elseif strcmp(counterfactual, 'no_saving_struct')
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_no_saving_struct';
else
    output_dir = '/Users/leon/Documents/SavingSIAM/Output/graph_perfect_foresight';
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

plot_xi_transition(P_2010, P_2100, T, TRANSITION_TYPE, Xi_trans, Xi_SS_2010, Xi_SS_2100, paths, output_dir);

plot_dynamic_moments(T, TRANSITION_TYPE,  agg_earn, agg_wealth, WY_t, agg_mobility, output_dir);

[delta_CE_mean, delta_CE_var, delta_W_var] = plot_welfare_variance(P_2010, P_2100, Lambda_SS_2010, Lambda_SS_2100, a_bar_SS_2010, a_bar_SS_2100, Xi_SS_2010, Xi_SS_2100, output_dir);

plot_welfare_transition(T, mean_earn_t, mean_wealth_t, mean_ce_t, output_dir);

% Spatial Decomposition
plot_spatial_decomposition(P_2010, P_2100, Xi_SS_2010, Xi_SS_2100, output_dir, sprintf('%s_no_counterfactual', counterfactual));
% Welfare Decomposition

%% ── STEP 9 : SAVE RESULTS ───────────────────────────────────────────────
fprintf('\n=== Step 9: Saving Results ===\n');

results = struct();
if strcmp(counterfactual, 'no_saving_struct')
    CE_1 = Xi_SS_2010;
else
    kappa = P_2010.r / (1 + P_2010.r);
    CE_1 = kappa * a_bar_SS_2010 + Xi_SS_2010;
end
results.agg_CE_SS1 = sum( sum(CE_1 .* Lambda_SS_2010) );
results.agg_earn_SS1 = sum(mean_e_2010 .* sum(Lambda_SS_2010, 2));
results.agg_wealth_SS1 = sum( sum(a_bar_SS_2010 .* Lambda_SS_2010) );
results.Xi_trans = Xi_trans;
results.Lambda_trans = Lambda_trans;
results.agg_earn = agg_earn;
results.agg_wealth = agg_wealth;
results.agg_cons = agg_cons;
results.agg_mobility = agg_mobility;
results.WY_t = WY_t;
results.mean_earn_t = mean_earn_t;
results.mean_wealth_t = mean_wealth_t;
results.mean_ce_t = mean_ce_t;
results.paths = paths;
results.output_dir = output_dir;
results.T = T;
results.TRANSITION_TYPE = TRANSITION_TYPE;

% Save variables for global map scaling
results.map_P_2010 = P_2010;
results.map_P = P_2100;
results.map_Lambda1 = Lambda_SS_2010;
results.map_Lambda2 = Lambda_SS_2100;
results.map_W1 = mean_w_2010;
results.map_W2 = mean_w_2100;
results.map_wye1 = wye_2010;
results.map_wye2 = wye_2100;
results.map_mean_w_2010 = mean_w_2010;
results.map_delta_CE_mean = delta_CE_mean;
results.map_delta_CE_var = delta_CE_var;
results.map_delta_W_var = delta_W_var;

save(fullfile(output_dir, 'dynamic_results.mat'), '-struct', 'results');

end
