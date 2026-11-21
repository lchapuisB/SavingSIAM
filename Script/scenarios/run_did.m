function results = run_did(P, counterfactual)
% ============================================================
% No-Migration Counterfactual via Difference-in-Differences
%
% Both arms start from the calibrated 2010 SS (with normal migration).
%   Transition 1 (control) : migration ban at t=1, NO climate change
%   Transition 2 (treatment): migration ban at t=1, WITH climate change
%   Diff (1 - 2) = pure CC effect when migration is unavailable
% ============================================================
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'plots'), ...
        fullfile(script_dir, '..', 'dynamics'), ...
        fullfile(script_dir, '..', 'helpers'), ...
        fullfile(script_dir, '..', 'static'));

if P.uncertainty == 1 & strcmp(counterfactual, 'no_migration')
    output_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_uncertainty_damage_path_no_migration_did');
elseif P.uncertainty == 2 & strcmp(counterfactual, 'no_migration')
    output_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_uncertainty_variance_shocks_no_migration_did');
elseif strcmp(counterfactual, 'no_migration')
    output_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_no_migration_did');
elseif strcmp(counterfactual, 'big_migration')
    output_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_big_migration_did');
elseif strcmp(counterfactual, 'small_migration')
    output_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_small_migration_did');
else
    error('Invalid counterfactual option. Choose from: no_migration, big_migration, small_migration.');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ── STEP 1 : BASELINE SS (shared, normal migration) ─────────────────────
fprintf('\n============================================================\n');
fprintf('NO-MIGRATION DID COUNTERFACTUAL\n');
fprintf('============================================================\n');
fprintf('\n=== Step 1: Computing Baseline SS (normal migration) ===\n');

P_base = P;  % calibrated P with normal psi

Xi_SS_2010 = compute_Xi(P_base);
[~, ~, ~, mean_w_2010, mean_e_2010, ~, Lambda_SS_2010, a_bar_SS_2010] = ...
    compute_SS_moments(P_base, Xi_SS_2010);

fprintf('  Baseline W/Y = %.3f\n', ...
    sum(mean_w_2010 .* sum(Lambda_SS_2010, 2)) / ...
    sum(mean_e_2010 .* sum(Lambda_SS_2010, 2)));

%% ── STEP 2 : BUILD MIGRATION BAN ────────────────────────────────────────
fprintf('\n=== Step 2: Building migration ban psi matrix ===\n');
if strcmp(counterfactual, 'no_migration')
    psi_new = 1e8 * ones(P.J, P.J);
    psi_new(logical(eye(P.J))) = 0;
elseif strcmp(counterfactual, 'big_migration')
    psi_new = P.psi * 0.9;
elseif strcmp(counterfactual, 'small_migration')
    psi_new = P.psi * 1.1;
end

%% ── STEP 3 : LOAD CC SHOCKS (needed for the treatment arm) ──────────────
fprintf('\n=== Step 3: Loading climate shocks ===\n');

wage_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_wage_2080.csv');
opts = detectImportOptions(wage_damages_path);
damages_wage = readtable(wage_damages_path, opts);

amen_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_amenity_2080.csv');
opts = detectImportOptions(amen_damages_path);
damages_amen = readtable(amen_damages_path, opts);

housing_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_housing_2080.csv');
opts = detectImportOptions(housing_damages_path);
damages_housing = readtable(housing_damages_path, opts);

shock_housing = damages_housing.shock_housing;
shock_amen    = damages_amen.shock_amenity;
shock_wage    = damages_wage.shock_wage;

%% ── STEP 4 : TRANSITION SETTINGS ────────────────────────────────────────
TRANSITION_TYPE = 'linear';
T_shock = 90;
T_sim   = 500;

%% ── STEP 5 : CONTROL ARM vhange in migration cost, NO CC) ─────────────────────────
fprintf('\n=== Step 5: Control arm — migration ban, no CC ===\n');

P_ctrl = P_base;
P_ctrl.psi = psi_new;


P_2100_ctrl = P_ctrl;

fprintf('  Computing terminal SS (autarky, no CC)...\n');
Xi_SS_term_ctrl = compute_Xi(P_2100_ctrl);

paths_ctrl_temp = build_paths(P_ctrl, P_2100_ctrl, T_shock, TRANSITION_TYPE);

keys = {'w', 'w_tilde', 'y_tilde', 'sig2_eps', 'sig2_eta', 'A', 'q'};
paths_ctrl = struct();
for k = 1:numel(keys)
    terminal_state = paths_ctrl_temp.(keys{k})(:, end);
    steady_state_padding = repmat(terminal_state, 1, T_sim - T_shock);
    paths_ctrl.(keys{k}) = [paths_ctrl_temp.(keys{k}), steady_state_padding];
end

T = T_sim;

fprintf('  Backward induction (control)...\n');
Xi_trans_ctrl = solve_transition_backward(P_ctrl, Xi_SS_term_ctrl, paths_ctrl, T);

fprintf('  Forward simulation (control)...\n');
Lambda_trans_ctrl = simulate_Lambda_forward(P_ctrl, Lambda_SS_2010, Xi_trans_ctrl, paths_ctrl, T);

fprintf('  Computing moments (control)...\n');
[agg_earn_ctrl, agg_wealth_ctrl, ~, agg_mobility_ctrl, mean_earn_ctrl, ...
 mean_wealth_ctrl, WY_ctrl, mean_ce_ctrl, ~] = ...
    compute_moments_dyn(Lambda_trans_ctrl, Xi_trans_ctrl, a_bar_SS_2010, paths_ctrl, P_ctrl);

%% ── STEP 6 : TREATMENT ARM (migration ban, WITH CC) ─────────────────────
fprintf('\n=== Step 6: Treatment arm — migration ban, with CC ===\n');

P_treat = P_base;
P_treat.psi = psi_new;

P_2100_treat = P_treat;
P_2100_treat.w = P_treat.w .* (1 + shock_wage);
P_2100_treat.q = P_treat.q .* (1 + shock_housing);
P_2100_treat.A = P_treat.A - shock_amen;
P_2100_treat.w_tilde = P_treat.w_tilde .* (1 + shock_wage);
P_2100_treat.y_tilde = P_2100_treat.w_tilde - P_2100_treat.q;


P_trans_start_treat = P_treat;
if isfield(P_treat, 'uncertainty')
    if P_treat.uncertainty == 1
        P_trans_start_treat.sig2_eta = P_treat.sig2_eta .* (1 + 50 * abs(shock_wage));
        P_2100_treat.sig2_eta = P_treat.sig2_eta;
    elseif P_treat.uncertainty == 2
        P_trans_start_treat.sig2_eta = P_treat.sig2_eta;
        P_2100_treat.sig2_eta = P_treat.sig2_eta .* (1 + 50 * abs(shock_wage));
    end
end


fprintf('  Computing terminal SS (autarky, with CC)...\n');
Xi_SS_term_treat = compute_Xi(P_2100_treat);


paths_treat_temp = build_paths(P_trans_start_treat, P_2100_treat, T_shock, TRANSITION_TYPE);

paths_treat = struct();
for k = 1:numel(keys)
    terminal_state = paths_treat_temp.(keys{k})(:, end);
    steady_state_padding = repmat(terminal_state, 1, T_sim - T_shock);
    paths_treat.(keys{k}) = [paths_treat_temp.(keys{k}), steady_state_padding];
end

% Backward induction
fprintf('  Backward induction (treatment)...\n');
Xi_trans_treat = solve_transition_backward(P_treat, Xi_SS_term_treat, paths_treat, T);

% Forward simulation
fprintf('  Forward simulation (treatment)...\n');
Lambda_trans_treat = simulate_Lambda_forward(P_treat, Lambda_SS_2010, Xi_trans_treat, paths_treat, T);

% Dynamic moments
fprintf('  Computing moments (treatment)...\n');
[agg_earn_treat, agg_wealth_treat, ~, agg_mobility_treat, mean_earn_treat, ...
 mean_wealth_treat, WY_treat, mean_ce_treat, ~] = ...
    compute_moments_dyn(Lambda_trans_treat, Xi_trans_treat, a_bar_SS_2010, paths_treat, P_treat);

%% ── STEP 7 : DIFFERENCE-IN-DIFFERENCES ──────────────────────────────────
fprintf('\n=== Step 7: Computing DID (Control - Treatment) ===\n');

% delta_CE > 0 means CC caused welfare loss (control is higher than treatment)
delta_ce_t = mean_ce_ctrl - mean_ce_treat;  

% Population-weighted aggregate CE at each t
pop_ctrl = squeeze(sum(Lambda_trans_ctrl, 2));   
pop_treat = squeeze(sum(Lambda_trans_treat, 2)); 

agg_CE_ctrl  = sum(mean_ce_ctrl .* pop_ctrl, 1);
agg_CE_treat = sum(mean_ce_treat .* pop_treat, 1);

fprintf('  Aggregate CE (control,  t=1): %.4f\n', agg_CE_ctrl(1));
fprintf('  Aggregate CE (treatment, t=1): %.4f\n', agg_CE_treat(1));
fprintf('  Aggregate CE (control,  t=%d): %.4f\n', T_shock, agg_CE_ctrl(T_shock));
fprintf('  Aggregate CE (treatment, t=%d): %.4f\n', T_shock, agg_CE_treat(T_shock));
fprintf('  Agg welfare loss from CC at t=%d: %.4f%%\n', T_shock, ...
    (agg_CE_treat(T_shock) - agg_CE_ctrl(T_shock)) / abs(agg_CE_ctrl(T_shock)) * 100);

%% ── STEP 8 : PLOTS ──────────────────────────────────────────────────────
fprintf('\n=== Step 8: Plotting DID results ===\n');

fprintf('  Computing explicit terminal SS moments for maps...\n');
[~, ~, ~, mean_wealth_term_ctrl, mean_e_term_ctrl, wye_ctrl, Lambda_SS_term_ctrl, a_bar_SS_term_ctrl] = ...
    compute_SS_moments(P_2100_ctrl, Xi_SS_term_ctrl);
[~, ~, ~, mean_wealth_term_treat, mean_e_term_treat, wye_treat, Lambda_SS_term_treat, a_bar_SS_term_treat] = ...
    compute_SS_moments(P_2100_treat, Xi_SS_term_treat);

% Counterfactual wealth for structural decomposition (Treat wages/amenities, Ctrl expectations)
delta_w = P_2100_treat.w_tilde(:) - P_2100_ctrl.w_tilde(:);
kappa_stat = P_2100_ctrl.r / (1 + P_2100_ctrl.r);
Xi_cf_wage_did = Xi_SS_term_ctrl + kappa_stat * (delta_w * P_2100_treat.z_grid(:)');

[~, ~, ~, ~, ~, ~, ~, a_bar_cf_did] = ...
    compute_SS_moments(P_2100_treat, Xi_cf_wage_did);
    
Xi_cf_local_did = compute_Xi_cf_local(P_2100_ctrl, P_2100_treat, Xi_SS_term_ctrl);
[~, ~, ~, ~, ~, ~, ~, a_bar_cf_local_did] = ...
    compute_SS_moments(P_2100_treat, Xi_cf_local_did);

%  Spatial Decomposition (DiD: Treatment vs Control in 2100)
plot_spatial_decomposition(P_2100_ctrl, P_2100_treat, Xi_SS_term_ctrl, Xi_SS_term_treat, ...
    output_dir, sprintf('did_%s_Treat_vs_Ctrl', counterfactual));


% 1. Variance and Inequality Maps (Treatment vs Control)
[delta_CE_mean, delta_CE_var, delta_W_var] = plot_welfare_variance(P_2100_ctrl, P_2100_treat, ...
    Lambda_SS_term_ctrl, Lambda_SS_term_treat, ...
    a_bar_SS_term_ctrl, a_bar_SS_term_treat, ...
    Xi_SS_term_ctrl, Xi_SS_term_treat, output_dir);

% 3. Dynamic Moments (Level Diffs and Pure DiD)
plot_did_dynamic_moments(T, T_shock, agg_earn_ctrl, agg_earn_treat, ...
    agg_wealth_ctrl, agg_wealth_treat, WY_ctrl, WY_treat, ...
    agg_mobility_ctrl, agg_mobility_treat, output_dir);

plot_did_welfare_transition(T, T_shock, mean_earn_ctrl, mean_earn_treat, ...
    mean_wealth_ctrl, mean_wealth_treat, mean_ce_ctrl, mean_ce_treat, output_dir);

% 4. Welfare Comparisons and Diagnostics
plot_did_comparison(P_base, T, T_shock, mean_ce_ctrl, mean_ce_treat, delta_ce_t, ...
    agg_CE_ctrl, agg_CE_treat, Lambda_trans_ctrl, Lambda_trans_treat, ...
    mean_wealth_ctrl, mean_wealth_treat, agg_mobility_ctrl, agg_mobility_treat, output_dir);

fprintf('\n=== Step 9: Saving results ===\n');

results = struct();
kappa_base = P_base.r / (1 + P_base.r);
CE_1_base = kappa_base * a_bar_SS_2010 + Xi_SS_2010;
agg_CE_SS1 = sum( sum(CE_1_base .* Lambda_SS_2010) );
agg_earn_SS1 = sum(mean_e_2010 .* sum(Lambda_SS_2010, 2));
agg_wealth_SS1 = sum( sum(a_bar_SS_2010 .* Lambda_SS_2010) );

results.agg_CE_SS1_ctrl = agg_CE_SS1;
results.agg_CE_SS1_treat = agg_CE_SS1;
results.agg_earn_SS1_ctrl = agg_earn_SS1;
results.agg_earn_SS1_treat = agg_earn_SS1;
results.agg_wealth_SS1_ctrl = agg_wealth_SS1;
results.agg_wealth_SS1_treat = agg_wealth_SS1;
results.Xi_trans_ctrl     = Xi_trans_ctrl;
results.Xi_trans_treat    = Xi_trans_treat;
results.Lambda_trans_ctrl = Lambda_trans_ctrl;
results.Lambda_trans_treat = Lambda_trans_treat;
results.mean_ce_ctrl      = mean_ce_ctrl;
results.mean_ce_treat     = mean_ce_treat;
results.delta_ce_t        = delta_ce_t;
results.agg_CE_ctrl       = agg_CE_ctrl;
results.agg_CE_treat      = agg_CE_treat;
results.mean_wealth_ctrl  = mean_wealth_ctrl;
results.mean_wealth_treat = mean_wealth_treat;
results.WY_ctrl           = WY_ctrl;
results.WY_treat          = WY_treat;
results.agg_earn_ctrl     = agg_earn_ctrl;
results.agg_earn_treat    = agg_earn_treat;
results.agg_wealth_ctrl   = agg_wealth_ctrl;
results.agg_wealth_treat  = agg_wealth_treat;
results.agg_mobility_ctrl = agg_mobility_ctrl;
results.agg_mobility_treat = agg_mobility_treat;
results.output_dir        = output_dir;

% Save variables for global map scaling
results.map_P_2010 = P_2100_ctrl;
results.map_P = P_2100_treat;
results.map_Lambda1 = Lambda_SS_term_ctrl;
results.map_Lambda2 = Lambda_SS_term_treat;
results.map_W1 = mean_wealth_term_ctrl;
results.map_W2 = mean_wealth_term_treat;
results.map_wye1 = wye_ctrl;
results.map_wye2 = wye_treat;
results.map_delta_CE_mean = delta_CE_mean;
results.map_delta_CE_var = delta_CE_var;
results.map_delta_W_var = delta_W_var;

save(fullfile(output_dir, 'did_results.mat'), '-struct', 'results');
fprintf('  Results saved to %s\n', output_dir);
fprintf('\n============================================================\n');
fprintf('DID COUNTERFACTUAL COMPLETE\n');
fprintf('============================================================\n');

end
