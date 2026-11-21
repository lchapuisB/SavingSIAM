% ============================================================
% Main Orchestrator for Climate Shock Scenarios
% ============================================================
%% -- Initializing baseline parameters
clear all; clc; close all force; %#ok<CLALL>
clear global my_sub_figs;
set(0, 'DefaultFigureVisible', 'off');

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'scenarios'), fullfile(script_dir, 'static'), fullfile(script_dir, 'helpers'), fullfile(script_dir, 'plots'), fullfile(script_dir, 'dynamics'), fullfile(script_dir, 'no_saving'));

model_static;
RUN_SCENARIO = false; 
P.uncertainty = 0;

%% -- Run Scenarios
if RUN_SCENARIO
    %% -- SCENARIO A: PERFECT FORESIGHT
    res_A = run_scenario('perfect_foresight', P, 'no_counterfactual'); %#ok<*UNRCH>
    %% -- SCENARIO B: MYOPIA (STATIC EXPECTATIONS)
    res_B = run_scenario('myopia', P, 'no_counterfactual');

    %% --  COUNTERFACTUAL SCENARIO: NO MIGRATION (DID)
    res_C = run_did(P, 'no_migration');

    %% --  COUNTERFACTUAL SCENARIO: LOW COST OF MIGRATION (DID)
    res_D = run_did(P, 'big_migration');

    %% --  COUNTERFACTUAL SCENARIO: HIGHCOST OF MIGRATION (DID)
    res_E = run_did(P, 'small_migration');

    %% --  COUNTERFACTUAL SCENARIO: NO SAVING (STRUCTURAL)
    res_F = run_scenario('perfect_foresight', P, 'no_saving_struct');

    %% -- UNCERTAINTY SCENARIO 1: DAMAGE PATH UNCERTAINTY (res_G)
    P_unc1 = P; P_unc1.uncertainty = 1;
    res_G = run_scenario('perfect_foresight', P_unc1, 'no_counterfactual');

    %% -- UNCERTAINTY SCENARIO 2: VARIANCE SHOCKS (res_H)
    P_unc2 = P; P_unc2.uncertainty = 2;
    res_H = run_scenario('perfect_foresight', P_unc2, 'no_counterfactual');

    %% -- UNCERTAINTY SCENARIO 3: VARIANCE SHOCKS + NO MIGRATION (res_I)
    P_unc3 = P; P_unc3.uncertainty = 2;
    res_I = run_did(P_unc3, 'no_migration');
else
    %% -- Load pre-computed results from Output directory
    fprintf('Loading pre-computed results from Output directory...\n');
    output_base = fullfile(fileparts(script_dir),'..', 'Output', 'big_nu');
    
    res_A = load(fullfile(output_base, 'graph_perfect_foresight', 'dynamic_results.mat'));
    res_B = load(fullfile(output_base, 'graph_myopia', 'dynamic_results.mat'));
    res_C = load(fullfile(output_base, 'graph_no_migration_did', 'did_results.mat'));
    res_D = load(fullfile(output_base, 'graph_big_migration_did', 'did_results.mat'));
    res_E = load(fullfile(output_base, 'graph_small_migration_did', 'did_results.mat'));
    res_F = load(fullfile(output_base, 'graph_no_saving_struct', 'dynamic_results.mat'));
    res_G = load(fullfile(output_base, 'graph_uncertainty_damage_path', 'dynamic_results.mat'));
    res_H = load(fullfile(output_base, 'graph_uncertainty_variance_shocks', 'dynamic_results.mat'));
    res_I = load(fullfile(output_base, 'graph_uncertainty_variance_shocks_no_migration_did', 'did_results.mat'));
end

%% -- GLOBAL MAP SYNCHRONIZATION
fprintf('\n--- Re-generating Maps with Global Bounds ---\n');
scenarios = {res_A, res_B, res_C, res_D, res_E, res_F, res_G, res_H, res_I};

all_W_pct = []; all_pop_pct = []; all_wy = []; all_ce = [];
all_ce_var = []; all_w_var = [];
all_dw = []; all_dA = []; all_dq = [];
all_W_peak = []; all_pop_peak = [];

for i = 1:length(scenarios)
    res = scenarios{i};
    if ~isfield(res, 'map_Lambda1')
        continue;
    end
    pop_SS1 = sum(res.map_Lambda1, 2);
    pop_SS2 = sum(res.map_Lambda2, 2);
    denom_W = abs(res.map_W2) + abs(res.map_W1);
    W_pct = zeros(size(res.map_W1));
    valid_W = denom_W > 1e-8;
    W_pct(valid_W) = 2 * (res.map_W2(valid_W) - res.map_W1(valid_W)) ./ denom_W(valid_W) * 100;
    
    denom_pop = abs(pop_SS2) + abs(pop_SS1);
    pop_pct = zeros(size(pop_SS1));
    valid_pop = denom_pop > 1e-8;
    pop_pct(valid_pop) = 2 * (pop_SS2(valid_pop) - pop_SS1(valid_pop)) ./ denom_pop(valid_pop) * 100;
    wy = res.map_wye2 - res.map_wye1;
    ce = res.map_delta_CE_mean;
    
    all_W_pct = [all_W_pct; W_pct(:)];
    all_pop_pct = [all_pop_pct; pop_pct(:)];
    all_wy = [all_wy; wy(:)];
    all_ce = [all_ce; ce(:)];

    % Inequality Maps
    if isfield(res, 'map_delta_CE_var')
        all_ce_var = [all_ce_var; res.map_delta_CE_var(:)];
        all_w_var = [all_w_var; res.map_delta_W_var(:)];
    end

    % Damage Maps
    if isfield(res, 'map_P_2010') && isfield(res, 'map_P')
        dw = (res.map_P.w - res.map_P_2010.w) ./ res.map_P_2010.w * 100;
        dA = (res.map_P.A - res.map_P_2010.A) ./ res.map_P_2010.gamma * 100;
        dq = (res.map_P.q - res.map_P_2010.q) ./ res.map_P_2010.q * 100;
        all_dw = [all_dw; dw(:)]; %#ok<*AGROW>
        all_dA = [all_dA; dA(:)];
        all_dq = [all_dq; dq(:)];
    end

    % Peak Transition Maps (run_scenario only)
    if isfield(res, 'WY_t') && isfield(res, 'mean_wealth_t') && isfield(res, 'map_mean_w_2010')
        [max_wy, t_peak] = max(res.WY_t);
        delta_W_peak = res.mean_wealth_t(:, t_peak) - res.map_mean_w_2010;
        
        if ndims(res.Lambda_trans) == 3
            pop_peak = squeeze(sum(res.Lambda_trans(:, :, t_peak), 2));
        else
            pop_peak = res.Lambda_trans(:, t_peak);
        end
        delta_pop_peak = pop_peak - pop_SS1;

        all_W_peak = [all_W_peak; delta_W_peak(:)];
        all_pop_peak = [all_pop_peak; delta_pop_peak(:)];
    end
end

clean_data = @(x) x(~isnan(x) & ~isinf(x));

all_W_pct = clean_data(all_W_pct);
all_pop_pct = clean_data(all_pop_pct);
all_wy = clean_data(all_wy);
all_ce = clean_data(all_ce);
all_ce_var = clean_data(all_ce_var);
all_w_var = clean_data(all_w_var);
all_dw = clean_data(all_dw);
all_dA = clean_data(all_dA);
all_dq = clean_data(all_dq);
all_W_peak = clean_data(all_W_peak);
all_pop_peak = clean_data(all_pop_peak);

bounds.W_pct = [prctile(all_W_pct, 2), prctile(all_W_pct, 98)];
bounds.pop_pct = [prctile(all_pop_pct, 2), prctile(all_pop_pct, 98)];
bounds.wy = [prctile(all_wy, 0), prctile(all_wy, 98)];
bounds.ce = [prctile(all_ce, 0), prctile(all_ce, 100)];

bounds_ineq = struct();
if ~isempty(all_ce_var)
    bounds_ineq.ce_var = [prctile(all_ce_var, 1), prctile(all_ce_var, 99)];
    bounds_ineq.w_var = [prctile(all_w_var, 2), prctile(all_w_var, 98)];
end

bounds_dmg = struct();
if ~isempty(all_dw)
    bounds_dmg.w = [prctile(all_dw, 0), prctile(all_dw, 100)];
    bounds_dmg.A = [prctile(all_dA, 0), prctile(all_dA, 100)];
    bounds_dmg.q = [prctile(all_dq, 0), prctile(all_dq, 100)];
end

bounds_peak = struct();
if ~isempty(all_W_peak)
    bounds_peak.w = [prctile(all_W_peak, 2), prctile(all_W_peak, 98)];
    bounds_peak.pop = [prctile(all_pop_peak, 2), prctile(all_pop_peak, 98)];
end

for i = 1:length(scenarios)
    res = scenarios{i};
    if ~isfield(res, 'map_Lambda1')
        continue;
    end
    
    curr_bounds = bounds;
    curr_bounds_ineq = bounds_ineq;
    curr_bounds_dmg = bounds_dmg;
    curr_bounds_peak = bounds_peak;
    
    if contains(res.output_dir, 'big_migration')
        curr_bounds = [];
        curr_bounds_ineq = [];
        curr_bounds_dmg = [];
        curr_bounds_peak = [];
    end

    plot_transition_maps(res.map_P, res.map_Lambda1, res.map_Lambda2, res.map_W1, res.map_W2, res.map_wye1, res.map_wye2, res.map_delta_CE_mean, res.output_dir, curr_bounds);
    
    if isfield(res, 'map_delta_CE_var')
        plot_inequality_maps(res.map_P.state_names, res.map_delta_CE_var, res.map_delta_W_var, res.output_dir, curr_bounds_ineq);
    end

    if isfield(res, 'map_P_2010') && isfield(res, 'map_P')
        plot_damage_maps(res.map_P_2010, res.map_P, res.output_dir, curr_bounds_dmg);
    end

    if isfield(res, 'WY_t') && isfield(res, 'mean_wealth_t') && isfield(res, 'map_mean_w_2010') && isfield(res, 'map_P_2010')
        plot_peak_transition_maps(res.map_P_2010, res.WY_t, res.mean_wealth_t, res.Lambda_trans, res.map_mean_w_2010, res.map_Lambda1, res.output_dir, curr_bounds_peak);
    end
end

%% -- other GRAPHS
fprintf('\n--- Generating Other Graphs  ---\n');
plot_ushape_t500_split(P, '');
plot_florida_dynamics(P, '');
plot_florida_dynamics_did(P, '');
%% -- CROSS-SCENARIO AUDIT
fprintf('\n======================================\n');
fprintf('Running Cross-Scenario Audit Analysis...\n');
compare_all_audit(res_A, res_B, res_C, res_D, res_E, res_F, res_G, res_H, res_I);
