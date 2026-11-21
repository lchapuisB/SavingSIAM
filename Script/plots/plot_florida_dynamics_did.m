function plot_florida_dynamics_did(P_base, suffix)
    set_aesthetics('subplot');
% plot_florida_dynamics_did.m
if nargin < 2
    suffix = '';
end

close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'scenarios'), fullfile(script_dir, '..', 'static'), fullfile(script_dir, '..','helpers'), fullfile(script_dir, '..','plots'), fullfile(script_dir, '..','dynamics'), fullfile(script_dir, 'no_saving'));

idx_fl = find(strcmpi(P_base.state_names, 'Florida'));

%% 2. Compute TRUE Baseline SS1 (with normal psi)
fprintf('Computing True Baseline SS1...\n');
Xi_SS1_base = compute_Xi(P_base);
[~, ~, ~, ~, ~, ~, Lambda_SS1_base, ~] = compute_SS_moments(P_base, Xi_SS1_base);
pop_FL_base = sum(Lambda_SS1_base(idx_fl, :));

%% 3. Load Shocks
wage_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_wage_2080.csv');
opts = detectImportOptions(wage_damages_path);
damages_wage = readtable(wage_damages_path, opts);
amen_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_amenity_2080.csv');
opts = detectImportOptions(amen_damages_path);
damages_amen = readtable(amen_damages_path, opts);
housing_damages_path = fullfile(script_dir,'..', '..', 'Output', 'data', 'shock_housing_2080.csv');
opts = detectImportOptions(housing_damages_path);
damages_housing = readtable(housing_damages_path, opts);

shock_housing = damages_housing.shock_housing;
shock_amen = damages_amen.shock_amenity;
shock_wage = damages_wage.shock_wage;

%% 4. Run Scenarios (DiD)
fprintf('\n--- Running Scenario 1: Normal Psi ---\n');
P_normal = P_base;
[pop_treat_norm, pop_ctrl_norm] = sim_did_scenario(P_normal, Lambda_SS1_base, shock_wage, shock_housing, shock_amen, idx_fl);

fprintf('\n--- Running Scenario 2: Psi -90%% ---\n');
P_low = P_base;
P_low.psi = P_base.psi * 0.9;
[pop_treat_low, pop_ctrl_low] = sim_did_scenario(P_low, Lambda_SS1_base, shock_wage, shock_housing, shock_amen, idx_fl);

fprintf('\n--- Running Scenario 3: Psi +10%% ---\n');
P_high = P_base;
P_high.psi = P_base.psi * 1.1;
[pop_treat_high, pop_ctrl_high] = sim_did_scenario(P_high, Lambda_SS1_base, shock_wage, shock_housing, shock_amen, idx_fl);


%% 5. Plot the DiD Results
T = length(pop_treat_norm);
t_vec = 1:T;

% Difference in Levels (Treat - Ctrl)
diff_norm = pop_treat_norm - pop_ctrl_norm;
diff_low  = pop_treat_low  - pop_ctrl_low;
diff_high = pop_treat_high - pop_ctrl_high;

% Percentage Difference relative to the Ctrl population at each t
pct_norm = diff_norm ./ pop_ctrl_norm * 100;
pct_low  = diff_low  ./ pop_ctrl_low * 100;
pct_high = diff_high ./ pop_ctrl_high * 100;

fig = figure('Visible', 'off', 'Position', [100 100 1200 500]);

subplot(1,2,1);
set(gcf, 'Tag', 'population_share');
plot(t_vec, diff_norm, 'k-', 'LineWidth', 3.0, 'DisplayName', 'Normal Cost'); hold on;
plot(t_vec, diff_low, 'LineStyle', '--', 'LineWidth', 3.0, 'DisplayName', 'Cost -10%');
plot(t_vec, diff_high, 'LineStyle', ':', 'LineWidth', 3.0, 'DisplayName', 'Cost +10%');
xline(90, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off');
yline(0, ':k', 'HandleVisibility', 'off');
xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
xlabel('Year'); ylabel('\Delta Population Share (Treat - Ctrl)');
legend('Location', 'best'); grid on;
yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

subplot(1,2,2);
set(gcf, 'Tag', 'population_pct_change');
plot(t_vec, pct_norm, 'k-', 'LineWidth', 3.0, 'DisplayName', 'Normal Cost'); hold on;
plot(t_vec, pct_low, 'LineStyle', '--', 'LineWidth', 3.0, 'DisplayName', 'Cost -10%');
plot(t_vec, pct_high, 'LineStyle', ':', 'LineWidth', 3.0, 'DisplayName', 'Cost +10%');
xline(90, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off');
yline(0, ':k', 'HandleVisibility', 'off');
xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
xlabel('Year'); ylabel('DiD Effect (%)');
legend('Location', 'best'); grid on;
yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

% sg% title('Florida Population Dynamics (DiD Isolation of Climate Shock)', 'FontSize', 20, 'FontWeight', 'bold');

output_dir = fullfile(script_dir, '..', 'Output', ['graph_perfect_foresight', suffix]);
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
save_fig(fig, output_dir, 'florida_did_comparison.png');
fprintf('\nDone! Graph saved to %s\n', fullfile(output_dir, 'florida_did_comparison.png'));
end


%% Helper Function
function [pop_treat, pop_ctrl] = sim_did_scenario(P_local, Lambda_SS1_base, shock_wage, shock_housing, shock_amen, idx_fl)
    
    T_shock = 90;  
    T_sim   = 500; 

    % --- CONTROL ARM (New Psi, No Climate Shock) ---
    P_2100_ctrl = P_local;
    Xi_SS2_ctrl = compute_Xi(P_2100_ctrl);
    
    paths_temp_ctrl = build_paths(P_local, P_2100_ctrl, T_shock, 'linear');
    keys = {'w','w_tilde','y_tilde','sig2_eps','sig2_eta', 'A', 'q'};
    paths_ctrl = struct(); 
    for k = 1:numel(keys)
        terminal_state = paths_temp_ctrl.(keys{k})(:, end);
        steady_state_padding = repmat(terminal_state, 1, T_sim - T_shock);
        paths_ctrl.(keys{k}) = [paths_temp_ctrl.(keys{k}), steady_state_padding];
    end
    
    Xi_trans_ctrl = solve_transition_backward(P_2100_ctrl, Xi_SS2_ctrl, paths_ctrl, T_sim);
    Lambda_trans_ctrl = simulate_Lambda_forward(P_local, Lambda_SS1_base, Xi_trans_ctrl, paths_ctrl, T_sim);
    
    pop_ctrl = squeeze(sum(Lambda_trans_ctrl(idx_fl, :, :), 2))';

    % --- TREATMENT ARM (New Psi, WITH Climate Shock) ---
    P_2100_treat = P_local;
    P_2100_treat.w = P_local.w .* (1 + shock_wage);
    P_2100_treat.q = P_local.q .* (1 + shock_housing);
    P_2100_treat.A = P_local.A - shock_amen .* P_local.w_tilde;
    P_2100_treat.w_tilde = P_local.w_tilde .* (1 + shock_wage);
    P_2100_treat.y_tilde = P_2100_treat.w_tilde - P_2100_treat.q;

    P_trans_start_treat = P_local;
    if isfield(P_local, 'uncertainty')
        if P_local.uncertainty == 1
            P_trans_start_treat.sig2_eta = P_local.sig2_eta .* (1 + 50 * abs(shock_wage));
            P_2100_treat.sig2_eta = P_local.sig2_eta;
        elseif P_local.uncertainty == 2
            P_trans_start_treat.sig2_eta = P_local.sig2_eta;
            P_2100_treat.sig2_eta = P_local.sig2_eta .* (1 + 50 * abs(shock_wage));
        elseif P_local.uncertainty == 0
            P_trans_start_treat.sig2_eta = P_local.sig2_eta;
            P_2100_treat.sig2_eta = P_local.sig2_eta;
        end
    end
    Xi_SS2_treat = compute_Xi(P_2100_treat);
    
    paths_temp_treat = build_paths(P_trans_start_treat, P_2100_treat, T_shock, 'linear');
    paths_treat = struct(); 
    for k = 1:numel(keys)
        terminal_state = paths_temp_treat.(keys{k})(:, end);
        steady_state_padding = repmat(terminal_state, 1, T_sim - T_shock);
        paths_treat.(keys{k}) = [paths_temp_treat.(keys{k}), steady_state_padding];
    end
    
    Xi_trans_treat = solve_transition_backward(P_2100_treat, Xi_SS2_treat, paths_treat, T_sim);
    Lambda_trans_treat = simulate_Lambda_forward(P_local, Lambda_SS1_base, Xi_trans_treat, paths_treat, T_sim);
    
    pop_treat = squeeze(sum(Lambda_trans_treat(idx_fl, :, :), 2))';
end
