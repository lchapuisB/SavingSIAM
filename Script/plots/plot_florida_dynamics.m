function plot_florida_dynamics(P, suffix)
    set_aesthetics('subplot');
% plot_florida_dynamics.m
if nargin < 2
    suffix = '';
end

close all;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'scenarios'), fullfile(script_dir, 'static'), fullfile(script_dir, 'helpers'), fullfile(script_dir, 'plots'), fullfile(script_dir, 'dynamics'), fullfile(script_dir, 'no_saving'));

idx_fl = find(strcmpi(P.state_names, 'Florida'));

%% 2. Load Shocks
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
shock_amen = damages_amen.shock_amenity;
shock_wage = damages_wage.shock_wage;

%% 3. Run Scenarios
fprintf('Running Normal Psi...\n');
P_normal = P;
Lambda_trans_normal = sim_scenario(P_normal, shock_wage, shock_housing, shock_amen);

fprintf('Running Psi -90%%...\n');
P_low = P;
P_low.psi = P.psi * 0.9;
Lambda_trans_low = sim_scenario(P_low, shock_wage, shock_housing, shock_amen);

fprintf('Running Psi +10%%...\n');
P_high = P;
P_high.psi = P.psi * 1.1;
Lambda_trans_high = sim_scenario(P_high, shock_wage, shock_housing, shock_amen);

%% 4. Extract and Plot Florida
pop_normal = squeeze(sum(Lambda_trans_normal(idx_fl, :, :), 2));
pop_low    = squeeze(sum(Lambda_trans_low(idx_fl, :, :), 2));
pop_high   = squeeze(sum(Lambda_trans_high(idx_fl, :, :), 2));

T = length(pop_normal);
t_vec = 1:T;

% Variation from SS1
pct_normal = (pop_normal - pop_normal(1)) / pop_normal(1) * 100;
pct_low    = (pop_low - pop_low(1)) / pop_low(1) * 100;
pct_high   = (pop_high - pop_high(1)) / pop_high(1) * 100;

fig = figure('Visible', 'off', 'Position', [100 100 1200 500]);
subplot(1,2,1);
set(gcf, 'Tag', 'population_share');
plot(t_vec, pop_normal, 'k-', 'LineWidth', 3.0, 'DisplayName', 'Normal Cost'); hold on;
plot(t_vec, pop_low, 'LineStyle', '--', 'LineWidth', 3.0, 'DisplayName', 'Cost -10%');
plot(t_vec, pop_high, 'LineStyle', ':', 'LineWidth', 3.0, 'DisplayName', 'Cost +10%');
xline(90, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off');
xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
xlabel('Year'); ylabel('Population Share of Florida');
legend('Location', 'best'); grid on;
yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

subplot(1,2,2);
set(gcf, 'Tag', 'population_pct_change');
plot(t_vec, pct_normal, 'k-', 'LineWidth', 3.0, 'DisplayName', 'Normal Cost'); hold on;
plot(t_vec, pct_low, 'LineStyle', '--', 'LineWidth', 3.0, 'DisplayName', 'Cost -10%');
plot(t_vec, pct_high, 'LineStyle', ':', 'LineWidth', 3.0, 'DisplayName', 'Cost +10%');
xline(90, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off');
yline(0, ':k', 'HandleVisibility', 'off');
xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
xlabel('Year'); ylabel('\Delta Population (% Change from SS1)');
legend('Location', 'best'); grid on;
yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

% sgtitle('Florida Population Dynamics (Perfect Foresight)', 'FontSize', 20, 'FontWeight', 'bold');

output_dir = fullfile(script_dir, '..', '..', 'Output', ['graph_perfect_foresight', suffix]);
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
save_fig(fig, output_dir, 'florida_dynamics_comparison.png');
fprintf('\nDone! Split graphs saved to %s\n', output_dir);
end


%% Helper Function
function Lambda_trans = sim_scenario(P_local, shock_wage, shock_housing, shock_amen)
    % SS1
    Xi_SS1 = compute_Xi(P_local);
    [~, ~, ~, ~, ~, ~, Lambda_SS1, ~] = compute_SS_moments(P_local, Xi_SS1);
    
    % SS2 (Terminal)
    P_2100 = P_local;
    P_2100.w = P_local.w .* (1 + shock_wage);
    P_2100.q = P_local.q .* (1 + shock_housing);
    P_2100.A = P_local.A - shock_amen .* P_local.w_tilde;
    P_2100.w_tilde = P_local.w_tilde .* (1 + shock_wage);
    P_2100.y_tilde = P_2100.w_tilde - P_2100.q;

    P_trans_start = P_local;
    if isfield(P_local, 'uncertainty')
        if P_local.uncertainty == 1
            P_trans_start.sig2_eta = P_local.sig2_eta .* (1 + 50 * abs(shock_wage));
            P_2100.sig2_eta = P_local.sig2_eta;
        elseif P_local.uncertainty == 2
            P_trans_start.sig2_eta = P_local.sig2_eta;
            P_2100.sig2_eta = P_local.sig2_eta .* (1 + 50 * abs(shock_wage));
        elseif P_local.uncertainty == 0
            P_trans_start.sig2_eta = P_local.sig2_eta;
            P_2100.sig2_eta = P_local.sig2_eta;
        end
    end
    Xi_SS2 = compute_Xi(P_2100);
    
    % Transition Path
    T_shock = 90;  
    T_sim   = 500; 
    paths_temp = build_paths(P_trans_start, P_2100, T_shock, 'linear');
    keys = {'w','w_tilde','y_tilde','sig2_eps','sig2_eta', 'A', 'q'};
    paths = struct(); 
    for k = 1:numel(keys)
        terminal_state = paths_temp.(keys{k})(:, end);
        steady_state_padding = repmat(terminal_state, 1, T_sim - T_shock);
        paths.(keys{k}) = [paths_temp.(keys{k}), steady_state_padding];
    end
    
    % Backward induction
    Xi_trans = solve_transition_backward(P_2100, Xi_SS2, paths, T_sim);
    
    % Forward simulation
    Lambda_trans = simulate_Lambda_forward(P_local, Lambda_SS1, Xi_trans, paths, T_sim);
end
