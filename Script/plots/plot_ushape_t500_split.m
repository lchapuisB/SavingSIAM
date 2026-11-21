function plot_ushape_t500_split(P, suffix)
    set_aesthetics('subplot');
% plot_ushape_t500_split.m
if nargin < 2
    suffix = '';
end

close all;
set(0, 'DefaultFigureVisible', 'off');

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'static'), fullfile(script_dir, '..', 'helpers'), fullfile(script_dir, '..', 'dynamics'));

%% 1. BASELINE SETUP
fprintf('\n--- Running Dynamic U-Shape Sweep at t=500 (Split Graphs with Decomposition) ---\n');
idx_FL = 10; % Florida
idx_WA = 48; % Washington
idx_CA = 5;  % California
J = P.J;

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

%% 2. SWEEP
psi_mults = [0.01, 0.1:0.1:4.0]; 
N_sweep = length(psi_mults);

fl_drops_dyn = zeros(N_sweep, 1);
fl_rel_stayers = zeros(N_sweep, 1);
fl_rel_inflow = zeros(N_sweep, 1);

wa_drops_dyn = zeros(N_sweep, 1);
wa_rel_stayers = zeros(N_sweep, 1);
wa_rel_inflow = zeros(N_sweep, 1);

ca_drops_dyn = zeros(N_sweep, 1);
ca_rel_stayers = zeros(N_sweep, 1);
ca_rel_inflow = zeros(N_sweep, 1);

T_shock = 90;
T_eval = 500;
keys = {'w', 'w_tilde', 'y_tilde', 'sig2_eps', 'sig2_eta', 'A', 'q'};

out_dir = fullfile(script_dir, '..', '..', 'Output', ['graph_ushape', suffix]);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

data_file = fullfile(out_dir, 'ushape_t500_data.mat');

if exist(data_file, 'file')
    load(data_file);
    fprintf('Loaded existing data from %s\n', data_file);
else
    for i = 1:N_sweep
        mult = psi_mults(i);
        fprintf('Computing Dynamic Transition for psi multiplier: %.2f\n', mult);
        
        % --- Control ---
        P_ctrl = P;
        P_ctrl.psi = P.psi * mult;
        Xi_SS_ctrl = compute_Xi(P_ctrl);
        [~, ~, ~, ~, ~, ~, Lam_SS_ctrl, ~] = compute_SS_moments(P_ctrl, Xi_SS_ctrl);
        
        paths_ctrl_temp = build_paths(P_ctrl, P_ctrl, T_shock, 'linear');
        paths_ctrl = struct();
        for k = 1:numel(keys)
            terminal_state = paths_ctrl_temp.(keys{k})(:, end);
            steady_state_padding = repmat(terminal_state, 1, T_eval - T_shock);
            paths_ctrl.(keys{k}) = [paths_ctrl_temp.(keys{k}), steady_state_padding];
        end
        
        Xi_trans_ctrl = solve_transition_backward(P_ctrl, Xi_SS_ctrl, paths_ctrl, T_eval);
        Lam_trans_ctrl = simulate_Lambda_forward(P_ctrl, Lam_SS_ctrl, Xi_trans_ctrl, paths_ctrl, T_eval);
        pop_ctrl_teval = sum(Lam_trans_ctrl(:,:,end), 2);
        
        % --- Treatment ---
        P_treat = P_ctrl;
        P_treat.w = P_treat.w .* (1 + shock_wage);
        P_treat.q = P_treat.q .* (1 + shock_housing);
        P_treat.A = P_treat.A - shock_amen .* P_treat.w_tilde;
        P_treat.w_tilde = P_treat.w_tilde .* (1 + shock_wage);
        P_treat.y_tilde = P_treat.w_tilde - P_treat.q;

        P_trans_start = P_ctrl;
        if isfield(P_treat, 'uncertainty')
            if P_treat.uncertainty == 1
                P_trans_start.sig2_eta = P_ctrl.sig2_eta .* (1 + 50 * abs(shock_wage));
                P_treat.sig2_eta = P_ctrl.sig2_eta;
            elseif P_treat.uncertainty == 2
                P_trans_start.sig2_eta = P_ctrl.sig2_eta;
                P_treat.sig2_eta = P_ctrl.sig2_eta .* (1 + 50 * abs(shock_wage));
            elseif P_treat.uncertainty == 0
                P_trans_start.sig2_eta = P_ctrl.sig2_eta;
                P_treat.sig2_eta = P_ctrl.sig2_eta;
            end
        end
        Xi_SS_treat = compute_Xi(P_treat);
        
        paths_treat_temp = build_paths(P_trans_start, P_treat, T_shock, 'linear');
        paths_treat = struct();
        for k = 1:numel(keys)
            terminal_state = paths_treat_temp.(keys{k})(:, end);
            steady_state_padding = repmat(terminal_state, 1, T_eval - T_shock);
            paths_treat.(keys{k}) = [paths_treat_temp.(keys{k}), steady_state_padding];
        end
        
        Xi_trans_treat = solve_transition_backward(P_ctrl, Xi_SS_treat, paths_treat, T_eval);
        Lam_trans_treat = simulate_Lambda_forward(P_ctrl, Lam_SS_ctrl, Xi_trans_treat, paths_treat, T_eval);
        pop_treat_teval = sum(Lam_trans_treat(:,:,end), 2);
        
        % --- Flows at T_eval ---
        Pi_all_ctrl = build_Pi_all(P_ctrl);
        B1hat_ctrl = compute_B1hat(P_ctrl);
        mu_ctrl = compute_mu_all(Xi_SS_ctrl, B1hat_ctrl, P_ctrl, Pi_all_ctrl);
        
        Flow_ctrl_t500 = zeros(J, J);
        for j = 1:J
            for l = 1:J
                Flow_ctrl_t500(j, l) = sum(Lam_trans_ctrl(j, :, end) .* squeeze(mu_ctrl(j, l, :))');
            end
        end
        
        Pi_all_treat = build_Pi_all(P_treat);
        B1hat_treat = compute_B1hat(P_treat);
        mu_treat = compute_mu_all(Xi_SS_treat, B1hat_treat, P_treat, Pi_all_treat);
        
        Flow_treat_t500 = zeros(J, J);
        for j = 1:J
            for l = 1:J
                Flow_treat_t500(j, l) = sum(Lam_trans_treat(j, :, end) .* squeeze(mu_treat(j, l, :))');
            end
        end

        % Drops & Decomposition (Relative %)
        fl_drops_dyn(i) = (pop_treat_teval(idx_FL) - pop_ctrl_teval(idx_FL)) / pop_ctrl_teval(idx_FL) * 100;
        fl_rel_stayers(i) = (Flow_treat_t500(idx_FL, idx_FL) - Flow_ctrl_t500(idx_FL, idx_FL)) / pop_ctrl_teval(idx_FL) * 100;
        fl_rel_inflow(i) = ((sum(Flow_treat_t500(:, idx_FL)) - Flow_treat_t500(idx_FL, idx_FL)) - (sum(Flow_ctrl_t500(:, idx_FL)) - Flow_ctrl_t500(idx_FL, idx_FL))) / pop_ctrl_teval(idx_FL) * 100;

        wa_drops_dyn(i) = (pop_treat_teval(idx_WA) - pop_ctrl_teval(idx_WA)) / pop_ctrl_teval(idx_WA) * 100;
        wa_rel_stayers(i) = (Flow_treat_t500(idx_WA, idx_WA) - Flow_ctrl_t500(idx_WA, idx_WA)) / pop_ctrl_teval(idx_WA) * 100;
        wa_rel_inflow(i) = ((sum(Flow_treat_t500(:, idx_WA)) - Flow_treat_t500(idx_WA, idx_WA)) - (sum(Flow_ctrl_t500(:, idx_WA)) - Flow_ctrl_t500(idx_WA, idx_WA))) / pop_ctrl_teval(idx_WA) * 100;

        ca_drops_dyn(i) = (pop_treat_teval(idx_CA) - pop_ctrl_teval(idx_CA)) / pop_ctrl_teval(idx_CA) * 100;
        ca_rel_stayers(i) = (Flow_treat_t500(idx_CA, idx_CA) - Flow_ctrl_t500(idx_CA, idx_CA)) / pop_ctrl_teval(idx_CA) * 100;
        ca_rel_inflow(i) = ((sum(Flow_treat_t500(:, idx_CA)) - Flow_treat_t500(idx_CA, idx_CA)) - (sum(Flow_ctrl_t500(:, idx_CA)) - Flow_ctrl_t500(idx_CA, idx_CA))) / pop_ctrl_teval(idx_CA) * 100;
    end
    save(data_file, 'psi_mults', 'fl_drops_dyn', 'fl_rel_stayers', 'fl_rel_inflow', ...
         'wa_drops_dyn', 'wa_rel_stayers', 'wa_rel_inflow', ...
         'ca_drops_dyn', 'ca_rel_stayers', 'ca_rel_inflow');
    fprintf('Saved data to %s\n', data_file);
end

%% 3. PLOT SEPARATE GRAPHS
    function plot_state(mults, drops, stayers, inflow, state_name, out_path)
        fig = figure('Position', [100 100 1000 600]);
        N = length(mults);
        
        b = bar(1:N, [stayers, inflow], 'stacked');

        
        hold on;
        plot(1:N, drops, 'o', 'MarkerSize', 1, 'DisplayName', 'Total Change');
        grid on;
        yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
        
        tick_idx = 1:5:N;
        if tick_idx(end) ~= N, tick_idx = [tick_idx, N]; end
        xticks(tick_idx);
        xticklabels(string(mults(tick_idx)));
        
        xlabel('Migration Cost Multiplier (\psi)');
        ylabel('Relative Pop Change (%)');
        legend({'Locals Fleeing', 'Others Not Coming', 'Total Change'}, 'Location', 'best', 'FontSize', 14);
        
       
        exportgraphics(fig, out_path, 'Resolution', 150);
    end

plot_state(psi_mults, fl_drops_dyn, fl_rel_stayers, fl_rel_inflow, 'Florida', fullfile(out_dir, 'ushape_t500_florida.png'));
plot_state(psi_mults, wa_drops_dyn, wa_rel_stayers, wa_rel_inflow, 'Washington', fullfile(out_dir, 'ushape_t500_washington.png'));
plot_state(psi_mults, ca_drops_dyn, ca_rel_stayers, ca_rel_inflow, 'California', fullfile(out_dir, 'ushape_t500_california.png'));

fprintf('\nDone saving separate graphs.\n');
end
