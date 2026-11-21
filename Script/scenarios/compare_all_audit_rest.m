function compare_all_audit(res_A, res_B, res_C, res_D, res_E, res_F, res_G, res_H, res_I)
% ============================================================
% Cross-Scenario Macro & Welfare Comparison (Audit)
% Combined version with new aesthetics
% ============================================================

    set_aesthetics('subplot');
    output_dir_macro = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Output', 'graph_macro_audit');
    if ~exist(output_dir_macro, 'dir'); mkdir(output_dir_macro); end
    output_dir_welfare = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Output', 'graph_welfare_audit');
    if ~exist(output_dir_welfare, 'dir'); mkdir(output_dir_welfare); end

    %% ── 1. COMPUTE TRAJECTORIES ─────────────────────────────────────────────
    scenarios_E = {'A: Perfect Foresight', 'B: Myopia', 'F: No Saving'};
    all_res_E = {res_A, res_B, res_F};

    scenarios_W = {'A: Perfect Foresight', 'B: Myopia'};
    all_res_W = {res_A, res_B};
    
    did_res = {res_C};
    did_names = {'C: DID No Mig'};

    T = length(res_A.agg_earn);
    T_shock = 90;
    t_vec = 1:T;

    % Earning, Mobility, CE preallocate
    agg_E = zeros(length(all_res_E), T);
    base_E = zeros(length(all_res_E), 1);
    agg_M = zeros(length(all_res_E), T);
    agg_CE = zeros(length(all_res_E), T);

    for s = 1:length(all_res_E)
        res = all_res_E{s};
        agg_E(s, :) = res.agg_earn';
        base_E(s) = res.agg_earn(1);
        
        if isfield(res, 'agg_mobility')
            agg_M(s, :) = res.agg_mobility';
        else
            agg_M(s, :) = zeros(1, T); 
        end
        
        pop = squeeze(sum(res.Lambda_trans, 2)); 
        if isfield(res, "mean_ce_t"); agg_CE(s, :) = sum(res.mean_ce_t .* pop, 1); else; agg_CE(s, :) = NaN(1, T); end
    end

    % Wealth & Saving Rate preallocate
    agg_W = zeros(length(all_res_W), T);
    agg_SR = zeros(length(all_res_W), T);
    base_W = zeros(length(all_res_W), 1);
    for s = 1:length(all_res_W)
        res = all_res_W{s};
        agg_W(s, :) = res.agg_wealth';
        saving = [0, diff(res.agg_wealth')];
        agg_SR(s, :) = saving ./ res.agg_earn' * 100;
        base_W(s) = res.agg_wealth(1);
    end

    % DID arms
    T_did = length(res_C.agg_earn_ctrl);
    agg_E_did = zeros(length(did_res), T_did);
    pct_E_did = zeros(length(did_res), T_did);
    agg_W_did = zeros(length(did_res), T_did);
    pct_W_did = zeros(length(did_res), T_did);
    agg_SR_did = zeros(length(did_res), T_did);
    agg_M_did = zeros(length(did_res), T_did);
    agg_CE_did_diff = zeros(length(did_res), T_did);
    pct_CE_did_diff = zeros(length(did_res), T_did);

    for s = 1:length(did_res)
        res = did_res{s};
        
        agg_E_did(s, :) = (res.agg_earn_treat' - res.agg_earn_ctrl');
        pct_E_did(s, :) = (res.agg_earn_treat' - res.agg_earn_ctrl') / res.agg_earn_ctrl(1) * 100;
        
        agg_W_did(s, :) = (res.agg_wealth_treat' - res.agg_wealth_ctrl');
        pct_W_did(s, :) = (res.agg_wealth_treat' - res.agg_wealth_ctrl') / res.agg_wealth_ctrl(1) * 100;
       
        saving_treat = [0, diff(res.agg_wealth_treat')];
        SR_treat = saving_treat ./ res.agg_earn_treat' * 100;
        saving_ctrl = [0, diff(res.agg_wealth_ctrl')];
        SR_ctrl = saving_ctrl ./ res.agg_earn_ctrl' * 100;
        agg_SR_did(s, :) = SR_treat - SR_ctrl;
        
        if isfield(res, 'agg_mobility_treat') && isfield(res, 'agg_mobility_ctrl')
            agg_M_did(s, :) = res.agg_mobility_treat' - res.agg_mobility_ctrl';
        end
        
        pop_ctrl  = squeeze(sum(res.Lambda_trans_ctrl, 2));
        pop_treat = squeeze(sum(res.Lambda_trans_treat, 2));
        ce_ctrl  = sum(res.mean_ce_ctrl  .* pop_ctrl,  1);
        ce_treat = sum(res.mean_ce_treat .* pop_treat, 1);
        
        agg_CE_did_diff(s, :) = ce_treat - ce_ctrl;
        % Denominator: aggregate baseline material consumption (no amenities)
        kappa_did = res.map_P_2010.r / (1 + res.map_P_2010.r);
        base_flow_did = res.map_P_2010.y_tilde + kappa_did * res.mean_wealth_ctrl(:, 1);
        den_did = sum(base_flow_did .* pop_ctrl(:, 1));
        pct_CE_did_diff(s, :) = (ce_treat - ce_ctrl) / abs(den_did) * 100;
    end

    colors = get(groot, 'defaultAxesColorOrder');
    N_colors = size(colors, 1);
    
    get_color = @(idx) colors(mod(idx - 1, N_colors) + 1, :);

    %% ── 2. MACRO: EARNING TRAJECTORIES ─────────────────────────────────────
    fig1 = figure('Visible', 'off', 'Position', [100 100 1400 600]);

    subplot(1, 2, 1); hold on; set(gcf, 'Tag', 'macro_earning_delta_earning');
    for s = 1:length(scenarios_E)
        plot(t_vec, agg_E(s, :) - base_E(s), 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, agg_E_did(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_E)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('\Delta Earning (Level)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    subplot(1, 2, 2); hold on; set(gcf, 'Tag', 'macro_earning_delta_earning_pct');
    for s = 1:length(scenarios_E)
        plot(t_vec, (agg_E(s, :) - base_E(s)) / base_E(s) * 100, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, pct_E_did(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_E)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('\Delta Earning (%)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig1, output_dir_macro, 'macro_earning.png');

    %% ── 3. MACRO: WEALTH TRAJECTORIES ─────────────────────────────────────
    fig2 = figure('Visible', 'off', 'Position', [100 100 1400 600]);

    subplot(1, 2, 1); hold on; set(gcf, 'Tag', 'macro_wealth_delta_wealth');
    for s = 1:length(scenarios_W)
        plot(t_vec, agg_W(s, :) - base_W(s), 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_W{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, agg_W_did(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_W)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('\Delta Wealth (Level)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    subplot(1, 2, 2); hold on; set(gcf, 'Tag', 'macro_wealth_delta_wealth_pct');
    for s = 1:length(scenarios_W)
        plot(t_vec, (agg_W(s, :) - base_W(s)) / base_W(s) * 100, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_W{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, pct_W_did(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_W)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('\Delta Wealth (%)');  legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig2, output_dir_macro, 'macro_wealth.png');

    %% ── 3.5 MACRO: WEALTH-TO-EARNINGS (W/Y) TRAJECTORIES ───────────────────
    fig2b = figure('Visible', 'off', 'Position', [100 100 1400 600]);

    subplot(1, 2, 1); hold on; set(gcf, 'Tag', 'macro_wy_delta_wy');
    for s = 1:length(scenarios_W)
        res = all_res_W{s};
        wy = res.agg_wealth' ./ res.agg_earn';
        base_wy = res.agg_wealth(1) / res.agg_earn(1);
        plot(t_vec, wy - base_wy, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_W{s});
    end
    for s = 1:length(did_res)
        res = did_res{s};
        wy_treat = res.agg_wealth_treat' ./ res.agg_earn_treat';
        wy_ctrl = res.agg_wealth_ctrl' ./ res.agg_earn_ctrl';
        plot(1:T_did, wy_treat - wy_ctrl, '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_W)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('\Delta W/Y (Level)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    subplot(1, 2, 2); hold on; set(gcf, 'Tag', 'macro_wy_delta_wy_pct');
    for s = 1:length(scenarios_W)
        res = all_res_W{s};
        wy = res.agg_wealth' ./ res.agg_earn';
        base_wy = res.agg_wealth(1) / res.agg_earn(1);
        plot(t_vec, (wy - base_wy) / base_wy * 100, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_W{s});
    end
    for s = 1:length(did_res)
        res = did_res{s};
        wy_treat = res.agg_wealth_treat' ./ res.agg_earn_treat';
        wy_ctrl = res.agg_wealth_ctrl' ./ res.agg_earn_ctrl';
        base_wy_ctrl = res.agg_wealth_ctrl(1) / res.agg_earn_ctrl(1);
        plot(1:T_did, (wy_treat - wy_ctrl) / base_wy_ctrl * 100, '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_W)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('\Delta W/Y (%)');  legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig2b, output_dir_macro, 'macro_wy.png');

    %% ── 4. MACRO: SAVING RATE TRAJECTORIES ─────────────────────────────────
    fig3 = figure('Visible', 'off', 'Position', [100 100 700 600]); hold on;

    for s = 1:length(scenarios_W)
        plot(t_vec, agg_SR(s, :), 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_W{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, agg_SR_did(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_W)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('Saving Rate (%)'); 
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig3, output_dir_macro, 'macro_saving.png');

    %% ── 5. MACRO: MIGRATION (MOBILITY) TRAJECTORIES ─────────────────────────
    fig4 = figure('Visible', 'off', 'Position', [100 100 700 600]); hold on;

    for s = 1:length(scenarios_E)
        plot(t_vec, agg_M(s, :), 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, agg_M_did(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_E)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time'); ylabel('Total Mass of Migrants'); 
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig4, output_dir_macro, 'macro_migration.png');

    %% ── 6. WELFARE TRAJECTORIES ────────────────────────────
    fig_ce_level = figure('Visible', 'off', 'Position', [100 100 600 600]);
    hold on; set(gcf, 'Tag', 'welfare_audit_all_scenarios_delta_aggregate_ce_level');
    for s = 1:length(scenarios_E)
        res = all_res_E{s};
        if isfield(res, "agg_CE_SS1"); level_change = agg_CE(s, :) - res.agg_CE_SS1; else; level_change = agg_CE(s, :) - agg_CE(s, 1); end
        plot(t_vec, level_change, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
    end
    for s = 1:length(did_res)
        plot(1:T_did, agg_CE_did_diff(s, :), '--', 'LineWidth', 3.0, 'Color', get_color(length(scenarios_E)+s), 'DisplayName', did_names{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time (year)'); ylabel('\Delta Aggregate CE (Level)'); legend('Location', 'best', 'FontSize', 10, 'Box', 'off');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig_ce_level, output_dir_welfare, 'welfare_audit_all_scenarios_ce_level.png');

    % Graph 2: Consumption-Equivalent Welfare Change (%)
    % Denominator is aggregate material consumption (y_tilde + kappa*a),
    % NOT total CE (which includes unidentified amenities A_j).
    fig_ce_pct = figure('Visible', 'off', 'Position', [100 100 600 600]);
    hold on; set(gcf, 'Tag', 'welfare_audit_all_scenarios_delta_pct_aggregate_ce');
    for s = 1:length(scenarios_E)
        res = all_res_E{s};
        if isfield(res, "agg_CE_SS1"); base_ce = res.agg_CE_SS1; else; base_ce = agg_CE(s, 1); end
        
        % Compute empirically-identified denominator: agg material consumption at t=1
        kappa_s = res.map_P_2010.r / (1 + res.map_P_2010.r);
        base_flow_cons = res.paths.y_tilde(:, 1) + kappa_s * res.mean_wealth_t(:, 1); % J x 1
        pop_0 = sum(res.Lambda_trans(:, :, 1), 2); % J x 1
        den_mat_cons = sum(base_flow_cons .* pop_0);
        
        pct_change = zeros(size(agg_CE(s, :)));
        if den_mat_cons ~= 0
            pct_change = (agg_CE(s, :) - base_ce) ./ abs(den_mat_cons) * 100;
        end
        plot(t_vec, pct_change, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time (year)'); ylabel('\Delta Aggregate CE (% of material consumption)'); legend('Location', 'best', 'FontSize', 10, 'Box', 'off');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig_ce_pct, output_dir_welfare, 'welfare_audit_all_scenarios_ce_pct.png');

    % Graph 3: Atkinson Social Welfare
    fig_atkinson = figure('Visible', 'off', 'Position', [100 100 600 600]);
    hold on; set(gcf, 'Tag', 'welfare_audit_all_scenarios_atkinson_welfare');
    epsilon = 1.0; % Log utilitarian
    for s = 1:length(scenarios_E)
        res = all_res_E{s};
        if isfield(res, 'paths') && isfield(res, 'map_P_2010')
            J = res.map_P_2010.J;
            Nz = res.map_P_2010.Nz;
            kappa = res.map_P_2010.r / (1 + res.map_P_2010.r);
            z_grid = res.map_P_2010.z_grid(:)';
            
            W_atk = zeros(1, T);
            for t = 1:T
                y_til_t = res.paths.y_tilde(:, t);
                w_til_t = res.paths.w_tilde(:, t);
                a_t = res.mean_wealth_t(:, t); % J x 1 approximation
                
                c_flow_mat = repmat(y_til_t, 1, Nz) + w_til_t * z_grid + repmat(kappa * a_t, 1, Nz);
                
                % Floor to small positive value to avoid log(0) if any negative
                c_flow_mat = max(c_flow_mat, 1e-4);
                
                W_atk(t) = compute_atkinson_welfare(c_flow_mat, res.Lambda_trans(:,:,t), epsilon);
            end
            % Plot percentage change relative to t=1
            plot(t_vec, (W_atk - W_atk(1)) / W_atk(1) * 100, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
        end
    end
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Time (year)'); ylabel('\Delta Atkinson Welfare (%)'); legend('Location', 'best', 'FontSize', 10, 'Box', 'off');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig_atkinson, output_dir_welfare, 'welfare_audit_all_scenarios_atkinson.png');

    %% ── 7. ATKINSON WELFARE SENSITIVITY (EPSILON) ────────────────────────
    fig6 = figure('Visible', 'off', 'Position', [100 100 1800 600]);
    epsilons = [0, 1.0, 2.0];
    titles = {'Linear Dollar (\epsilon = 0)', 'Log Utilitarian (\epsilon = 1)', 'High Inequality Aversion (\epsilon = 2)'};
    
    for e = 1:3
        eps_val = epsilons(e);
        subplot(1, 3, e); hold on; set(gcf, 'Tag', sprintf('atkinson_welfare_eps_%d', e));
        for s = 1:length(scenarios_E)
            res = all_res_E{s};
            if isfield(res, 'paths') && isfield(res, 'map_P_2010')
                J = res.map_P_2010.J;
                Nz = res.map_P_2010.Nz;
                kappa = res.map_P_2010.r / (1 + res.map_P_2010.r);
                z_grid = res.map_P_2010.z_grid(:)';
                
                W_atk = zeros(1, T);
                for t = 1:T
                    y_til_t = res.paths.y_tilde(:, t);
                    w_til_t = res.paths.w_tilde(:, t);
                    a_t = res.mean_wealth_t(:, t);
                    
                    c_flow_mat = repmat(y_til_t, 1, Nz) + w_til_t * z_grid + repmat(kappa * a_t, 1, Nz);
                    c_flow_mat = max(c_flow_mat, 1e-4);
                    
                    W_atk(t) = compute_atkinson_welfare(c_flow_mat, res.Lambda_trans(:,:,t), eps_val);
                end
                plot(t_vec, (W_atk - W_atk(1)) / W_atk(1) * 100, 'LineWidth', 3.0, 'Color', get_color(s), 'DisplayName', scenarios_E{s});
            end
        end
        xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
        grid on; xlabel('Time (year)'); ylabel('\Delta Atkinson Welfare (%)'); title(titles{e}); legend('Location', 'best');
        yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    end

    save_fig(fig6, output_dir_welfare, 'welfare_audit_atkinson_sensitivity.png');

    fprintf('\nAudit figures saved to %s and %s\n', output_dir_macro, output_dir_welfare);
end
