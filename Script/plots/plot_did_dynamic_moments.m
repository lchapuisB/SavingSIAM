function plot_did_dynamic_moments(T, T_shock, agg_earn_ctrl, agg_earn_treat, agg_wealth_ctrl, agg_wealth_treat, WY_ctrl, WY_treat, agg_mobility_ctrl, agg_mobility_treat, output_dir)
    set_aesthetics('subplot');
% ============================================================
% Plots Dynamic Moments for DiD (Levels and Differences)
% ============================================================
    t_arr = 1:T;
    
    % 1. Compute deviations from initial level for plotting
    earn_ctrl_dev = agg_earn_ctrl - agg_earn_ctrl(1);
    earn_treat_dev = agg_earn_treat - agg_earn_treat(1);
    
    wealth_ctrl_dev = agg_wealth_ctrl - agg_wealth_ctrl(1);
    wealth_treat_dev = agg_wealth_treat - agg_wealth_treat(1);
    
    wy_ctrl_dev = WY_ctrl - WY_ctrl(1);
    wy_treat_dev = WY_treat - WY_treat(1);
    
    mob_ctrl_dev = agg_mobility_ctrl - agg_mobility_ctrl(1);
    mob_treat_dev = agg_mobility_treat - agg_mobility_treat(1);

    % 2. Compute Pure DiD (Treatment - Control)
    did_earn = agg_earn_treat - agg_earn_ctrl;
    did_wealth = agg_wealth_treat - agg_wealth_ctrl;
    did_wy = WY_treat - WY_ctrl;
    did_mob = agg_mobility_treat - agg_mobility_ctrl;

    %% FIGURES (Deviations from t=1 Level)
    % Earnings
    fig1a = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, earn_ctrl_dev, 'b-', 'LineWidth', 3.0, 'DisplayName', 'Control (No CC)');
    plot(t_arr, earn_treat_dev, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Treatment (With CC)');
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Year'); ylabel('\Delta Earnings (from t=1)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig1a, output_dir, 'did_dynamic_moments_levels_earnings.png');

    % Wealth
    fig1b = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, wealth_ctrl_dev, 'b-', 'LineWidth', 3.0, 'DisplayName', 'Control');
    plot(t_arr, wealth_treat_dev, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Year'); ylabel('\Delta Wealth (from t=1)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig1b, output_dir, 'did_dynamic_moments_levels_wealth.png');

    % W/Y
    fig1c = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, wy_ctrl_dev, 'b-', 'LineWidth', 3.0, 'DisplayName', 'Control');
    plot(t_arr, wy_treat_dev, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Year'); ylabel('\Delta W/Y (from t=1)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig1c, output_dir, 'did_dynamic_moments_levels_wy.png');

    % Migration
    fig1d = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, mob_ctrl_dev*100, 'b-', 'LineWidth', 3.0, 'DisplayName', 'Control');
    plot(t_arr, mob_treat_dev*100, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'LineWidth', 3.0, 'HandleVisibility', 'off'); yline(0, ':k', 'HandleVisibility', 'off');
    grid on; xlabel('Year'); ylabel('\Delta Migration Rate (from t=1, %)'); legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig1d, output_dir, 'did_dynamic_moments_levels_migration.png');

    %% FIGURES: Pure DiD (Treatment - Control)
    % Earnings
    fig2a = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_earn, 'LineWidth', 3.0);
    xline(T_shock, '--k', 'LineWidth', 3.0); yline(0, ':k');
    grid on; xlabel('Year'); ylabel('Treat - Ctrl');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig2a, output_dir, 'did_dynamic_moments_diff_earnings.png');

    % Wealth
    fig2b = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_wealth, 'LineWidth', 3.0);
    xline(T_shock, '--k', 'LineWidth', 3.0); yline(0, ':k');
    grid on; xlabel('Year'); ylabel('Treat - Ctrl');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig2b, output_dir, 'did_dynamic_moments_diff_wealth.png');

    % W/Y
    fig2c = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_wy, 'LineWidth', 3.0);
    xline(T_shock, '--k', 'LineWidth', 3.0); yline(0, ':k');
    grid on; xlabel('Year'); ylabel('Treat - Ctrl');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig2c, output_dir, 'did_dynamic_moments_diff_wy.png');

    % Migration
    fig2d = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_mob*100, 'LineWidth', 3.0);
    xline(T_shock, '--k', 'LineWidth', 3.0); yline(0, ':k');
    grid on; xlabel('Year'); ylabel('Treat - Ctrl (%)');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig2d, output_dir, 'did_dynamic_moments_diff_migration.png');

    fprintf('  Dynamic moment plots saved to %s\n', output_dir);
end
