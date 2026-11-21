function plot_did_welfare_transition(T, T_shock, mean_earn_ctrl, mean_earn_treat, ...
    mean_wealth_ctrl, mean_wealth_treat, mean_ce_ctrl, mean_ce_treat, output_dir)
    set_aesthetics('subplot');
    t_arr = 1:T;
    
    % Between-State Std Dev (Spatial Inequality)
    var_earn_ctrl = std(mean_earn_ctrl, 0, 1);
    var_earn_treat = std(mean_earn_treat, 0, 1);
    
    var_wealth_ctrl = std(mean_wealth_ctrl, 0, 1);
    var_wealth_treat = std(mean_wealth_treat, 0, 1);
    
    var_ce_ctrl = std(mean_ce_ctrl, 0, 1);
    var_ce_treat = std(mean_ce_treat, 0, 1);
    
    % Variation in % relative to t=1
    var_earn_pct_ctrl = (var_earn_ctrl - var_earn_ctrl(1)) / var_earn_ctrl(1) * 100;
    var_earn_pct_treat = (var_earn_treat - var_earn_treat(1)) / var_earn_treat(1) * 100;
    
    var_wealth_pct_ctrl = (var_wealth_ctrl - var_wealth_ctrl(1)) / var_wealth_ctrl(1) * 100;
    var_wealth_pct_treat = (var_wealth_treat - var_wealth_treat(1)) / var_wealth_treat(1) * 100;
    
    var_ce_pct_ctrl = (var_ce_ctrl - var_ce_ctrl(1)) / var_ce_ctrl(1) * 100;
    var_ce_pct_treat = (var_ce_treat - var_ce_treat(1)) / var_ce_treat(1) * 100;

    % Pure DiD Diff
    did_var_earn = var_earn_pct_treat - var_earn_pct_ctrl;
    did_var_wealth = var_wealth_pct_treat - var_wealth_pct_ctrl;
    did_var_ce = var_ce_pct_treat - var_ce_pct_ctrl;
    
    % Plot 1: Levels Comparison (Control vs Treatment)
    fig_levels = figure('Visible', 'off', 'Position', [100 100 1200 400]);

    % Spatial Welfare Inequality
    subplot(1,3,1); hold on;
    plot(t_arr, var_ce_pct_ctrl, 'DisplayName', 'Control (No CC)');
    plot(t_arr, var_ce_pct_treat, 'DisplayName', 'Treatment (CC)');
    xline(T_shock, '--k', 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std Dev (%)'); title('Welfare Inequality');
    set(gcf, 'Tag', 'welfare_inequality');
    grid on; legend('Location', 'best');

    % Between-State Income Variance
    subplot(1,3,2); hold on;
    plot(t_arr, var_earn_pct_ctrl, 'DisplayName', 'Control');
    plot(t_arr, var_earn_pct_treat, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std Dev (%)'); title('Earnings Inequality');
    set(gcf, 'Tag', 'income_inequality');
    grid on;

    % Between-State Wealth Variance
    subplot(1,3,3); hold on;
    plot(t_arr, var_wealth_pct_ctrl, 'DisplayName', 'Control');
    plot(t_arr, var_wealth_pct_treat, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std Dev (%)'); title('Wealth Inequality');
    set(gcf, 'Tag', 'wealth_inequality');
    grid on;

    save_fig(fig_levels, output_dir, 'did_welfare_transition_levels.png');

    % Plot 2: Diff (DiD)
    % Spatial Welfare Inequality
    fig_diff_ce = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_var_ce);
    xline(T_shock, '--k', 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('Treat - Ctrl (%)'); title('Welfare Inequality Diff');
    set(gcf, 'Tag', 'did_welfare_inequality');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig_diff_ce, output_dir, 'did_welfare_transition_diff_ce.png');

    % Between-State Income Variance
    fig_diff_earn = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_var_earn);
    xline(T_shock, '--k', 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('Treat - Ctrl (%)'); title('Earnings Inequality Diff');
    set(gcf, 'Tag', 'did_income_inequality');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig_diff_earn, output_dir, 'did_welfare_transition_diff_earn.png');

    % Between-State Wealth Variance
    fig_diff_wealth = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, did_var_wealth);
    xline(T_shock, '--k', 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('Treat - Ctrl (%)'); title('Wealth Inequality Diff');
    set(gcf, 'Tag', 'did_wealth_inequality');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig_diff_wealth, output_dir, 'did_welfare_transition_diff_wealth.png');
end
