function plot_welfare_transition(T, mean_earn_t, mean_wealth_t, mean_ce_t, output_dir)
    set_aesthetics('subplot');
    t_arr = 1:T;
    
    % Between-State Std Dev (Spatial Inequality)
    var_earn_between = std(mean_earn_t, 0, 1);
    var_wealth_between = std(mean_wealth_t, 0, 1);
    var_ce_between = std(mean_ce_t, 0, 1);
    
    % Variation in % relative to t=1
    var_earn_pct = (var_earn_between - var_earn_between(1)) / var_earn_between(1) * 100;
    var_wealth_pct = (var_wealth_between - var_wealth_between(1)) / var_wealth_between(1) * 100;
    var_ce_pct = (var_ce_between - var_ce_between(1)) / var_ce_between(1) * 100;
    
    fig_welf = figure('Visible', 'off', 'Position', [100 100 1200 400]);

    % Spatial Welfare Inequality
    subplot(1,3,1);
    plot(t_arr, var_ce_pct, 'LineWidth', 3.0);
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std Dev (%)');
    set(gcf, 'Tag', 'welfare_inequality');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    % Between-State Income Variance
    subplot(1,3,2);
    plot(t_arr, var_earn_pct, 'LineWidth', 3.0);
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std Dev (%)');
    set(gcf, 'Tag', 'income_inequality');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    % Between-State Wealth Variance
    subplot(1,3,3);
    plot(t_arr, var_wealth_pct, 'LineWidth', 3.0);
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std Dev (%)');
    set(gcf, 'Tag', 'wealth_inequality');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig_welf, output_dir, 'welfare_transition.png');
end
