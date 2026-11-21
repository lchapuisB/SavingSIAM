function plot_dynamic_moments(T, TRANSITION_TYPE, agg_earn, agg_wealth, WY_t, agg_mobility, output_dir)
    set_aesthetics('subplot');
    t_arr = 1:T;
    col = [0.25, 0.52, 0.72];  % Cobalt Blue (first palette color)

    sub_dir = fullfile(output_dir, 'dynamic_moments');
    if ~exist(sub_dir, 'dir'), mkdir(sub_dir); end

    % Earnings
    fig1 = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, agg_earn, 'Color', col, 'LineWidth', 3.0);
    yline(agg_earn(1), ':k', 'HandleVisibility', 'off');
    grid on; xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); xlabel('Year'); ylabel('National Earnings');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig1, sub_dir, 'dynamic_moments_national_earnings.png');

    % Wealth
    fig2 = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, agg_wealth, 'Color', col, 'LineWidth', 3.0);
    yline(agg_wealth(1), ':k', 'HandleVisibility', 'off');
    grid on; xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); xlabel('Year'); ylabel('National Wealth');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig2, sub_dir, 'dynamic_moments_national_wealth.png');

    % W/Y
    fig3 = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, WY_t, 'Color', col, 'LineWidth', 3.0);
    yline(WY_t(1), ':k', 'HandleVisibility', 'off');
    grid on; xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); xlabel('Year'); ylabel('National W/Y');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig3, sub_dir, 'dynamic_moments_national_wy.png');

    % Migration
    fig4 = figure('Visible', 'off', 'Position', [150 150 600 400]); hold on;
    plot(t_arr, agg_mobility*100, 'Color', col, 'LineWidth', 3.0);
    yline(agg_mobility(1)*100, ':k', 'HandleVisibility', 'off');
    grid on; xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); xlabel('Year'); ylabel('Migration Rate (%)');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    save_fig(fig4, sub_dir, 'dynamic_moments_national_mobility.png');

    fprintf('  Dynamic moment plots saved to %s\n', sub_dir);
end