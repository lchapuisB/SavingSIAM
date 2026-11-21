function plot_xi_transition(P_2010, P_2100, T, TRANSITION_TYPE, Xi_trans, Xi_SS_2010, Xi_SS_2100, paths, output_dir)
    set_aesthetics('subplot');
    [~, idx_sort] = sort(P_2010.w);
    lo = idx_sort(1:5);
    hi = idx_sort(end-4:end);
    t_arr = 1:T;
    iz_mid = round(P_2010.Nz / 2);

    fig1 = figure('Visible', 'off', 'Position', [100 100 1200 800]);

    % A. Average Xi Evolution
    subplot(2,2,1); hold on;
    set(gcf, 'Tag', 'average_evolution');
    avg_xi = squeeze(mean(Xi_trans(:, iz_mid, :), 1));
    avg_xi_lo = squeeze(mean(Xi_trans(lo, iz_mid, :), 1));
    avg_xi_hi = squeeze(mean(Xi_trans(hi, iz_mid, :), 1));
    
    plot(t_arr, avg_xi, 'LineWidth', 3.0);
    plot(t_arr, avg_xi_lo, 'LineWidth', 3.0, 'LineStyle', '--');
    plot(t_arr, avg_xi_hi, 'LineWidth', 3.0, 'LineStyle', '--');
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Xi(z=0)');
    legend({'National Avg', 'Low-Wage States', 'High-Wage States'}, 'Location', 'best');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    % B. Spatial Std Dev of Xi Over Time
    subplot(2,2,2); hold on;
    set(gcf, 'Tag', 'spatial_std_dev');
    std_xi = squeeze(std(Xi_trans(:, iz_mid, :), 0, 1));
    std_xi_pct = (std_xi - std_xi(1)) / std_xi(1) * 100;
    plot(t_arr, std_xi_pct, 'LineWidth', 3.0);
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Delta Std(\Xi) (%)');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    % C. Xi as a function of z
    subplot(2,2,3); hold on;
    set(gcf, 'Tag', 'gradient_productivity');
    % Plot for a representative low and high state in 2010 and 2100
    state_lo = lo(1);
    state_hi = hi(1);
    plot(P_2010.z_grid, Xi_SS_2010(state_lo, :), 'LineStyle', '-');
    plot(P_2010.z_grid, Xi_SS_2100(state_lo, :), 'LineStyle', '--');
    plot(P_2010.z_grid, Xi_SS_2010(state_hi, :), 'LineStyle', '-');
    plot(P_2010.z_grid, Xi_SS_2100(state_hi, :), 'LineStyle', '--');
    
    xlabel('Productivity (z)'); ylabel('\Xi(z)');
    legend({'Low-Wage (2010)', 'Low-Wage (2100)', 'High-Wage (2010)', 'High-Wage (2100)'}, 'Location', 'best');
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    % D. Map of Delta Xi at z=0
    ax4 = subplot(2,2,4); hold on;
    set(gcf, 'Tag', 'map_delta_xi');
    delta_xi = Xi_SS_2100(:, iz_mid) - Xi_SS_2010(:, iz_mid);
    plot_choropleth(P_2100.state_names, delta_xi, ax4);
    
    % sg% title(sprintf('Transition 2010→2100 — T=%d, type=%s', T, TRANSITION_TYPE));
    save_fig(fig1, output_dir, 'xi_transition.png');
end