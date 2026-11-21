function plot_transition_dynamics(P_1990, P_2010, T, TRANSITION_TYPE, ...
    Xi_trans, Xi_SS_1990, Xi_SS_2010, paths, ...
    Lambda_SS_1990, Lambda_SS_2010, ...
    agg_earn, agg_wealth, WY_t, agg_mobility, mean_earn_t, ...
    mean_e_1990, mean_e_2010, wye_1990, wye_2010, ...
    Flow_jl, Z_mean_jl, output_dir)
    
    set_aesthetics('subplot');

    [~, idx_sort] = sort(P_2010.w);
    lo = idx_sort(1:5);
    hi = idx_sort(end-4:end);
    t_arr = 1:T;
    iz_mid = round(P_2010.Nz / 2);

    fig1 = figure('Visible', 'off', 'Position', [100 100 1200 800]);

    subplot(2,2,1); hold on;
    set(gcf, 'Tag', 'xi_transition_lo_hi');
    for i = 1:numel(lo)
        plot(t_arr, squeeze(Xi_trans(lo(i), iz_mid, :)), ...
             'Color',[0.25 0.52 0.72],'LineWidth',3.0);
    end
    for i = 1:numel(hi)
        plot(t_arr, squeeze(Xi_trans(hi(i), iz_mid, :)), ...
             'Color',[0.82 0.42 0.35],'LineWidth',3.0);
    end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Xi(z=0)');
    legend([plot(nan,nan), ...
            plot(nan,nan)], ...
           {'Bas salaire','Haut salaire'}, 'Location','best');
    grid on;

    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    subplot(2,2,2); hold on;
    set(gcf, 'Tag', 'xi_ss_scatter');
    scatter(P_2010.w, Xi_SS_1990(:,iz_mid), 40, [0.25 0.52 0.72], 'filled');
    scatter(P_2010.w, Xi_SS_2010(:,iz_mid), 40, [0.82 0.42 0.35], 's','filled');
    xlabel('w_{ell}'); ylabel('\Xi(z=0)');
    legend({'SS 1990','SS 2010'},'Location','best'); grid on;

    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    subplot(2,2,3); hold on;
    set(gcf, 'Tag', 'xi_convergence');
    diff_path = squeeze(Xi_trans(:, iz_mid, :)) - Xi_SS_2010(:, iz_mid);
    for j = 1:5:P_2010.J
        plot(t_arr, diff_path(j,:),'LineWidth',3.0);
    end
    yline(0,'k--','LineWidth',3.0);
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('\Xi_t - \Xi^{SS}_{2010}');

    subplot(2,2,4); hold on;
    set(gcf, 'Tag', 'wage_path');
    for i = 1:min(3,numel(lo))
        plot(t_arr, paths.w(lo(i),:),'LineWidth',3.0);
    end
    for i = 1:min(3,numel(hi))
        plot(t_arr, paths.w(hi(i),:),'LineWidth',3.0);
    end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('w_{ell}(t)');

    save_fig(fig1, output_dir, 'xi_transition.png');

    fig2 = figure('Visible', 'off', 'Position', [150 150 1400 400]);

    subplot(1,4,1);
    set(gcf, 'Tag', 'agg_earnings');
    plot(t_arr, agg_earn, 'LineWidth', 3.0); grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('Earnings aggrégées'); 

    subplot(1,4,2);
    set(gcf, 'Tag', 'agg_wealth');
    plot(t_arr, agg_wealth, 'LineWidth', 3.0); grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('Richesse agrégée'); % title('Richesse');

    subplot(1,4,3);
    set(gcf, 'Tag', 'national_wy');
    plot(t_arr, WY_t,'LineWidth',3.0);
    grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('W/Y national'); % title('Ratio Richesse/Revenu');

    subplot(1,4,4);
    set(gcf, 'Tag', 'mobility_rate');
    plot(t_arr, agg_mobility*100, 'LineWidth', 3.0);
     grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('Taux de mobilité (%)'); % title('Mobilité');

    % sg% title(sprintf('Moments dynamiques — %s', TRANSITION_TYPE));
    save_fig(fig2, output_dir, 'dynamic_moments.png');

    fig3 = figure('Visible', 'off', 'Position', [200 200 900 400]);
    subplot(1,2,1);
    set(gcf, 'Tag', 'std_log_earnings');
    plot(t_arr, std(log(mean_earn_t),0,1), 'LineWidth', 3.0); grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('std(log earnings) entre états');
    % title('Dispersion des salaires');

    subplot(1,2,2);
    set(gcf, 'Tag', 'cv_earnings');
    plot(t_arr, std(log(mean_earn_t),0,1) ./ mean(log(mean_earn_t),1), ...
         'Color', [0.82, 0.42, 0.35], 'LineWidth', 3.0); grid on;
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end
    xline(90, '--k', 'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    xlabel('Year'); ylabel('CV des earnings'); 
    save_fig(fig3, output_dir, 'spatial_convergence.png');

    delta_w    = P_2010.w - P_1990.w;
    delta_pop  = sum(Lambda_SS_2010, 2) - sum(Lambda_SS_1990, 2);
    delta_earn = mean_e_2010 - mean_e_1990;
    delta_wy   = wye_2010 - wye_1990;

    fig_map = figure('Visible', 'off', 'Position', [100 100 1400 900]);
    n_half = 128;
    t_val = linspace(0, 1, n_half)';
    gamma = 3; 
    
    r_red = 0.75 + 0.23 * (t_val.^gamma);
    g_red = 0.30 + 0.68 * (t_val.^gamma);
    b_red = 0.20 + 0.70 * (t_val.^gamma);
    r_to_w = [r_red, g_red, b_red];
    
    r_green = 0.25 + 0.73 * ((1-t_val).^gamma);
    g_green = 0.45 + 0.53 * ((1-t_val).^gamma);
    b_green = 0.15 + 0.75 * ((1-t_val).^gamma);
    w_to_g = [r_green, g_green, b_green];
    
    colormap_choice = [r_to_w; w_to_g(2:end, :)]; 
        
    ax1 = subplot(2,2,1);
    set(gcf, 'Tag', 'map_delta_w');
    plot_choropleth(P_2010.state_names, delta_earn, ax1, colormap_choice);
        
    ax2 = subplot(2,2,2);
    set(gcf, 'Tag', 'map_delta_pop');
    plot_choropleth(P_2010.state_names, delta_pop * 100, ax2, colormap_choice);
        
    ax3 = subplot(2,2,3);
    set(gcf, 'Tag', 'map_delta_earnings');
    plot_choropleth(P_2010.state_names, delta_earn, ax3, colormap_choice);
        
    ax4 = subplot(2,2,4);
    set(gcf, 'Tag', 'map_delta_wy');
    plot_choropleth(P_2010.state_names, delta_wy, ax4, colormap_choice);
        
    save_fig(fig_map, output_dir, 'transition_maps.png');

    J = P_2010.J;
    mean_z_resident = zeros(J, 1);
    for l = 1:J
        pop_l = sum(Lambda_SS_2010(l, :));
        if pop_l > 1e-15
            mean_z_resident(l) = sum(Lambda_SS_2010(l, :) .* P_2010.z_grid') / pop_l;
        end
    end

    in_migrants_z = zeros(J, 1);
    in_flow_total = zeros(J, 1);

    for l = 1:J
        for j = 1:J
            if j ~= l
                in_flow_total(l) = in_flow_total(l) + Flow_jl(j, l, T);
                in_migrants_z(l) = in_migrants_z(l) + Flow_jl(j, l, T) * Z_mean_jl(j, l, T);
            end
        end
        if in_flow_total(l) > 1e-15
            in_migrants_z(l) = in_migrants_z(l) / in_flow_total(l);
        end
    end

    fig_sorting = figure('Visible', 'off', 'Position', [150 150 1000 450]);

    subplot(1,2,1);
    set(gcf, 'Tag', 'std_log_earnings');
    set(gcf, 'Tag', 'std_log_earnings'); hold on;
    scatter(P_2010.w, mean_z_resident, 50, 'filled', 'MarkerFaceColor', [0.25 0.52 0.72]);
    p1 = polyfit(P_2010.w, mean_z_resident, 1);
    x_fit = linspace(min(P_2010.w), max(P_2010.w), 100);
    plot(x_fit, polyval(p1, x_fit), 'k--', 'LineWidth', 3.0);
    xlabel('Salaire local (w_l)'); ylabel('Productivité moyenne des résidents E[z | l]');
    

    subplot(1,2,2);
    set(gcf, 'Tag', 'cv_earnings');
    set(gcf, 'Tag', 'cv_earnings'); hold on;
    scatter(P_2010.w, in_migrants_z, 50, 'filled', 'MarkerFaceColor', [0.82 0.42 0.35]);
    p2 = polyfit(P_2010.w, in_migrants_z, 1);
    plot(x_fit, polyval(p2, x_fit), 'k--', 'LineWidth', 3.0);
    xlabel('Salaire local (w_l)'); ylabel('Productivité moyenne des immigrants vers l');
   
    save_fig(fig_sorting, output_dir, 'spatial_sorting.png');

    sorting_1990 = zeros(P_1990.J, 1);
    sorting_2010 = zeros(P_2010.J, 1);

    pop_1990 = sum(Lambda_SS_1990, 2);
    for l = 1:P_1990.J
        if pop_1990(l) > 1e-15
            sorting_1990(l) = sum(P_1990.z_grid' .* Lambda_SS_1990(l, :)) / pop_1990(l);
        end
    end

    pop_2010 = sum(Lambda_SS_2010, 2);
    for l = 1:P_2010.J
        if pop_2010(l) > 1e-15
            sorting_2010(l) = sum(P_2010.z_grid' .* Lambda_SS_2010(l, :)) / pop_2010(l);
        end
    end

    delta_sorting = sorting_2010 - sorting_1990;

    fig_sorting_map = figure('Visible', 'off', 'Position', [200 200 900 600]);
    plot_choropleth(P_2010.state_names, delta_sorting, ...
        '\Delta Tri Spatial (\Delta E[z|l]) entre 1990 et 2010', gca, colormap_choice);
    save_fig(fig_sorting_map, output_dir, 'map_delta_sorting.png');
end