function plot_transition_maps( P_SS2, Lambda_SS1, Lambda_SS2, mean_W1, mean_W2, wye_SS1, wye_SS2, delta_CE_mean, output_dir, bounds)
    set_aesthetics('subplot');
    pop_SS1 = sum(Lambda_SS1, 2);
    pop_SS2 = sum(Lambda_SS2, 2);
    denom_W = abs(mean_W2) + abs(mean_W1);
    delta_W_pct = zeros(size(mean_W1));
    valid_W = denom_W > 1e-8;
    delta_W_pct(valid_W) = 2 * (mean_W2(valid_W) - mean_W1(valid_W)) ./ denom_W(valid_W) * 100;
    
    denom_pop = abs(pop_SS2) + abs(pop_SS1);
    delta_pop_pct = zeros(size(pop_SS1));
    valid_pop = denom_pop > 1e-8;
    delta_pop_pct(valid_pop) = 2 * (pop_SS2(valid_pop) - pop_SS1(valid_pop)) ./ denom_pop(valid_pop) * 100;
    delta_wy   = wye_SS2 - wye_SS1;
    
    b_W_pct = []; b_pop_pct = []; b_wy = []; b_ce = [];
    if nargin >= 10 && ~isempty(bounds)
        if isfield(bounds, 'W_pct')
            b_W_pct = bounds.W_pct;
        end
        b_pop_pct = bounds.pop_pct;
        b_wy = bounds.wy;
        b_ce = bounds.ce;
    end

    fig_map = figure('Visible', 'off', 'Position', [100 100 2000 1400]);
        
    ax1 = subplot(2,2,1);

    set(gcf, 'Tag', 'pop_changes_pp');
    if any(mean_W1 ~= 0) 
        plot_choropleth(P_SS2.state_names, delta_W_pct, ax1, b_W_pct);
    else
        delta_W = mean_W2 - mean_W1;
        plot_choropleth(P_SS2.state_names, delta_W, ax1);
    end

    ax2 = subplot(2,2,2);
    set(gcf, 'Tag', 'local_pop_changes_pct');
    plot_choropleth(P_SS2.state_names, delta_pop_pct, ax2, b_pop_pct);

    ax3 = subplot(2,2,3);
    set(gcf, 'Tag', 'wealth_income_ratio_changes');
    plot_choropleth(P_SS2.state_names, delta_wy, ax3, b_wy);
    
    ax4 = subplot(2,2,4);
    set(gcf, 'Tag', 'delta_mean_welfare_ce');
    plot_choropleth(P_SS2.state_names, delta_CE_mean, ax4, b_ce);
        
    % sg% title('Spatial Comparative Statics: Long-Run Steady State Differences (SS2 vs SS1)', 'FontSize', 20, 'FontWeight', 'bold');
    save_fig(fig_map, output_dir, 'transition_maps.png');
end