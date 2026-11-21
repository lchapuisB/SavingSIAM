function plot_peak_transition_maps(P_SS1, WY_t, mean_wealth_t, Lambda_trans, mean_w_SS1, Lambda_SS1, output_dir, bounds)
    set_aesthetics('subplot');
    [max_wy, t_peak] = max(WY_t);
    
    b_w = []; b_pop = [];
    if nargin >= 8 && ~isempty(bounds)
        b_w = bounds.w;
        b_pop = bounds.pop;
    end
    
    delta_W_peak = mean_wealth_t(:, t_peak) - mean_w_SS1;
    
    % Check dimensionality of Lambda_trans (could be J x Nz x T)
    if ndims(Lambda_trans) == 3
        pop_peak = squeeze(sum(Lambda_trans(:, :, t_peak), 2));
    else
        pop_peak = Lambda_trans(:, t_peak);
    end
    pop_SS1 = sum(Lambda_SS1, 2);
    
    delta_pop_peak = pop_peak - pop_SS1;

    fig_peak = figure('Visible', 'off', 'Position', [150 150 1000 450]);

    ax1 = subplot(1,2,1);
    set(gcf, 'Tag', 'wealth_accumulation_peak');
    plot_choropleth(P_SS1.state_names, delta_W_peak, ax1, b_w);

    ax2 = subplot(1,2,2);
    set(gcf, 'Tag', 'population_shift_peak');
    plot_choropleth(P_SS1.state_names, delta_pop_peak * 100, ax2, b_pop);

    % sg% title(sprintf('Precautionary Saving Peak (W/Y = %.2f at t=%d)', max_wy, t_peak), 'FontSize', 20, 'FontWeight', 'bold');
    save_fig(fig_peak, output_dir, 'peak_transition_maps.png');
end
