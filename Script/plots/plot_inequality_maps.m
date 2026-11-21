function plot_inequality_maps(state_names, delta_CE_var, delta_W_var, output_dir, bounds)
    set_aesthetics('subplot');
    
    b_ce_var = []; b_w_var = [];
    if nargin >= 5 && ~isempty(bounds)
        b_ce_var = bounds.ce_var;
        b_w_var = bounds.w_var;
    end

    % Combined figure: Welfare Inequality and Wealth Inequality
    fig_inequality = figure('Visible', 'off', 'Position', [100 100 1600 600]);
    
    % Figure 2: Welfare Inequality
    ax1 = subplot(1, 2, 1);
    set(gcf, 'Tag', 'welfare_inequality');
    plot_choropleth(state_names, delta_CE_var, ax1, b_ce_var);
    
    % Figure 3: Wealth Inequality
    ax2 = subplot(1, 2, 2);
    set(gcf, 'Tag', 'wealth_inequality');
    plot_choropleth(state_names, delta_W_var, ax2, b_w_var);
    
    % sg% title('Spatial Comparative Statics: Inequality Changes', 'FontSize', 20, 'FontWeight', 'bold');
    save_fig(fig_inequality, output_dir, 'welfare_wealth_inequality_maps.png');
end
