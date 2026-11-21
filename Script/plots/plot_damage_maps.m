function plot_damage_maps(P_SS1, P_SS2, output_dir, bounds)
    set_aesthetics('subplot');
    
    % Wages and housing as percentage changes
    delta_w = (P_SS2.w - P_SS1.w) ./ P_SS1.w * 100;
    delta_q = (P_SS2.q - P_SS1.q) ./ P_SS1.q * 100;
    
    % Map absolute util drop into Consumption-Equivalent percentage points
    delta_A_ce = (P_SS2.A - P_SS1.A) ./ P_SS1.gamma * 100;
    
    b_w = []; b_A = []; b_q = [];
    if nargin >= 4 && ~isempty(bounds)
        b_w = bounds.w;
        b_A = bounds.A;
        b_q = bounds.q;
    end

    fig_damage = figure('Visible', 'off', 'Position', [100 100 1400 450]);

    ax1 = subplot(1,3,1);
    set(gcf, 'Tag', 'damage_wage');
    plot_choropleth(P_SS2.state_names, delta_w, ax1, b_w);

    ax2 = subplot(1,3,2);
    set(gcf, 'Tag', 'damage_amenity');
    plot_choropleth(P_SS2.state_names, delta_A_ce, ax2, b_A);

    ax3 = subplot(1,3,3);
    set(gcf, 'Tag', 'damage_housing');
    plot_choropleth(P_SS2.state_names, delta_q, ax3, b_q);

    % sg% title('Climate Change Damages by State (SS2 vs SS1)', 'FontSize', 20, 'FontWeight', 'bold');
    save_fig(fig_damage, output_dir, 'damage_maps.png');
end