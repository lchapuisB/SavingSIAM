function plot_welfare_decomposition(P_SS1, Lambda_SS1, Lambda_SS2, a_bar_SS1, a_bar_cf, a_bar_cf_local, a_bar_SS2, Xi_SS1, Xi_cf, Xi_cf_local, Xi_SS2, output_dir, scenario_name, counterfactual)
    set_aesthetics();
% ============================================================
% Structural Shift-Share Decomposition of Welfare (CE)
% Isolates Mechanical, Local Behavioral, Spatial Behavioral, and Compositional
% ============================================================

    J = P_SS1.J;
    state_names = P_SS1.state_names;

    CE_baseline = zeros(J,1);
    CE_mech = zeros(J,1);
    CE_local = zeros(J,1);
    CE_spatial = zeros(J,1);
    CE_terminal = zeros(J,1);

    if strcmp(counterfactual, 'no_saving_struct')
        kappa = 0;
    else
        kappa = P_SS1.r / (1 + P_SS1.r);
    end

    CE1 = Xi_SS1;
    CE2_mech = Xi_cf;
    CE2_local =  Xi_cf_local;
    CE2_spatial = Xi_SS2;
    CE2_term =  Xi_SS2;

    for j = 1:J
        pop_1 = sum(Lambda_SS1(j, :));
        pop_2 = sum(Lambda_SS2(j, :));

        if pop_1 > 1e-15, CE_baseline(j) = sum(Lambda_SS1(j, :) .* CE1(j, :)) / pop_1; end

        if pop_1 > 1e-15, CE_mech(j) = sum(Lambda_SS1(j, :) .* CE2_mech(j, :)) / pop_1; end

        if pop_1 > 1e-15, CE_local(j) = sum(Lambda_SS1(j, :) .* CE2_local(j, :)) / pop_1; end

        if pop_1 > 1e-15, CE_spatial(j) = sum(Lambda_SS1(j, :) .* CE2_spatial(j, :)) / pop_1; end

        if pop_2 > 1e-15, CE_terminal(j) = sum(Lambda_SS2(j, :) .* CE2_term(j, :)) / pop_2; end
    end

    total_change = CE_terminal - CE_baseline;
    mech_margin = CE_mech - CE_baseline;
    local_margin = CE_local - CE_mech;
    spatial_margin = CE_spatial - CE_local;
    comp_margin = CE_terminal - CE_spatial;

    %% --- 1. Export Table ---
    T = table(state_names(:), CE_baseline(:), CE_terminal(:), total_change(:), mech_margin(:), local_margin(:), spatial_margin(:), comp_margin(:), ...
        'VariableNames', {'State', 'CE_Baseline', 'CE_Terminal', 'Total_Change', 'Mech_Margin', 'Local_Margin', 'Spatial_Margin', 'Comp_Margin'});
    writetable(T, fullfile(output_dir, sprintf('welfare_decomposition_%s.csv', scenario_name)));

    %% --- 2. Plot Stacked Bar Chart ---
    [~, sort_idx] = sort(total_change, 'descend');
    
    idx_pos = sort_idx(1:min(5, J));
    
    idx_neg = sort_idx(max(1, J-4):J);
    
    idx_top = unique([idx_pos; idx_neg], 'stable');
    top_n = length(idx_top);
    
    fig = figure('Visible', 'off', 'Position', [150 150 1000 600]);
    
    bar_data = [mech_margin(idx_top), local_margin(idx_top), spatial_margin(idx_top), comp_margin(idx_top)];
    
    b = bar(bar_data, 'stacked');
    b(1).FaceColor = [0.82, 0.42, 0.35]; % Terracotta (Mechanical)
    b(2).FaceColor = [0.85, 0.70, 0.30]; % Muted Gold (Local Behavioral)
    b(3).FaceColor = [0.35, 0.68, 0.48]; % Emerald Muted (Spatial Contagion)
    b(4).FaceColor = [0.25, 0.52, 0.72]; % Cobalt Blue (Compositional)
    
    xticks(1:top_n);
    xticklabels(state_names(idx_top));
    xtickangle(45);
    
    % Add total change markers
    hold on;
    plot(1:top_n, total_change(idx_top), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'DisplayName', 'Total \Delta CE');
    
    ylabel('\Delta Welfare (CE)');
    legend({'Mechanical', 'Local Behavioral', 'Spatial Contagion', 'Compositional', 'Total'}, ...
           'Location', 'southwest', 'NumColumns', 2, 'FontSize', 16, 'Box', 'off');
    
    grid on;
    
    save_fig(fig, output_dir, sprintf('welfare_decomposition_bar_%s.png', scenario_name));
end
