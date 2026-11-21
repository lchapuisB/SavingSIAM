function plot_spatial_sorting(P, Lambda, out_dir)
    % P is a struct containing w and z_grid
    % Lambda is a J x Nz matrix of population distribution over z in each state
    % out_dir is the directory to save the plot
    
    J = P.J;
    z_grid = P.z_grid(:)'; % ensure it is a row vector
    
    % Compute expected z for each state
    E_z = zeros(J, 1);
    for j = 1:J
        mass_j = sum(Lambda(j, :));
        if mass_j > 1e-10
            E_z(j) = sum(Lambda(j, :) .* z_grid) / mass_j;
        end
    end
    
    % Log wages
    log_w = log(P.w(:));
    
    % Regression E_z on log_w
    X = [ones(J, 1), log_w];
    Y = E_z;
    beta = (X'*X) \ (X'*Y);
    
    % Create figure
    fig = figure('Name', 'Spatial Sorting', 'Visible', 'off');
    scatter(log_w, E_z, 100, 'o', 'filled', 'MarkerFaceColor', [0.7 0.2 0.2], 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Plot regression line
    w_min = min(log_w);
    w_max = max(log_w);
    w_grid = linspace(w_min, w_max, 100)';
    z_fit = beta(1) + beta(2) * w_grid;
    p1 = plot(w_grid, z_fit, '--', 'Color', [0.1 0.3 0.6], 'LineWidth', 2);
    
    % Plot target line (slope = 0.11 passing through means)
    target_slope = 0.11;
    mean_w = mean(log_w);
    mean_z = mean(E_z);
    z_target = mean_z + target_slope * (w_grid - mean_w);
    p2 = plot(w_grid, z_target, ':', 'Color', [0.2 0.8 0.2], 'LineWidth', 2.5);
    
    % Add correlation and slope
    corr_matrix = corrcoef(log_w, E_z);
    corr_val = corr_matrix(1, 2);
    
    legend(['Model States (\nu = ' num2str(P.nu) ')'], ...
           sprintf('Model Fit (Slope = %.3f)', beta(2)), ...
           'Dvorkin Target (Slope = 0.110)', 'Location', 'best');
    
    xlabel('Log Wage');
    ylabel('Average Productivity (E(z))');
    title('Spatial Sorting (Model vs. Empirical Target)');
    
    grid on;
    box on;
    set(gca, 'FontSize', 12);
    
    % Save figure
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    saveas(fig, fullfile(out_dir, 'spatial_sorting.png'));
    close(fig);
    fprintf('Spatial sorting plot saved to %s\\n', fullfile(out_dir, 'spatial_sorting.png'));
end
