function plot_model_vs_empirical_WY(P, wealth_to_earn, output_dir)
    set_aesthetics();
% PLOT_MODEL_VS_EMPIRICAL_WY  Scatter of model W/Y vs empirical W/Y
%   (Geowealth-US 2010) across states, with OLS regression line.
%
%   Inputs:
%       P              – parameter structure (needs P.state_names, P.J)
%       wealth_to_earn – [J x 1] model wealth-to-earnings ratio per state
%       output_dir     – path for saving the figure


    data_path = fullfile(fileparts(mfilename('fullpath')), ...
                         '..', '..', 'Data', 'state_wealth_inequality.csv');
    T = readtable(data_path);
    T2010 = T(T.year == 2010, :);

    %% 2. Compute empirical national average income for normalization
    % (In the model, population-weighted average wage is normalized to 1)
    national_income = sum(T2010.income_mean .* T2010.households) / sum(T2010.households);

    %% 3. Match states and compute empirical W/Y
    J = P.J;
    emp_WY = NaN(J, 1);
    for j = 1:J
        idx = strcmpi(T2010.state, P.state_names{j});
        if any(idx)
            % Normalize both wealth and income by national average income
            norm_W = T2010.wealth_mean(idx) / national_income;
            norm_Y = T2010.income_mean(idx) / national_income;
            emp_WY(j) = norm_W / norm_Y;
        end
    end

    valid = ~isnan(emp_WY);
    x = emp_WY(valid);              
    y = wealth_to_earn(valid);      
    names = P.state_names(valid);

    X_reg = [ones(numel(x), 1), x];
    b = X_reg \ y;
    x_fit = linspace(min(x), max(x), 200)';
    y_fit = b(1) + b(2) * x_fit;

    y_hat = X_reg * b;
    SS_res = sum((y - y_hat).^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = 1 - SS_res / SS_tot;

    %% 4. Plot
    fig = figure('Visible', 'off', 'Position', [100 100 700 550]);
    hold on;

    scatter(x, y, 50, [0.20 0.55 0.45], 'filled', ...
            'MarkerFaceAlpha', 0.7);

    plot(x_fit, y_fit, '-', ...
         'LineWidth', 3.0);

    % 45-degree line for reference
    lim_lo = min([min(x), min(y)]) * 0.9;
    lim_hi = max([max(x), max(y)]) * 1.1;
    plot([lim_lo lim_hi], [lim_lo lim_hi], '--', ...
         'Color', [0.6 0.6 0.6], 'LineWidth', 3.0);

    % Label a few states
    label_states = {'California', 'New York', 'Hawaii', ...
                    'Mississippi', 'West Virginia', 'District of Columbia', ...
                    'Connecticut', 'New Jersey'};
    for k = 1:numel(names)
        if ismember(names{k}, label_states)
            text(x(k), y(k), ['  ' names{k}], ...
                 'FontSize', 18);
        end
    end

    xlabel('Empirical W/Y (Geowealth-US 2010)', 'FontSize', 18);
    ylabel('Model W/Y', 'FontSize', 18);
    title(sprintf('Model vs Empirical Wealth-to-Earnings   (R^2 = %.3f)', R2), ...
          'FontSize', 18);
    legend({'States', sprintf('OLS fit  (slope = %.3f)', b(2)), '45° line'}, ...
           'Location', 'northwest', 'FontSize', 18);

    grid on; box on;
    set(gca, 'FontSize', 18);
    hold off;

    %% 5. Save
    save_fig(fig, output_dir, 'model_vs_empirical_WY.png');
    fprintf('Saved model_vs_empirical_WY.png\n');
end
