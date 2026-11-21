function plot_model_vs_empirical_W(P, mean_wealth, output_dir)
    set_aesthetics();
% PLOT_MODEL_VS_EMPIRICAL_W  Scatter of model mean wealth vs empirical mean
%   wealth (Geowealth-US 2010) across states, with OLS regression line.
%
%   Inputs:
%       P            – parameter structure (needs P.state_names, P.J)
%       mean_wealth  – [J x 1] model mean wealth per state
%       output_dir   – path for saving the figure

    %% 1. Load empirical data (Geowealth-US, year 2010)
    data_path = fullfile(fileparts(mfilename('fullpath')), ...
                         '..', '..', 'Data', 'state_wealth_inequality.csv');
    T = readtable(data_path);
    T2010 = T(T.year == 2020, :);

    national_income = sum(T2010.income_mean .* T2010.households) / sum(T2010.households);

    %% 3. Match states
    J = P.J;
    emp_W = NaN(J, 1);
    for j = 1:J
        idx = strcmpi(T2010.state, P.state_names{j});
        if any(idx)
            emp_W(j) = T2010.wealth_mean(idx) / national_income;
        end
    end

    valid = ~isnan(emp_W);
    x = emp_W(valid);         
    y = mean_wealth(valid);    
    names = P.state_names(valid);

    %% 3. OLS regression
    X_reg = [ones(numel(x), 1), x];
    b = X_reg \ y;
    x_fit = linspace(min(x), max(x), 200)';
    y_fit = b(1) + b(2) * x_fit;

    % R-squared
    y_hat = X_reg * b;
    SS_res = sum((y - y_hat).^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = 1 - SS_res / SS_tot;

    %% 4. Plot
    fig = figure('Visible', 'off', 'Position', [100 100 700 550]);
    hold on;

    scatter(x, y, 50, [0.20 0.40 0.70], 'filled', ...
            'MarkerFaceAlpha', 0.7);

    plot(x_fit, y_fit, '-', ...
         'LineWidth', 3.0);

    % Label a few states for readability
    label_states = {'California', 'New York', 'Hawaii', ...
                    'Mississippi', 'West Virginia', 'District of Columbia', ...
                    'Connecticut', 'New Jersey'};
    for k = 1:numel(names)
        if ismember(names{k}, label_states)
            text(x(k), y(k), ['  ' names{k}], ...
                 'FontSize', 18);
        end
    end

    xlabel('Empirical Mean Wealth (Normalized, Geowealth-US 2010)', 'FontSize', 18);
    ylabel('Model Mean Wealth (model units)', 'FontSize', 18);
    title(sprintf('Model vs Empirical Mean Wealth   (R^2 = %.3f)', R2), ...
          'FontSize', 18);
    legend({'States', sprintf('OLS fit  (slope = %.2f)', b(2))}, ...
           'Location', 'northwest', 'FontSize', 18);

    grid on; box on;
    set(gca, 'FontSize', 18);
    hold off;

    %% 5. Save
    save_fig(fig, output_dir, 'model_vs_empirical_W.png');
    fprintf('Saved model_vs_empirical_W.png\n');
end
