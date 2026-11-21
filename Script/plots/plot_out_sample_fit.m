function plot_out_sample_fit(P, Lambda, mean_wealth, Flow_jl, F_data, output_dir)
    set_aesthetics();
    
    fig = figure('Visible', 'off', 'Position', [100 100 1500 450]);
    
    %% --- Panel 1: W comparisons ---
    subplot(1, 3, 1);
    hold on;
    % Load empirical data (Geowealth-US, year 2010)
    data_path = fullfile(fileparts(mfilename('fullpath')), ...
                         '..', '..', 'Data', 'state_wealth_inequality.csv');
    T = readtable(data_path);
    T2010 = T(T.year == 2020, :);

    national_income = sum(T2010.income_mean .* T2010.households) / sum(T2010.households);

    J = P.J;
    emp_W = NaN(J, 1);
    for j = 1:J
        idx = strcmpi(T2010.state, P.state_names{j});
        if any(idx)
            emp_W(j) = T2010.wealth_mean(idx) / national_income;
        end
    end

    valid = ~isnan(emp_W);
    x_w = emp_W(valid);         
    y_w = mean_wealth(valid);    
    names_w = P.state_names(valid);

    X_reg = [ones(numel(x_w), 1), x_w];
    b_w = X_reg \ y_w;
    x_fit = linspace(min(x_w), max(x_w), 200)';
    y_fit = b_w(1) + b_w(2) * x_fit;

    y_hat = X_reg * b_w;
    SS_res = sum((y_w - y_hat).^2);
    SS_tot = sum((y_w - mean(y_w)).^2);
    R2_w = 1 - SS_res / SS_tot;

    scatter(x_w, y_w, 150, [0.20 0.40 0.70], 'filled', 'MarkerFaceAlpha', 0.7);
    plot(x_fit, y_fit, '-', 'LineWidth', 3.0, 'Color', [0.20 0.40 0.70]);

    label_states = {'California', 'New York', 'Mississippi', 'District of Columbia'};
    for k = 1:numel(names_w)
        if ismember(names_w{k}, label_states)
            text(x_w(k), y_w(k), ['  ' names_w{k}], 'FontSize', 16);
        end
    end

    xlabel('Empirical Mean Wealth', 'FontSize', 18);
    ylabel('Model Mean Wealth', 'FontSize', 18);
    grid on; box on;
    set(gca, 'FontSize', 18);
    hold off;

    %% --- Panel 2: Migration Sorting ---
    subplot(1, 3, 2);
    hold on;
    
    z_grid = P.z_grid(:)';
    E_z = zeros(J, 1);
    for j = 1:J
        mass_j = sum(Lambda(j, :));
        if mass_j > 1e-10
            E_z(j) = sum(Lambda(j, :) .* z_grid) / mass_j;
        end
    end
    
    log_w = log(P.w(:));
    
    X_z = [ones(J, 1), log_w];
    beta_z = (X_z'*X_z) \ (X_z'*E_z);
    
    scatter(log_w, E_z, 150, [0.7 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.6);
    
    w_min = min(log_w);
    w_max = max(log_w);
    w_grid = linspace(w_min, w_max, 100)';
    z_fit = beta_z(1) + beta_z(2) * w_grid;
    plot(w_grid, z_fit, '--', 'Color', [0.1 0.3 0.6], 'LineWidth', 2);
    
    target_slope = 0.11;
    mean_w = mean(log_w);
    mean_z = mean(E_z);
    z_target = mean_z + target_slope * (w_grid - mean_w);
    plot(w_grid, z_target, ':', 'Color', [0.2 0.8 0.2], 'LineWidth', 2.5);
    
    legend({'States', sprintf('Fit (Slope = %.3f)', beta_z(2)), 'Target (Slope = 0.110)'}, 'Location', 'best', 'FontSize', 16);
    xlabel('Log Wage', 'FontSize', 18);
    ylabel('E(z)', 'FontSize', 18);
    grid on; box on;
    set(gca, 'FontSize', 18);
    hold off;

    %% --- Panel 3: Migration Flows ---
    subplot(1, 3, 3);
    hold on;
    
    Flow_movers_model = Flow_jl;
    Flow_movers_model(logical(eye(J))) = 0;
    
    F_movers_emp = F_data;
    F_movers_emp(logical(eye(J))) = 0;
    
    emp_shares = F_movers_emp(:) / sum(F_movers_emp(:));
    mod_shares = Flow_movers_model(:) / sum(Flow_movers_model(:));
    
    valid_f = (emp_shares > 0) & (mod_shares > 0);
    x_f = log10(emp_shares(valid_f));
    y_f = log10(mod_shares(valid_f));
    
    scatter(x_f, y_f, 90, [0.4 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.5);
    
    min_val = min([x_f; y_f]);
    max_val = max([x_f; y_f]);
    plot([min_val, max_val], [min_val, max_val], '--k', 'LineWidth', 2);
    
    X_f = [ones(numel(x_f), 1), x_f];
    b_f = X_f \ y_f;
    x_fit_f = linspace(min_val, max_val, 100)';
    y_fit_f = b_f(1) + b_f(2) * x_fit_f;
    plot(x_fit_f, y_fit_f, '-', 'Color', [0.4 0.6 0.2], 'LineWidth', 2);
    
    y_hat_f = X_f * b_f;
    SS_res_f = sum((y_f - y_hat_f).^2);
    SS_tot_f = sum((y_f - mean(y_f)).^2);
    R2_f = 1 - SS_res_f / SS_tot_f;
    
    legend({'Movers', '45 Degree Line', sprintf('Fit (R^2 = %.2f)', R2_f)}, 'Location', 'best', 'FontSize', 16);
    xlabel('Log_{10}(Empirical Flow Share)', 'FontSize', 18);
    ylabel('Log_{10}(Model Flow Share)', 'FontSize', 18);
    grid on; box on;
    set(gca, 'FontSize', 18);
    hold off;

    %% --- Save ---
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save_fig(fig, output_dir, 'out_sample_fit.png');
    fprintf('Saved out_sample_fit.png to %s\n', output_dir);
end
