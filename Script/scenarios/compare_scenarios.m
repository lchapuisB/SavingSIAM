function compare_scenarios(res_A, res_B, P)
    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_comparison');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    
    T = length(res_A.WY_t);
    t_vec = 1:T;
    
    % 2. Compare Aggregate Wealth
    figure('Visible', 'off');
    plot(t_vec, res_A.agg_wealth, 'LineWidth', 3.0, 'DisplayName', 'Perfect Foresight');
    hold on;
    plot(t_vec, res_B.agg_wealth, 'LineWidth', 3.0, 'DisplayName', 'Myopia');
    grid on;
    xlabel('Time');
    ylabel('Aggregate Wealth');
    title('Comparison of Aggregate Wealth');
    legend('Location', 'best');
    saveas(gcf, fullfile(out_dir, 'comparison_agg_wealth.png'));
    
    % Find indices for Washington and Florida
    idx_wa = find(strcmpi(P.state_names, 'Washington'));
    idx_fl = find(strcmpi(P.state_names, 'Florida'));
    
    % 4. State Mean Welfare (Wealth) Comparison
    figure('Visible', 'off');
    hold on;
    
    nat_mean_A = mean(res_A.mean_wealth_t, 1);
    nat_mean_B = mean(res_B.mean_wealth_t, 1);
    
    plot(t_vec, nat_mean_A, 'b-', 'LineWidth', 3.0, 'DisplayName', 'National Mean (Perfect Foresight)');
    plot(t_vec, nat_mean_B, 'r-', 'LineWidth', 3.0, 'DisplayName', 'National Mean (Myopia)');
    
    if ~isempty(idx_wa)
        plot(t_vec, res_A.mean_wealth_t(idx_wa, :), 'b--', 'LineWidth', 3.0, 'DisplayName', 'Washington (Perfect Foresight)');
        plot(t_vec, res_B.mean_wealth_t(idx_wa, :), 'r--', 'LineWidth', 3.0, 'DisplayName', 'Washington (Myopia)');
    end
    
    if ~isempty(idx_fl)
        plot(t_vec, res_A.mean_wealth_t(idx_fl, :), 'b:', 'LineWidth', 3.0, 'DisplayName', 'Florida (Perfect Foresight)');
        plot(t_vec, res_B.mean_wealth_t(idx_fl, :), 'r:', 'LineWidth', 3.0, 'DisplayName', 'Florida (Myopia)');
    end
    
    grid on;
    xlabel('Time');
    ylabel('Mean Wealth (Welfare proxy)');
    title('Mean Wealth Comparison: Perfect Foresight vs Myopia');
    legend('Location', 'best');
    saveas(gcf, fullfile(out_dir, 'comparison_welfare_wealth.png'));

       
    % Compute CE_mean for res_A (Perfect Foresight)
    kappa = P.r / (1 + P.r);
    CE_mean_A = zeros(P.J, T);
    for t = 1:T
        for l = 1:P.J
            pop = sum(res_A.Lambda_trans(l, :, t));
            if pop > 1e-15
                weights = res_A.Lambda_trans(l, :, t) / pop;
                mean_xi = sum(res_A.Xi_trans(l, :, t) .* weights);
                % CE = kappa * a_bar + Xi
                CE_mean_A(l, t) = kappa * res_A.mean_wealth_t(l, t) + mean_xi;
            end
        end
    end
    % Compute CE_mean for res_B (Myopia)
    CE_mean_B = zeros(P.J, T);
    for t = 1:T
        for l = 1:P.J
            pop = sum(res_B.Lambda_trans(l, :, t));
            if pop > 1e-15
                weights = res_B.Lambda_trans(l, :, t) / pop;
                mean_xi = sum(res_B.Xi_trans(l, :, t) .* weights);
                CE_mean_B(l, t) = kappa * res_B.mean_wealth_t(l, t) + mean_xi;
            end
        end
    end
    % 6. True Welfare Comparison (Certainty Equivalent)
    figure('Visible', 'off');
    hold on;
    
    nat_CE_A = mean(CE_mean_A, 1);
    nat_CE_B = mean(CE_mean_B, 1);
    
    plot(t_vec, nat_CE_A, 'b-', 'LineWidth', 3.0, 'DisplayName', 'National Mean (Perfect Foresight)');
    plot(t_vec, nat_CE_B, 'r-', 'LineWidth', 3.0, 'DisplayName', 'National Mean (Myopia)');
    
    if ~isempty(idx_wa)
        plot(t_vec, CE_mean_A(idx_wa, :), 'b--', 'LineWidth', 3.0, 'DisplayName', 'Washington (Perfect Foresight)');
        plot(t_vec, CE_mean_B(idx_wa, :), 'r--', 'LineWidth', 3.0, 'DisplayName', 'Washington (Myopia)');
    end
    
    if ~isempty(idx_fl)
        plot(t_vec, CE_mean_A(idx_fl, :), 'b:', 'LineWidth', 3.0, 'DisplayName', 'Florida (Perfect Foresight)');
        plot(t_vec, CE_mean_B(idx_fl, :), 'r:', 'LineWidth', 3.0, 'DisplayName', 'Florida (Myopia)');
    end
    
    grid on;
    xlabel('Time');
    ylabel('Mean Certainty Equivalent (Welfare)');
    title('Welfare (CE) Comparison: Perfect Foresight vs Myopia');
    legend('Location', 'best');
    saveas(gcf, fullfile(out_dir, 'comparison_true_welfare.png'));
    % 7. Welfare Difference between Washington and Florida
    if ~isempty(idx_wa) && ~isempty(idx_fl)
        figure('Visible', 'off');
        hold on;
        
        diff_CE_A = CE_mean_A(idx_wa, :) - CE_mean_A(idx_fl, :);
        diff_CE_B = CE_mean_B(idx_wa, :) - CE_mean_B(idx_fl, :);
        
        plot(t_vec, diff_CE_A, 'b-', 'LineWidth', 3.0, 'DisplayName', 'Diff (WA - FL) Perfect Foresight');
        plot(t_vec, diff_CE_B, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Diff (WA - FL) Myopia');
        
        grid on;
        xlabel('Time');
        ylabel('Welfare Difference (WA - FL)');
        title('Welfare Gap between Washington and Florida');
        legend('Location', 'best');
        saveas(gcf, fullfile(out_dir, 'comparison_welfare_diff_WA_FL.png'));
    end
    fprintf('  Comparison plots saved to %s\n', out_dir);
end
    


