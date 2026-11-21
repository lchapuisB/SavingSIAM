function [delta_CE_mean, delta_CE_var, delta_W_var] = plot_welfare_variance(P_SS1, P_SS2, Lambda_SS1, Lambda_SS2, a_bar_SS1, a_bar_SS2, Xi_SS1, Xi_SS2, output_dir)
    set_aesthetics();

    J = P_SS1.J;
    r = P_SS1.r;
    kappa = r / (1 + r);

    % 1. Initialization of State-Level Metrics
    CE_mean_SS1 = zeros(J, 1); CE_mean_SS2 = zeros(J, 1);
    CE_var_SS1  = zeros(J, 1); CE_var_SS2  = zeros(J, 1);
    W_var_SS1   = zeros(J, 1); W_var_SS2   = zeros(J, 1);

    % 2. Computation of Welfare and Variance
    for l = 1:J
        pop_SS1 = sum(Lambda_SS1(l, :));
        pop_SS2 = sum(Lambda_SS2(l, :));

        if pop_SS1 > 1e-15
            weights_1 = Lambda_SS1(l, :) / pop_SS1;
            CE_1 = kappa * a_bar_SS1(l, :) + Xi_SS1(l, :);
            
            CE_mean_SS1(l) = sum(CE_1 .* weights_1);
            CE_var_SS1(l)  = sum((CE_1 - CE_mean_SS1(l)).^2 .* weights_1);
            
            mean_a_1 = sum(a_bar_SS1(l, :) .* weights_1);
            W_var_SS1(l) = sum((a_bar_SS1(l, :) - mean_a_1).^2 .* weights_1);
        end

        if pop_SS2 > 1e-15
            weights_2 = Lambda_SS2(l, :) / pop_SS2;
            CE_2 = kappa * a_bar_SS2(l, :) + Xi_SS2(l, :);
            
            CE_mean_SS2(l) = sum(CE_2 .* weights_2);
            CE_var_SS2(l)  = sum((CE_2 - CE_mean_SS2(l)).^2 .* weights_2);
            
            mean_a_2 = sum(a_bar_SS2(l, :) .* weights_2);
            W_var_SS2(l) = sum((a_bar_SS2(l, :) - mean_a_2).^2 .* weights_2);
        end
    end

    % 3. Overall Macro Welfare
    agg_CE_SS1 = sum(CE_mean_SS1 .* sum(Lambda_SS1, 2));
    agg_CE_SS2 = sum(CE_mean_SS2 .* sum(Lambda_SS2, 2));
    fprintf('\n=== Welfare Analysis ===\n');
    fprintf('Aggregate Certainty Equivalent (SS1): %.4f\n', agg_CE_SS1);
    fprintf('Aggregate Certainty Equivalent (SS2): %.4f\n', agg_CE_SS2);
    fprintf('Welfare Variation (Aggregate):         %.4f%%\n', ((agg_CE_SS2 - agg_CE_SS1)/abs(agg_CE_SS1))*100);

    % 4. State-Level Deltas
    delta_CE_mean = CE_mean_SS2 - CE_mean_SS1;
    CE_std_SS1 = sqrt(CE_var_SS1); CE_std_SS2 = sqrt(CE_var_SS2);
    W_std_SS1 = sqrt(W_var_SS1); W_std_SS2 = sqrt(W_var_SS2);
    
    delta_CE_var = zeros(J, 1);
    delta_W_var = zeros(J, 1);
    
    valid_CE = CE_std_SS1 > 1e-8;
    delta_CE_var(valid_CE) = (CE_std_SS2(valid_CE) - CE_std_SS1(valid_CE)) ./ CE_std_SS1(valid_CE) * 100;
    
    valid_W = W_std_SS1 > 1e-8;
    delta_W_var(valid_W) = (W_std_SS2(valid_W) - W_std_SS1(valid_W)) ./ W_std_SS1(valid_W) * 100;

    % Plotting moved to plot_inequality_maps.m
end