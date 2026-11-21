function P = calibrate_param_nu(P_init, targets, options)

    if nargin < 3 || isempty(options)
        options.inner_tol = 1e-6;
        options.Display   = 'iter';
    end
    F_data = readmatrix('/Users/leon/Documents/SavingSIAM/Output/data/empirical_F_matrix.csv');
    dist_mat = state_distances(P_init); 
    F_movers = F_data;
    F_movers(logical(eye(size(F_movers)))) = 0;
    total_empirical_movers = sum(F_movers(:));
    log_dist_weighted = F_movers .* log(1 + dist_mat);
    target_avg_log_dist = sum(log_dist_weighted(:)) / total_empirical_movers;

    if nargin < 2 || isempty(targets)
        targets.mob_rate  = 0.0250;
        targets.WY_target = 4.0900;
        targets.avg_log_dist = target_avg_log_dist; 
        targets.sorting_slope = 0.110; 
        targets.pop_share = P_init.pop_share;
        targets.pop_share = targets.pop_share(:) / sum(targets.pop_share(:));
        targets.young_wealth = [];
    end

    if ~isfield(P_init, 'dist_mat')
        P_init.dist_mat = state_distances(P_init);
    end
    P = P_init;
    
    % On optimise les 4 paramètres conjointement
    theta_0 = [0.8; 0.4; log(P.lambda_bar); log(P.nu)];

    fsolve_options = optimoptions('fsolve', ...
        'Algorithm', 'trust-region-dogleg', ...
        'Display', options.Display, ...
        'StepTolerance', 1e-8, ...          
        'FunctionTolerance', 1e-8, ...     
        'OptimalityTolerance', 1e-7, ...    
        'MaxIterations', 500, ...        
        'MaxFunctionEvaluations', 2000, ... 
        'FiniteDifferenceType', 'central', ...
        'UseParallel', false); 

    fprintf('\n============================================================\n');
    fprintf('CALIBRATION JOINTE SIMULTANEE (Les 4 moments empiriques)\n');
    fprintf('============================================================\n');
    
    obj_fun = @(theta) macro_moment_errors(theta, P, targets, options.inner_tol);
    
    [theta_sol, ~] = fsolve(obj_fun, theta_0, fsolve_options);

    psi_fixed_sol = exp(theta_sol(1));
    psi_dist_sol  = exp(theta_sol(2));
    P.lambda_bar  = exp(theta_sol(3));
    P.nu          = exp(theta_sol(4));

    P.psi_fixed = psi_fixed_sol;
    P.psi_dist = psi_dist_sol;

    P.psi = zeros(P.J, P.J);
    for j = 1:P.J
        for l = 1:P.J
            if j ~= l
                P.psi(j, l) = psi_fixed_sol + psi_dist_sol * log(1 + P.dist_mat(j, l));
            end
        end
    end
    
    P.A_init_done = false;
    [~, P] = macro_moment_errors(theta_sol, P, targets, options.inner_tol);

    delta = P.delta; r = P.r;
    Xi_curr = compute_Xi(P, 1e-7);
    E_y   = mean(P.y_tilde);
    E_Xi  = mean(Xi_curr(:));

    if isfield(targets, 'young_wealth') && ~isempty(targets.young_wealth)
        P.a0 = max(targets.young_wealth - (1-delta)/delta * (E_y - (1+r)*E_Xi), 0);
    end
     
fprintf('\n============================================================\n');
fprintf('STRUCTURAL PARAMETER CALIBRATION SUMMARY\n');
fprintf('============================================================\n');
fprintf('Outer Loop Frictions (Gravity & Scale):\n');
fprintf('  psi_fixed   (Fixed Relocation Cost) : %8.4f\n', P.psi_fixed);
fprintf('  psi_dist    (Distance Elasticity)   : %8.4f\n', P.psi_dist);
fprintf('  lambda_bar  (Weibull Scale)         : %8.4f\n', P.lambda_bar);
fprintf('  nu          (Weibull Shape)         : %8.4f\n', P.nu);
fprintf('\nInner Loop Clearings (Spatial Amenities A_l):\n');
fprintf('  Cross-Sectional Mean Amenity        : %8.4f\n', mean(P.A));
fprintf('  Cross-Sectional Variance            : %8.4f\n', var(P.A));

[min_A, min_idx] = min(P.A);
[max_A, max_idx] = max(P.A);
fprintf('  Minimum Amenity                     : %8.4f  (%s)\n', min_A, P.state_names{min_idx});
fprintf('  Maximum Amenity                     : %8.4f  (%s)\n', max_A, P.state_names{max_idx});

fprintf('\nDiagnostic State Audit:\n');
audit_states = {'California', 'New York', 'Wyoming', 'Texas'};
for i = 1:length(audit_states)
    state_idx = strcmp(P.state_names, audit_states{i});
    fprintf('  %-20s                : %8.4f\n', audit_states{i}, P.A(state_idx));
end
fprintf('============================================================\n');

end

function [err,P] = macro_moment_errors(theta, P, targets, tol)
    psi_fixed = exp(theta(1));
    psi_dist  = exp(theta(2));
    P.lambda_bar = exp(theta(3)); 
    P.nu = exp(theta(4));
    J = P.J;
    P.psi = zeros(J, J);
    
    for j = 1:J
        for l = 1:J
            if j ~= l
                P.psi(j, l) = psi_fixed + psi_dist * log(1 + P.dist_mat(j, l));
            end
        end
    end
 
    % --- CORRECTION CRITIQUE  ---
    P.A_init_done = false; 
    
    P = calibrate_amenities_blp(P, targets.pop_share, tol);
    Xi = compute_Xi(P, tol);
    [Flow_jl, ~, out_migration_rate, mean_wealth, mean_earn, ~, Lambda, ~] = compute_SS_moments(P, Xi);
    
    pop_model = sum(Lambda, 2);
    pop_model = pop_model / sum(pop_model);
    
    mob_model = sum(out_migration_rate .* pop_model);
    WY_model = sum(mean_wealth .* pop_model) / sum(mean_earn .* pop_model);
    
    Flow_movers = Flow_jl;
    Flow_movers(logical(eye(J))) = 0;
    total_model_movers = sum(Flow_movers(:));
    model_log_dist_weighted = Flow_movers .* log(1 + P.dist_mat);
    avg_log_dist_model = sum(model_log_dist_weighted(:)) / total_model_movers;
    
    % Calcul du tri spatial
    z_grid = P.z_grid(:)';
    E_z = zeros(J, 1);
    for j = 1:J
        mass_j = sum(Lambda(j, :));
        if mass_j > 1e-10
            E_z(j) = sum(Lambda(j, :) .* z_grid) / mass_j;
        end
    end
    log_w = log(P.w(:));
    X = [ones(J, 1), log_w];
    beta = (X'*X) \ (X'*E_z);
    slope_model = beta(2);

    err = zeros(4, 1);
    err(1) = (mob_model - targets.mob_rate) / targets.mob_rate;
    err(2) = (avg_log_dist_model - targets.avg_log_dist) / targets.avg_log_dist;
    err(3) = (WY_model - targets.WY_target) / targets.WY_target;
    err(4) = (slope_model - targets.sorting_slope) / targets.sorting_slope;
end

function P = calibrate_amenities_blp(P, pop_target, tol)
    log_pop_target = log(pop_target);
    
    if ~isfield(P,'A_init_done') || ~P.A_init_done
        kap   = P.r / (1 + P.r);
        w_til = P.w;
        y_til = w_til - P.q;
        
        B1hat = P.gamma*kap*y_til - 0.5*(kap*P.gamma*w_til).^2 .* P.sig2_eta;
        
        P.A = (1/P.nu) * log_pop_target - B1hat;
        P.A = P.A - mean(P.A);
        P.A_init_done = true;
    end
    
    err = inf;
    it = 0;
    if ~isfield(P, 'Xi_warm')
        P.Xi_warm = [];
    end

    while err > tol && it < 100
        it = it + 1;
        
        Xi_curr = compute_Xi(P, 1e-7, P.Xi_warm);
        P.Xi_warm = Xi_curr;
        
        [~,~,~,~,~,~, Lam, ~] = compute_SS_moments(P, Xi_curr);
        
        L_raw   = sum(Lam, 2);
        L_model = L_raw / sum(L_raw);
      
        pop_model = L_model(:);
        
        resid = log_pop_target - log(max(pop_model, 1e-15));
        err = max(abs(resid));
        
        if err < tol
            break;
        end
        
        P.A = P.A + (1/P.nu) * resid;
        P.A = P.A - mean(P.A);
    end
end