function Xi = compute_Xi_no_saving(P, tol, Xi_guess)
% Solves the fixed point Xi for the NO-SAVING counterfactual.
%
% Recursion:
% exp(-gam*Xi_ell(z')) = exp(F_ell(z')) + beta * Omega_ell(z')
% where
% F_ell(z') = -gam*y_til(ell) - gam*w_til(ell)*z' - A(ell) + 0.5*(gam*w_til(ell))^2 * sig2_eta(ell)
% Omega_ell(z') = lambda_bar * [ (1/J) sum_m exp(-nu*psi_{ell,m} - nu*log(E_z'[exp(-gam*Xi_m)])) ]^{-1/nu}

if nargin < 2; tol = 1e-7; end

J   = P.J;
Nz  = P.Nz;
nu  = P.nu;
gam = P.gamma;
A   = P.A(:); 
lambda_bar = P.lambda_bar;
beta = P.beta;
w_til = P.w_tilde(:);
y_til = P.y_tilde(:);
sig2_eta = P.sig2_eta(:);
psi = P.psi;
z_grid = P.z_grid(:)'; 

if nargin < 3 || isempty(Xi_guess)
    Xi_init = zeros(J, Nz);
else
    Xi_init = Xi_guess;
end

% Flow utility exponent pre-computation
F_exp = zeros(J, Nz);
for ell = 1:J
    F_exp(ell, :) = -gam * y_til(ell) - gam * w_til(ell) * z_grid -  A(ell) + 0.5 * (gam * w_til(ell))^2 * sig2_eta(ell);
end

Pi_trans = zeros(Nz, Nz, J);
for m = 1:J
    [~ , PI_temp] = discretize_normal(P.Nz, 0, P.rho, sqrt(P.sig2_eps(m)));
    Pi_trans(:, :, m) = PI_temp; 
end

tolerance = tol;
maxiter = 2000;
iter = 0;
error_val = inf;
Xi = Xi_init;

% Transition cost matrix inside Weibull
nu_B_ns = -nu * psi; 

while iter < maxiter && error_val > tolerance
    Xi_old = Xi;
    Xi_new = zeros(J, Nz);
    log_Ez = zeros(J, Nz);

    % 1. Compute log( E_{z'|z,m}[exp(-gam*Xi_old_m(z'))] )
    for m = 1:J
        V_m = -gam * Xi_old(m, :);
        V_max = max(V_m);
        E_val_vec = Pi_trans(:, :, m) * exp(V_m - V_max)'; 
        log_Ez(m, :) = V_max + log(E_val_vec');
    end

    % 2. Compute the new Xi for each origin ell and past shock iz
    for iz = 1:Nz
        % arg_exp is JxJ: row ell (origin), col m (destination)
        arg_exp = nu_B_ns - nu * log_Ez(:, iz)'; 
        
        % Log-sum-exp over destinations m
        max_arg = max(arg_exp, [], 2);
        mean_exp = mean(exp(arg_exp - max_arg), 2);
        log_sum_m = max_arg + log(mean_exp);
        
        % Log of discounted continuation value: log(beta * Omega_ell(z'))
        log_beta_Omega = log(beta) + log(lambda_bar) - (1/nu) * log_sum_m; 
        
        % 3. Add flow payoff and continuation value (both in log-space)
        flow_val = F_exp(:, iz);
        cont_val = log_beta_Omega; 
        
        % log-sum-exp trick for (flow + continuation)
        max_val = max(flow_val, cont_val);
        log_sum_terms = max_val + log(exp(flow_val - max_val) + exp(cont_val - max_val));
        
        % Final value
        Xi_new(:, iz) = -(1/gam) * log_sum_terms;
    end

    Xi = Xi_new;
    error_val = max(abs(Xi(:) - Xi_old(:)));
    iter = iter + 1;
end
end
