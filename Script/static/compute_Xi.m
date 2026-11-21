function Xi = compute_Xi(P, tol, Xi_guess)
% Solves the fixed point Xi, using eq(9)
%
%  - (gam / r_tilde) * Xi_ell(z) = 
%                        D_{1,ell} + D_{2,ell}(z)
%                        -(1/nu)*log(\mean_{m =1}^J exp(nu*B_ellm
%                        - nu*log(E_{z'|z,ell}[exp(-gamma*Xi_m(z'))]))
%                       
% with D_{1,ell} = log(beta_tilde) + A_ell + log(lambda_bar)
%      D_{2,ell}(z) = -gamma*kappa*(1+r)*w_tilde_ell*z     
%      B_{ell,m} = A_m - psi{ell,m} + gamma*kappa*y_tilde_m 
%                     - 0.5*(kappa*gamma*w_tilde_ell*sigma_{eta_ell})^2
%      r_tilde = 1/(1+r) and beta_tilde = beta*(1+r)
%
% 
% WWe use an iterative method to find a fixed point 

%% -- INITIALIZATION
if nargin < 2; tol = 1e-7; end

J   = P.J;
Nz  = P.Nz;
r   = P.r;
nu  = P.nu;
gam = P.gamma;
A   = P.A(:); 
lambda_bar = P.lambda_bar;
beta = P.beta;
w_til = P.w_tilde(:);
y_til = P.y_tilde(:);
sig2_eta = P.sig2_eta(:);
r_til = 1/(1+r);  
b_til = beta * (1 + r);
psi = P.psi;
kappa = P.kappa;
z_grid = P.z_grid(:)'; 

if nargin < 3 || isempty(Xi_guess)
    Xi_init = zeros(J, Nz);
else
    Xi_init = Xi_guess;
end

D1 = log(b_til) + A + log(lambda_bar);
D2 = -gam * kappa * (1 + r) * (w_til * z_grid); 

B1hat = gam * kappa * y_til - 0.5 * (kappa * gam * w_til).^2 .* sig2_eta;
B = (A' + B1hat') - psi; 
nu_B = nu * B; 

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

while iter < maxiter && error_val > tolerance
    Xi_old = Xi;
    Xi_new = zeros(J, Nz);
    log_Ez = zeros(J, Nz);

    for m = 1:J
        V_m = -gam * Xi_old(m, :);
        V_max = max(V_m);
        E_val_vec = Pi_trans(:, :, m) * exp(V_m - V_max)'; 
        log_Ez(m, :) = V_max + log(E_val_vec');
    end

    for iz = 1:Nz
        arg_exp = nu_B - nu * log_Ez(:, iz)'; 
        
        max_arg = max(arg_exp, [], 2);
        mean_exp = mean(exp(arg_exp - max_arg), 2);
        log_sum_m = max_arg + log(mean_exp);
        
        RHS = D1 + D2(:, iz) - (1/nu) * log_sum_m;
        Xi_new(:, iz) = -(r_til / gam) * RHS;
    end

    Xi = Xi_new;
    error_val = max(abs(Xi(:) - Xi_old(:)));
    iter = iter + 1;
end
end