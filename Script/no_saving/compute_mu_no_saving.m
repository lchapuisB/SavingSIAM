function mu = compute_mu_no_saving(Xi, P, Pi_all)
% Computes the migration probabilities for the NO-SAVING counterfactual.
%
% mu_{jl}(z) is proportional to exp(-nu*psi_{jl} - nu*log(E_z'[exp(-gam*Xi_l(z'))]))

    J = P.J; Nz = P.Nz; nu = P.nu; g = P.gamma;

    % 1. Compute log( E_{z'|z,m}[exp(-gam*Xi_m(z'))] )
    log_Ez = zeros(J, Nz);
    for m = 1:J
        v    = -g * Xi(m, :);
        vmax = max(v);
        ev   = exp(v - vmax);
        log_Ez(m, :) = vmax + log(ev * Pi_all{m});
    end

    % 2. Compute migration probabilities
    mu = zeros(J, J, Nz);
    nu_B_ns = -nu * P.psi; 
    
    for iz = 1:Nz
        log_w = nu_B_ns - nu * log_Ez(:, iz)';  
        log_w = log_w - max(log_w, [], 2);
        w_exp = exp(log_w);
        mu(:, :, iz) = w_exp ./ sum(w_exp, 2);
    end
end
