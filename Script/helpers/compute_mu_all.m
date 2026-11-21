function mu = compute_mu_all(Xi, B1hat, P, Pi_all)
    J = P.J; Nz = P.Nz; nu = P.nu; g = P.gamma;

    log_Ez = zeros(J, Nz);
    for m = 1:J
        v    = -g * Xi(m, :);
        vmax = max(v);
        ev   = exp(v - vmax);
        log_Ez(m, :) = vmax + log(ev * Pi_all{m});
    end
    B1_mat = (P.A(:)' + B1hat(:)') - P.psi;      

    mu = zeros(J, J, Nz);
    for iz = 1:Nz
        log_w = nu * B1_mat - nu * log_Ez(:, iz)';  
        log_w = log_w - max(log_w, [], 2);
        w_exp = exp(log_w);
        mu(:, :, iz) = w_exp ./ sum(w_exp, 2);
    end
end
