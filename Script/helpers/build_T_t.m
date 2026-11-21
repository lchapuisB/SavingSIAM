function [T_full, T_surv] = build_T_t(mu, Pi_all, P)
    J = P.J; Nz = P.Nz; delta = P.delta;
    N_states = J * Nz;
    iz0 = round(Nz / 2);

    idx = @(j, iz) (j - 1) * Nz + iz;

    % pre-bu
    Pi3 = zeros(Nz, Nz, J);
    for m = 1:J
        Pi3(:, :, m) = Pi_all{m}';  %
    end

    nnz_est = J * Nz * J * Nz;
    rows = zeros(nnz_est, 1);
    cols = zeros(nnz_est, 1);
    vals_s = zeros(nnz_est, 1);
    vals_d = zeros(nnz_est, 1);

    k = 0;
    for j = 1:J
        for iz = 1:Nz
            s = idx(j, iz);
            for ell = 1:J
                for iz_new = 1:Nz
                    k = k + 1;
                    s_prime = idx(ell, iz_new);
                    rows(k) = s;
                    cols(k) = s_prime;
                    vals_s(k) = (1-delta) * mu(j, ell, iz) * Pi3(iz, iz_new, ell);
                    vals_d(k) = delta     * mu(j, ell, iz0) * Pi3(iz0, iz_new, ell);
                end
            end
        end
    end

    T_surv = sparse(rows, cols, vals_s, N_states, N_states);
    T_death = sparse(rows, cols, vals_d, N_states, N_states);
    T_full = T_surv + T_death;
end