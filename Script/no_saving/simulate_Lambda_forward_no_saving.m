function Lambda = simulate_Lambda_forward_no_saving(P_2010, Lambda_SS_2010, Xi_trans, paths, T)
    J   = P_2010.J;
    Nz  = P_2010.Nz;
    N_states = J * Nz;

    Lambda = zeros(J, Nz, T);
    Lambda(:, :, 1) = Lambda_SS_2010;

    fprintf('  Forward simulation (No-Saving)...\n');
    for t = 1:(T-1)
        P_t = P_2010;
        P_t.w       = paths.w(:, t);
        P_t.w_tilde = paths.w_tilde(:, t);
        P_t.y_tilde = paths.y_tilde(:, t);
        P_t.sig2_eps= paths.sig2_eps(:, t);
        P_t.sig2_eta= paths.sig2_eta(:, t);
        P_t.A       = paths.A(:, t);
        P_t.q       = paths.q(:, t);

        Pi_all = build_Pi_all(P_t);
        Xi_t   = Xi_trans(:, :, t);

        mu = compute_mu_no_saving(Xi_t, P_t, Pi_all);
        T_full = build_T_t(mu, Pi_all, P_t);

        Lam_vec = reshape(Lambda(:, :, t)', 1, N_states);  
        Lam_new = Lam_vec * T_full;                         

        tot = sum(Lam_new);
        if tot > 0; Lam_new = Lam_new / tot; end
        Lambda(:, :, t+1) = reshape(Lam_new, Nz, J)';

        if mod(t, 20) == 0
            fprintf('    t=%d: sum(Lambda)=%.6f\n', t, sum(Lam_new));
        end
    end
end
