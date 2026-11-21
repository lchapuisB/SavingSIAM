function Xi = solve_transition_backward_no_saving(P_base, Xi_terminal, paths, T)
    Xi = zeros(P_base.J, P_base.Nz, T);
    Xi(:, :, T) = Xi_terminal;

    fprintf('  Backward induction (No-Saving) over %d periods...\n', T);
    for t = (T-1):-1:1

        P_t = P_base;
        P_t.w         = paths.w(:, t);
        P_t.w_tilde   = paths.w_tilde(:, t);
        P_t.y_tilde   = paths.y_tilde(:, t);
        P_t.sig2_eps  = paths.sig2_eps(:, t);
        P_t.sig2_eta  = paths.sig2_eta(:, t);
        P_t.A         = paths.A(:, t);
        P_t.q         = paths.q(:,t);

        P_next = P_base;
        P_next.w         = paths.w(:, t+1);
        P_next.w_tilde   = paths.w_tilde(:, t+1);
        P_next.y_tilde   = paths.y_tilde(:, t+1);
        P_next.sig2_eps  = paths.sig2_eps(:, t+1);
        P_next.sig2_eta  = paths.sig2_eta(:, t+1);
        P_next.A         = paths.A(:, t+1);
        P_next.q         = paths.q(:,t+1);

        Pi_all_next = build_Pi_all(P_next);

        Xi(:, :, t) = backward_step_Xi_no_saving(Xi(:, :, t+1), P_t, P_next, Pi_all_next);

        if mod(t,20)==0 || t==1
            dxi = max(abs(Xi(:,:,t) - Xi(:,:,t+1)), [], 'all');
            fprintf('    t=%3d: max|Xi_t - Xi_{t+1}| = %.6f\n', t, dxi);
        end
    end
end
