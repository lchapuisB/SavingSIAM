function Xi = solve_transition_myopia(P_base, Xi_initial, paths, T)
    Xi = zeros(P_base.J, P_base.Nz, T);
    
    fprintf('  Forward myopia sequence over %d periods...\n', T);
    Xi_guess = Xi_initial;
    
    for t = 1:T
        P_t = P_base;
        P_t.w         = paths.w(:, t);
        P_t.w_tilde   = paths.w_tilde(:, t);
        P_t.y_tilde   = paths.y_tilde(:, t);
        P_t.sig2_eps  = paths.sig2_eps(:, t);
        P_t.sig2_eta  = paths.sig2_eta(:, t);
        P_t.A         = paths.A(:, t);
        P_t.q         = paths.q(:,t);
        
        Xi_t = compute_Xi(P_t, 1e-6, Xi_guess);
        Xi(:, :, t) = Xi_t;
        
        Xi_guess = Xi_t; 
        
        if mod(t,20)==0 || t==1
            fprintf('    t=%3d: solved myopia steady state.\n', t);
        end
    end
end
