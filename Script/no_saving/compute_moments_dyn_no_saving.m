function [agg_earn, agg_wealth, agg_cons, agg_mobility, mean_earn_t, mean_wealth_t, WY_t, mean_ce_t, Flow_jl, Z_mean_jl] = compute_moments_dyn_no_saving(Lambda, Xi, paths, P)
% Dynamic moments for the NO-SAVING counterfactual.
% Wealth is identically zero. CE is simply the population-weighted Xi.

    J  = P.J;
    Nz = P.Nz;
    T  = size(Lambda, 3);
    delta = P.delta;
    iz0 = round(Nz / 2);

    agg_earn     = zeros(T, 1);
    agg_wealth   = zeros(T, 1);
    agg_cons     = zeros(T, 1);
    agg_mobility = zeros(T, 1);
    mean_earn_t  = zeros(J, T);
    mean_wealth_t= zeros(J, T);
    WY_t         = zeros(T, 1);
    mean_ce_t    = zeros(J, T);
    Flow_jl      = zeros(J, J, T);
    Z_mass_jl    = zeros(J, J, T);
    Z_mean_jl    = zeros(J, J, T);

    for t = 1:T-1
        P_t = P;
        P_t.w       = paths.w(:, t);
        P_t.w_tilde = paths.w_tilde(:, t);
        P_t.y_tilde = paths.y_tilde(:, t);
        P_t.sig2_eps= paths.sig2_eps(:, t);
        P_t.sig2_eta= paths.sig2_eta(:, t);
        P_t.A       = paths.A(:, t);
        P_t.q       = paths.q(:, t);

        Xi_t   = Xi(:, :, t);
        Pi_all = build_Pi_all(P_t);

        mu = compute_mu_no_saving(Xi_t, P_t, Pi_all);

        Lam_t = Lambda(:, :, t);              % J × Nz

        pop_l = sum(Lam_t, 2);
        pop_l_safe = max(pop_l, 1e-15);
        for j = 1:J
            mean_ce_t(j, t) = sum(Lam_t(j, :) .* Xi_t(j, :)) / pop_l_safe(j);
        end

        % ── aggregate mobility ──
        for j = 1:J
            stay_prob = 0;
            for iz = 1:Nz
                stay_prob = stay_prob + Lam_t(j, iz) * mu(j, j, iz);
            end
            agg_mobility(t) = agg_mobility(t) + (sum(Lam_t(j,:)) - stay_prob);
        end

        % ── flows & z-sorting ──
        for j = 1:J
            for ell = 1:J
                for iz = 1:Nz
                    mass = Lam_t(j, iz);
                    if mass < 1e-15; continue; end
                    f = mass * mu(j, ell, iz);
                    Flow_jl(j, ell, t) = Flow_jl(j, ell, t) + f;
                    Z_mass_jl(j, ell, t) = Z_mass_jl(j, ell, t) + f * P_t.z_grid(iz);
                end
            end
        end
    
        total_deaths = delta * pop_l;
        for j = 1:J
            for ell = 1:J
                f_nb = total_deaths(j) * mu(j, ell, iz0);
                Flow_jl(j, ell, t) = Flow_jl(j, ell, t) + f_nb;
                Z_mass_jl(j, ell, t) = Z_mass_jl(j, ell, t) + f_nb * P_t.z_grid(iz0);
            end
        end

        for j = 1:J
            for ell = 1:J
                if Flow_jl(j, ell, t) > 1e-15
                    Z_mean_jl(j, ell, t) = Z_mass_jl(j, ell, t) / Flow_jl(j, ell, t);
                end
            end
        end

        Lam_next = Lambda(:, :, t+1);
        pop_l_next = max(sum(Lam_next, 2), 1e-15);
        for ell = 1:J
            mean_earn_t(ell, t) = P_t.w(ell) * (1 + sum(P_t.z_grid(:) .* Lam_next(ell, :)') / pop_l_next(ell));
        end
        agg_earn(t) = sum(mean_earn_t(:, t) .* sum(Lam_next, 2));

        % Without saving, consumption equals income
        % The wealth logic is bypassed.
    end

    % ── terminal period ──
    agg_earn(T)     = agg_earn(T-1);
    agg_cons(T)     = agg_cons(T-1);
    agg_mobility(T) = agg_mobility(T-1);
    mean_earn_t(:,T) = mean_earn_t(:,T-1);

    Lam_T = Lambda(:,:,T);
    Xi_T = Xi(:,:,T);
    for j = 1:J
        mean_ce_t(j,T) = sum(Lam_T(j,:) .* Xi_T(j,:)) / max(sum(Lam_T(j,:)), 1e-15);
    end
    Flow_jl(:,:,T)   = Flow_jl(:,:,T-1);
    Z_mean_jl(:,:,T) = Z_mean_jl(:,:,T-1);
end
