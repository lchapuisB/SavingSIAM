function [agg_earn, agg_wealth, agg_cons, agg_mobility, mean_earn_t, mean_wealth_t, WY_t, mean_ce_t, Flow_jl, Z_mean_jl] = compute_moments_dyn(Lambda, Xi, a_bar_init, paths, P)
    J  = P.J;
    Nz = P.Nz;
    T  = size(Lambda, 3);
    r  = P.r;
    kappa = r / (1 + r);
    delta = P.delta;
    N_states = J * Nz;
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

    a_bar = a_bar_init;  

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
        B1hat  = compute_B1hat(P_t);
        Pi_all = build_Pi_all(P_t);

        mu = compute_mu_all(Xi_t, B1hat, P_t, Pi_all);
        [~, T_surv] = build_T_t(mu, Pi_all, P_t);

        Lam_t = Lambda(:, :, t);              % J × Nz
        Lam_vec = reshape(Lam_t', 1, N_states);

        % ── aggregate wealth (before transition) ──
        a_bar_vec = reshape(a_bar', 1, N_states);
        agg_wealth(t) = Lam_vec * a_bar_vec';

        pop_l = sum(Lam_t, 2);
        pop_l_safe = max(pop_l, 1e-15);
        for j = 1:J
            mean_wealth_t(j, t) = sum(Lam_t(j, :) .* a_bar(j, :)) / pop_l_safe(j);
            mean_ce_t(j, t) = sum(Lam_t(j, :) .* (kappa * a_bar(j, :) + Xi_t(j, :))) / pop_l_safe(j);
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

        %for ell = 1:J
         %   c_bar_ell = sum(Xi_t(ell, :)' .* Lam_next(ell, :)') / pop_l_next(ell) ...
         %             + kappa * P_t.y_tilde(ell);
        %end

        C_income_vec = zeros(1, N_states);
        for ell = 1:J
            for iz = 1:Nz
                s = (ell - 1) * Nz + iz;
                C_income_vec(s) = P_t.y_tilde(ell) + (1+r) * (P_t.w_tilde(ell) * P_t.z_grid(iz) - Xi_t(ell, iz));
            end
        end
        wealth_mass_in = (a_bar_vec .* Lam_vec) * T_surv;  
        Lam_new_vec = reshape(Lam_next', 1, N_states);

  
        income_mass = Lam_new_vec .* C_income_vec;

        nb_wealth = zeros(1, N_states);
        if P.a0 ~= 0
            for j = 1:J
                for ell = 1:J
                    mass_nb = total_deaths(j) * mu(j, ell, iz0);
                    Pi_ell = Pi_all{ell};
                    for iz_new = 1:Nz
                        s = (ell-1)*Nz + iz_new;
                        nb_wealth(s) = nb_wealth(s) + mass_nb * Pi_ell(iz_new, iz0) * P.a0;
                    end
                end
            end
        end

        A_new_vec = (1+r)*(1-kappa) * wealth_mass_in + income_mass + nb_wealth;

        % recover per-capita wealth
        for ell = 1:J
            for iz = 1:Nz
                s = (ell-1)*Nz + iz;
                if Lam_new_vec(s) > 1e-15
                    a_bar(ell, iz) = A_new_vec(s) / Lam_new_vec(s);
                end
            end
        end

        WY_t(t) = agg_wealth(t) / (agg_earn(t) + 1e-30);
    end

    % ── terminal period ──
    agg_earn(T)     = agg_earn(T-1);
    agg_wealth(T)   = sum(a_bar(:,:) .* Lambda(:,:,T), 'all');
    agg_cons(T)     = agg_cons(T-1);
    agg_mobility(T) = agg_mobility(T-1);
    mean_earn_t(:,T) = mean_earn_t(:,T-1);

    Lam_T = Lambda(:,:,T);
    Xi_T = Xi(:,:,T);
    for j = 1:J
        mean_wealth_t(j,T) = sum(Lam_T(j,:) .* a_bar(j,:)) / max(sum(Lam_T(j,:)), 1e-15);
        mean_ce_t(j,T) = sum(Lam_T(j,:) .* (kappa * a_bar(j,:) + Xi_T(j,:))) / max(sum(Lam_T(j,:)), 1e-15);
    end
    WY_t(T) = agg_wealth(T) / (sum(agg_earn(T) .* sum(Lam_T,2)) + 1e-30);
    Flow_jl(:,:,T)   = Flow_jl(:,:,T-1);
    Z_mean_jl(:,:,T) = Z_mean_jl(:,:,T-1);
end