function plot_did_comparison(P, T, T_shock, ...
    mean_ce_ctrl, mean_ce_treat, delta_ce_t, ...
    agg_CE_ctrl, agg_CE_treat, ...
    Lambda_ctrl, Lambda_treat, ...
    mean_wealth_ctrl, mean_wealth_treat, ...
    agg_mobility_ctrl, agg_mobility_treat, ...
    output_dir)
    
    set_aesthetics('subplot');
% Plots for the No-Migration DiD Counterfactual
% ============================================================

    t_vec = 1:T;
    blue  = [0.25 0.52 0.72];
    red   = [0.82 0.42 0.35];
    plum  = [0.62 0.40 0.68];

    %% ── FIGURE 1: Aggregate CE over time (control vs treatment) ──────────
    fig1 = figure('Visible', 'off', 'Position', [100 100 1400 500]);

    subplot(1, 2, 1);
    set(gcf, 'Tag', 'aggregate_certainty_equivalent');
    plot(t_vec, agg_CE_ctrl, 'Color', blue, 'LineWidth', 3.0, 'DisplayName', 'Control (ban, no CC)');
    hold on;
    plot(t_vec, agg_CE_treat, 'Color', red, 'LineWidth', 3.0, 'DisplayName', 'Treatment (ban + CC)');
    xline(T_shock, '--k', 'LineWidth', 3.0, 'DisplayName', 'End of CC ramp');
    grid on;
    xlabel('Time (year)'); ylabel('Aggregate Certainty Equivalent');
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    subplot(1, 2, 2);
    set(gcf, 'Tag', 'delta_ce_control_treatment');
    did_agg = agg_CE_treat - agg_CE_ctrl;
    plot(t_vec, did_agg, 'Color', plum, 'LineWidth', 3.0);
    hold on;
    xline(T_shock, '--k', 'LineWidth', 3.0);
    yline(0, ':k');
    grid on;
    xlabel('Time (year)'); ylabel('\Delta CE (Control - Treatment)');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig1, output_dir, 'did_aggregate_welfare.png');

    %% ── FIGURE 2: State-level CE at key horizons ─────────────────────────
    fig2 = figure('Visible', 'off', 'Position', [100 100 1600 500]);
    
    horizons = [round(T_shock/2), T_shock, min(2*T_shock, T)];
    horizon_labels = arrayfun(@(h) sprintf('t = %d', h), horizons, 'UniformOutput', false);

    for p = 1:3
        t_h = horizons(p);
        ax = subplot(1, 3, p);
        set(gcf, 'Tag', sprintf('%d', p));
        plot_choropleth(P.state_names, -delta_ce_t(:, t_h), ax);
    end
    save_fig(fig2, output_dir, 'did_state_maps.png');

    %% ── FIGURE 3: Selected states CE trajectories ────────────────────────
    fig3 = figure('Visible', 'off', 'Position', [100 100 1400 500]);

    highlight_states = {'Florida', 'California', 'Texas', 'New York', 'Washington', 'Wyoming'};
    colors = get(groot, 'defaultAxesColorOrder');

    subplot(1, 2, 1);
    set(gcf, 'Tag', 'ce_treat_ce_ctrl');
    hold on;
    for s = 1:length(highlight_states)
        idx = find(strcmpi(P.state_names, highlight_states{s}));
        if ~isempty(idx)
            c_idx = mod(s-1, size(colors, 1)) + 1;
            plot(t_vec, mean_ce_treat(idx, :) - mean_ce_ctrl(idx, :), ...
                'Color', colors(c_idx, :), 'LineWidth', 3.0, ...
                'DisplayName', highlight_states{s});
        end
    end
    xline(T_shock, '--k', 'LineWidth', 3.0);
    yline(0, ':k');
    grid on;
    xlabel('Time (year)'); ylabel('CE_{treat} - CE_{ctrl}');
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    subplot(1, 2, 2);
    set(gcf, 'Tag', 'delta_population_share');
    hold on;
    pop_ctrl  = squeeze(sum(Lambda_ctrl, 2));   % J x T
    pop_treat = squeeze(sum(Lambda_treat, 2));  % J x T
    for s = 1:length(highlight_states)
        idx = find(strcmpi(P.state_names, highlight_states{s}));
        if ~isempty(idx)
            c_idx = mod(s-1, size(colors, 1)) + 1;
            plot(t_vec, pop_treat(idx, :) - pop_ctrl(idx, :), ...
                'Color', colors(c_idx, :), 'LineWidth', 3.0, ...
                'DisplayName', highlight_states{s});
        end
    end
    xline(T_shock, '--k', 'LineWidth', 3.0);
    yline(0, ':k');
    grid on;
    xlabel('Time (year)'); ylabel('\Delta Population Share');
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig3, output_dir, 'did_state_trajectories.png');

    %% ── FIGURE 4: Wealth comparison ──────────────────────────────────────
    fig4 = figure('Visible', 'off', 'Position', [100 100 1200 400]);

    nat_wealth_ctrl  = sum(mean_wealth_ctrl .* pop_ctrl, 1);
    nat_wealth_treat = sum(mean_wealth_treat .* pop_treat, 1);

    subplot(1, 2, 1);
    set(gcf, 'Tag', 'aggregate_wealth');
    plot(t_vec, nat_wealth_ctrl, 'Color', blue, 'LineWidth', 3.0, 'DisplayName', 'Control');
    hold on;
    plot(t_vec, nat_wealth_treat, 'Color', red, 'LineWidth', 3.0, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'LineWidth', 3.0);
    grid on;
    xlabel('Time (year)'); ylabel('Aggregate Wealth');
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    subplot(1, 2, 2);
    set(gcf, 'Tag', 'mobility_rate');
    plot(t_vec, agg_mobility_ctrl, 'Color', blue, 'LineWidth', 3.0, 'DisplayName', 'Control');
    hold on;
    plot(t_vec, agg_mobility_treat, 'Color', red, 'LineWidth', 3.0, 'DisplayName', 'Treatment');
    xline(T_shock, '--k', 'LineWidth', 3.0);
    grid on;
    xlabel('Time (year)'); ylabel('Mobility Rate');
    legend('Location', 'best');
    yl = ylim; dy = diff(yl); if dy > 0, ylim([yl(1) - 0.06*dy, yl(2) + 0.06*dy]); end

    save_fig(fig4, output_dir, 'did_wealth_mobility.png');

    fprintf('  DiD plots saved to %s\n', output_dir);
end
