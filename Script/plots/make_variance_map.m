script_dir = fullfile(pwd, 'Script');
addpath(fullfile(script_dir, 'scenarios'), fullfile(script_dir, 'static'), fullfile(script_dir, 'helpers'), fullfile(script_dir, 'plots'), fullfile(script_dir, 'dynamics'), fullfile(script_dir, 'no_saving'));

model_static;
state_names = P.state_names;

wage_damages_path = fullfile(pwd, 'Output', 'data', 'shock_wage_2080.csv');
opts = detectImportOptions(wage_damages_path);
damages_wage = readtable(wage_damages_path, opts);
shock_wage = damages_wage.shock_wage;

delta_sigma2 = 50 * abs(shock_wage) * 100;

fig_var = figure('Visible', 'off', 'Position', [100 100 600 450]);
ax1 = axes(fig_var);
plot_choropleth(state_names, delta_sigma2, ax1);

out_dir = fullfile(pwd, 'Output', 'graph_uncertainty_variance_shocks');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

save_fig(fig_var, out_dir, 'variance_map.png');
fprintf('Variance map generated.\n');
