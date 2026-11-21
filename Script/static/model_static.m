% ============================================================
% Static model 2010 : Dvorkin & Greaney replication
% ============================================================
clear; clc; close all;
clear macro_moment_errors;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'..', 'plots'), ...
        fullfile(script_dir, '..','dynamics'), ...
        fullfile(script_dir, '..', 'helpers'), ...
        fullfile(script_dir, '..','static'));
%% STEP 0: Setup
P = set_params();
synth = generate_synthetic_params();
P.w        = synth.w;
P.sig2_eps = synth.sig2_eps;
P.sig2_eta = synth.sig2_eta;
P.q        = P.q;
P.a0       = synth.a0;
P.w_tilde = P.w;
P.y_tilde = P.w_tilde - P.q;
P.kappa = P.r/(1+P.r);
P.pop_share = synth.pop_share;
kappa =  P.r/(1+P.r);
target_ratio = 4.09;

%% STEP 1: Discretize the z grid

P.Nz   = 15;
[P.z_grid] = discretize_normal(P.Nz,0,P.rho, sqrt(max(P.sig2_eps)));

%% STEP 2: Calibrate some parameters (moving cost, amenities, initial wealth and Weibull parameter)
cache_file = fullfile(script_dir, '..', '..', 'Output', 'big_nu', 'calibrated_P.mat');
if isfile(cache_file)
    fprintf('Loading calibrated parameters from cache...\n');
    load(cache_file, 'P');
else
    fprintf('Calibrating parameters\n');
    P = calibrate_param(P);
    save(cache_file, 'P');
    fprintf('Calibration saved to %s\n', cache_file);
end
[current_script_path, ~, ~] = fileparts(mfilename('fullpath'));
project_root = fullfile(current_script_path, '..', '..');
P.border_mat = readmatrix(fullfile(project_root, 'Data', 'border_matrix.csv'));

%% STEP 3: compute static moment at the steady state
fprintf('\n============================================================\n');
fprintf('Step 3: Computing Xi stationary and moments : 2010-2015 baseline\n');
fprintf('============================================================\n');

[Xi_SS] = compute_Xi(P, 1e-7);


[Flow_jl, Z_mean_jl, out_migration_rate, mean_wealth, mean_earn, wealth_to_earn, Lambda, a_bar] = compute_SS_moments(P, Xi_SS);

sd_ratio = sqrt(var(wealth_to_earn));

fprintf('\n============================================================\n');
fprintf('POST-CALIBRATION AUDIT: AGGREGATE MOMENTS\n');
fprintf('============================================================\n');

agg_mob_rate = sum(P.pop_share .* out_migration_rate);

agg_wealth = sum(P.pop_share .* mean_wealth);
agg_earn   = sum(P.pop_share .* mean_earn);
agg_WY_ratio = agg_wealth / agg_earn;

fprintf('Aggregate Mobility Rate      : %.4f  (Target: 0.0250)\n', agg_mob_rate);
fprintf('Wealth-to-Earnings Ratio     : %.4f  (Target: 4.0900)\n', agg_WY_ratio);

pop_model_raw = sum(Lambda, 2);
pop_model = pop_model_raw / sum(pop_model_raw);
pop_error = max(abs(pop_model - P.pop_share));
fprintf('Max Absolute Error in Pop Shares: %e  (Target: < 1e-5)\n', pop_error);
fprintf('\nFirst 5 States Population Share Check:\n');
fprintf('%-15s %-12s %-12s\n', 'State', 'Data', 'Model');
fprintf('-----------------------------------------\n');
for i = 1:5
    fprintf('%-15s %-12.4f %-12.4f\n', P.state_names{i}, P.pop_share(i), pop_model(i));
end

F_data = readmatrix('/Users/leon/Documents/SavingSIAM/Output/data/empirical_F_matrix.csv');

F_movers_emp = F_data;
F_movers_emp(logical(eye(P.J))) = 0;
total_empirical_movers = sum(F_movers_emp(:));
emp_log_dist_weighted = F_movers_emp .* log(1 + P.dist_mat);
target_avg_log_dist = sum(emp_log_dist_weighted(:)) / total_empirical_movers;

Flow_movers_model = Flow_jl;
Flow_movers_model(logical(eye(P.J))) = 0;
total_model_movers = sum(Flow_movers_model(:));
model_log_dist_weighted = Flow_movers_model .* log(1 + P.dist_mat);
avg_log_dist_model = sum(model_log_dist_weighted(:)) / total_model_movers;

dist_error_pct = abs(avg_log_dist_model - target_avg_log_dist) / target_avg_log_dist * 100;

fprintf('\n============================================================\n');
fprintf('MIGRATION DISTANCE VERIFICATION\n');
fprintf('============================================================\n');
fprintf('Empirical Avg Log Distance : %8.4f\n', target_avg_log_dist);
fprintf('Model Avg Log Distance     : %8.4f\n', avg_log_dist_model);
fprintf('Percentage Error           : %8.4f%%\n', dist_error_pct);
fprintf('============================================================\n');

% Compute sorting slope for audit
% Compute sorting slope for audit
J = P.J;
E_z = zeros(J, 1);
for j = 1:J
    mass_j = sum(Lambda(j, :));
    if mass_j > 1e-10
        E_z(j) = sum(Lambda(j, :) .* P.z_grid(:)') / mass_j;
    end
end
log_w = log(P.w(:));
X = [ones(J, 1), log_w];
beta = (X'*X) \ (X'*E_z);
slope_model = beta(2);
e = E_z - X * beta;
s2 = (e' * e) / (J - 2);
V_beta = s2 * inv(X'*X);
se_slope = sqrt(V_beta(2,2));
beta0_H0 = 0.11;
t_stat = (slope_model - beta0_H0) / se_slope;
p_val = 2 * (1 - tcdf(abs(t_stat), J - 2));
fprintf('\n============================================================\n');
fprintf('SPATIAL SORTING VERIFICATION (nu parameter)\n');
fprintf('============================================================\n');
fprintf('Empirical Target Slope     :   0.1100\n');
fprintf('Model Slope                : %8.4f\n', slope_model);
fprintf('Std. Error                 : %8.4f\n', se_slope);
fprintf('Percentage Error           : %8.4f%%\n', abs(slope_model - 0.11)/0.11 * 100);
fprintf('H0: slope = 0.11  t-stat  : %8.4f\n', t_stat);
fprintf('             p-value       : %8.4f\n', p_val);
if p_val > 0.05
    fprintf('=> Cannot reject H0 at 5%% level\n');
else
    fprintf('=> Reject H0 at 5%% level\n');
end
fprintf('============================================================\n');
%% ------------------------------------------------------------------------
% MODEL VS DATA: MIGRATION PATTERN ANALYSIS
% ------------------------------------------------------------------------

Flow_model_movers = Flow_jl;
Flow_model_movers(logical(eye(P.J))) = 0;

Flow_emp_movers = F_data;
Flow_emp_movers(logical(eye(P.J))) = 0;

total_model_movers = sum(Flow_model_movers(:));
total_emp_movers   = sum(Flow_emp_movers(:));

contig_model_moves = sum(sum(Flow_model_movers .* P.border_mat));
contig_emp_moves   = sum(sum(Flow_emp_movers .* P.border_mat));

share_contig_model = contig_model_moves / total_model_movers;
share_contig_emp   = contig_emp_moves / total_emp_movers;

vec_model = Flow_model_movers(:);
vec_emp   = Flow_emp_movers(:);
corr_matrix = corrcoef(vec_model, vec_emp);
global_corr = corr_matrix(1,2);

% 4. Console Report
fprintf('\n========================================================\n');
fprintf('   SPATIAL FRICTIONS: MODEL VS EMPIRICAL COMPARISON\n');
fprintf('========================================================\n');
fprintf('Share of Contiguous Moves (Data)  : %8.4f\n', share_contig_emp);
fprintf('Share of Contiguous Moves (Model) : %8.4f\n', share_contig_model);
fprintf('Absolute Difference               : %8.4f\n', abs(share_contig_model - share_contig_emp));
fprintf('--------------------------------------------------------\n');
fprintf('Global Migration Flow Correlation : %8.4f\n', global_corr);
fprintf('========================================================\n\n');

%% STEP 4: Model vs Empirical Wealth Comparison (Geowealth-US 2010)
graph_dir = fullfile(script_dir, '..', '..', 'Output', 'graph_comparison');
if ~isfolder(graph_dir), mkdir(graph_dir); end

plot_model_vs_empirical_W(P, mean_wealth, graph_dir);
plot_model_vs_empirical_WY(P, wealth_to_earn, graph_dir);
plot_out_sample_fit(P, Lambda, mean_wealth, Flow_jl, F_data, graph_dir);
%test_amenity_calibration(P, graph_dir);
