function save_data(P, suffix)
% save_data.m
if nargin < 2
    suffix = '';
end
script_dir = fileparts(mfilename('fullpath'));
% I will just re-run the loop for a few points and print them
addpath(fullfile(script_dir, '..', 'static'), fullfile(script_dir, '..', 'helpers'), fullfile(script_dir, '..', 'dynamics'));
FL_idx = 10;

wage_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_wage_2080.csv');
opts = detectImportOptions(wage_damages_path);
damages_wage = readtable(wage_damages_path, opts);
amen_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_amenity_2080.csv');
opts = detectImportOptions(amen_damages_path);
damages_amen = readtable(amen_damages_path, opts);
housing_damages_path = fullfile(script_dir, '..', '..', 'Output', 'data', 'shock_housing_2080.csv');
opts = detectImportOptions(housing_damages_path);
damages_housing = readtable(housing_damages_path, opts);

shock_housing = damages_housing.shock_housing;
shock_amen    = damages_amen.shock_amenity;
shock_wage    = damages_wage.shock_wage;

mults = [0.01, 0.5, 0.9, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0];
fprintf('mult, drop\n');
for i = 1:length(mults)
    mult = mults(i);
    P_ctrl = P; P_ctrl.psi = P.psi * mult;
    Xi_ctrl = compute_Xi(P_ctrl);
    [~, ~, ~, ~, ~, ~, Lam_ctrl, ~] = compute_SS_moments(P_ctrl, Xi_ctrl);
    
    P_treat = P_ctrl;
    P_treat.w = P_treat.w .* (1 + shock_wage);
    P_treat.q = P_treat.q .* (1 + shock_housing);
    P_treat.A = P_treat.A - shock_amen .* P_treat.w_tilde;
    P_treat.w_tilde = P_treat.w_tilde .* (1 + shock_wage);
    P_treat.y_tilde = P_treat.w_tilde - P_treat.q;

    if isfield(P_treat, 'uncertainty')
        if P_treat.uncertainty == 1
            P_treat.sig2_eta = P_ctrl.sig2_eta;
        elseif P_treat.uncertainty == 2
            P_treat.sig2_eta = P_ctrl.sig2_eta .* (1 + 50 * abs(shock_wage));
        end
    end
    Xi_treat = compute_Xi(P_treat);
    [~, ~, ~, ~, ~, ~, Lam_treat, ~] = compute_SS_moments(P_treat, Xi_treat);
    
    drop = (sum(Lam_treat(FL_idx,:)) - sum(Lam_ctrl(FL_idx,:))) * 100;
    fprintf('%.2f, %.6f\n', mult, drop);
end
end
