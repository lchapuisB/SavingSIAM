function [lambda_state, lambda_nat] = compute_consumption_equivalents(res_counterfactual, res_baseline, P_base)
% COMPUTE_CONSUMPTION_EQUIVALENTS Calculates percentage flow Certainty Equivalent
% gain/loss relative to empirical baseline flow consumption:
%
%   Delta % CE_{l,t} = (CE^1_{l,t} - CE^0_{l,t}) / (y_tilde^0_l + kappa * a_bar^0_{l,t})
%
% Inputs:
%   res_counterfactual - Struct from run_scenario for counterfactual (contains mean_ce_t, mean_wealth_t, Lambda_trans)
%   res_baseline       - Struct from run_scenario for baseline (contains mean_ce_t, mean_wealth_t, Lambda_trans)
%   P_base             - Parameter struct (contains J, r, y_tilde)
%
% Outputs:
%   lambda_state - J x T matrix of state-level % flow CE changes
%   lambda_nat   - 1 x T vector of national aggregate % flow CE changes

J = P_base.J;
T = size(res_baseline.mean_ce_t, 2);
kappa = P_base.r / (1 + P_base.r);

lambda_state = zeros(J, T);
lambda_nat   = zeros(1, T);

y_tilde_0 = P_base.y_tilde; % J x 1 baseline disposable wage vector

for t = 1:T
    CE_0 = res_baseline.mean_ce_t(:, t);        % J x 1 state mean CEs
    CE_1 = res_counterfactual.mean_ce_t(:, t); % J x 1 state mean CEs
    
    a_0  = res_baseline.mean_wealth_t(:, t);    % J x 1 state mean wealth
    a_1  = res_counterfactual.mean_wealth_t(:, t); % J x 1 state mean wealth
    
    Lam_0 = res_baseline.Lambda_trans(:, :, t); % J x Nz population tensor
    Lam_1 = res_counterfactual.Lambda_trans(:, :, t); % J x Nz population tensor
    
    pop_0 = sum(Lam_0, 2); % J x 1 state populations
    pop_1 = sum(Lam_1, 2); % J x 1 state populations
    
    % Baseline flow consumption: disposable wage + asset annuity
    base_flow_cons_0 = y_tilde_0 + kappa * a_0; % J x 1
    
    % State-level % flow CE change relative to baseline flow consumption
    lambda_state(:, t) = (CE_1 - CE_0) ./ abs(base_flow_cons_0);
    
    % National aggregate numerator and denominator
    num_nat = sum((CE_1 - CE_0) .* pop_1);
    den_nat = sum(base_flow_cons_0 .* pop_0);
    
    if den_nat > 0
        lambda_nat(t) = num_nat / den_nat;
    end
end
end
