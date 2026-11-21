function W_atkinson = compute_atkinson_welfare(c_flow_mat, Lambda_mat, epsilon)
% COMPUTE_ATKINSON_WELFARE Calculates inequality-adjusted social welfare directly
% over positive flow consumption levels c_flow = y_tilde + w_tilde*z + kappa*a > 0.
%
% Inputs:
%   c_flow_mat - J x Nz matrix of positive flow consumption levels ($/yr)
%   Lambda_mat - J x Nz matrix of population distribution
%   epsilon    - Inequality aversion parameter (epsilon >= 0)
%
% Output:
%   W_atkinson - Inequality-adjusted social welfare ($/yr)

pop_total = sum(Lambda_mat(:));
weights   = Lambda_mat(:) / pop_total;
c_vec     = c_flow_mat(:);

% Validate positive flow consumption domain
if any(c_vec <= 0)
    error('compute_atkinson_welfare:NonPositiveConsumption', ...
          'Flow consumption must be strictly positive (c_flow > 0).');
end

if abs(epsilon - 1.0) < 1e-6
    % Utilitarian Log (Bernoulli) form
    W_atkinson = exp(sum(weights .* log(c_vec)));
elseif epsilon == 0
    % Standard Linear Dollar Utilitarian form
    W_atkinson = sum(weights .* c_vec);
else
    % General Atkinson Power form
    power_sum  = sum(weights .* (c_vec .^ (1 - epsilon)));
    W_atkinson = (power_sum ^ (1 / (1 - epsilon)));
end
end
