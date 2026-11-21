function P = generate_synthetic_params()
% =========================================================================
%
% Generates synthetic state-level parameters for Dvorkin & Greaney (2024)
% replication WITHOUT using the original microdata.
%
% Mimics the empirical patterns described in the paper:
%   - Positive corr(w_j, sig2_eps_j) ~ 0.8  (Figure 5b)
%   - Weak   corr(w_j, sig2_eta_j)  ~ 0.2  (Figure 6b)
%   - Wages span roughly [0.887, 1.209] normalized to mean 1
%   - sig2_eps ranges [0.011, 0.023] (Table 9)
%   - sig2_eta ranges [0.052, 0.112] (Table 9)
%   - Mean wealth/earnings ~ 4.09 nationally
%
% OUTPUT:
%   P  — structure with all model parameters
%
% =========================================================================

rng(42);   % fix seed for reproducibility

%% -------------------------------------------------------------------------
%  1. STATE NAMES & LABELS
% -------------------------------------------------------------------------
P.state_names = { ...
    'Alabama','Alaska', 'Arizona','Arkansas','California', ...
    'Colorado','Connecticut','Delaware','District of Columbia','Florida', ...
    'Georgia','Hawaii','Idaho','Illinois','Indiana', ...
    'Iowa','Kansas','Kentucky','Louisiana','Maine', ...
    'Maryland','Massachusetts','Michigan','Minnesota','Mississippi', ...
    'Missouri','Montana','Nebraska','Nevada','New Hampshire', ...
    'New Jersey','New Mexico','New York','North Carolina','North Dakota', ...
    'Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island', ...
    'South Carolina','South Dakota','Tennessee','Texas','Utah', ...
    'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming'};

P.J = numel(P.state_names);   % 51

%% -------------------------------------------------------------------------
%  2. GLOBAL CALIBRATED PARAMETERS
% -------------------------------------------------------------------------
P.r          = 0.03;            % real interest rate
P.beta       = 1 / (1 + P.r);  % discount factor
P.gamma      = 3;               % CARA risk aversion
P.rho        = 0.934;           % AR(1) persistence (from PSID ratio)
P.nu         = 3.95;            % Weibull shape  (matched to sorting)
P.delta      = 0.025;           % death / replacement probability
P.tau        = 0;            % proportional income tax
P.lambda_bar = 0.372;           % Weibull scale (matched to W/Y = 4.09)
P.a0         = 0;               % initial assets of newborns

% Derived scalars
P.kappa      = P.r / (1 + P.r);         % stationary MPS (= r/(1+r))
P.r_tilde    = 1 / (1 + P.r);           % discount in Xi recursion
P.beta_tilde = P.beta * (1 + P.r);      % = 1 exactly when beta=1/(1+r)


RPP_2015 = [
    59.246;   % Alabama
    114.548;  % Alaska
    93.772;   % Arizona
    59.372;   % Arkansas
    158.044;  % California
    123.873;  % Colorado
    132.470;  % Connecticut
    95.175;   % Delaware
    178.359;  % District of Columbia
    109.263;  % Florida
    80.455;   % Georgia
    153.726;  % Hawaii
    74.196;   % Idaho
    103.287;  % Illinois
    73.773;   % Indiana
    70.031;   % Iowa
    75.063;   % Kansas
    66.279;   % Kentucky
    71.514;   % Louisiana
    77.047;   % Maine
    131.505;  % Maryland
    130.005;  % Massachusetts
    80.739;   % Michigan
    97.248;   % Minnesota
    58.588;   % Mississippi
    73.226;   % Missouri
    73.850;   % Montana
    75.410;   % Nebraska
    98.742;   % Nevada
    119.650;  % New Hampshire
    145.268;  % New Jersey
    78.836;   % New Mexico
    134.327;  % New York
    75.763;   % North Carolina
    78.073;   % North Dakota
    75.578;   % Ohio
    66.381;   % Oklahoma
    105.383;  % Oregon
    88.932;   % Pennsylvania
    113.131;  % Rhode Island
    74.555;   % South Carolina
    62.178;   % South Dakota
    72.763;   % Tennessee
    93.957;   % Texas
    93.617;   % Utah
    98.901;   % Vermont
    114.868;  % Virginia
    119.584;  % Washington
    61.693;   % West Virginia
    81.187;   % Wisconsin
    83.272    % Wyoming
];

alpha_housing = 0.30;
P.q = alpha_housing * (RPP_2015 / 100); 
   % housing costs q_j  (from BEA)
        % disposable income net of rent [J x 1]

%% -------------------------------------------------------------------------
%  8. MEAN PERSISTENT SHOCK BY STATE  E[z|j]
%     From Table 9: ranges [-0.028, 0.073], positively correlated with w
%     We use actual Table 9 values
% -------------------------------------------------------------------------

Ez_table9 = [ ...
   -0.007; 0.015;  -0.006; -0.010;  0.011; ...  % AL AK AZ AR CA
    0.006;  0.011;  0.005;  0.073; -0.011; ...  % CO CT DE DC FL
    0.002;  0.005; -0.011;  0.003; -0.008; ...  % GA HI ID IL IN
   -0.004; -0.011; -0.011; -0.002; -0.008; ...  % IA KS KY LA ME
    0.013;  0.017; -0.008;  0.001; -0.018; ...  % MD MA MI MN MS
   -0.009; -0.019; -0.006;  0.002;  0.010; ...  % MO MT NE NV NH
    0.016; -0.015;  0.012; -0.011; -0.015; ...  % NJ NM NY NC ND
   -0.007; -0.008; -0.008; -0.005; -0.006; ...  % OH OK OR PA RI
   -0.004; -0.010; -0.004;  0.004; -0.007; ...  % SC SD TN TX UT
   -0.005;  0.006;  0.012; -0.013; -0.007; ...  % VT VA WA WV WI
   -0.028];                                       % WY

P.Ez = Ez_table9;

%% -------------------------------------------------------------------------
%  10. POPULATION SHARES  pop_j
%     From actual US state populations (normalized)
% -------------------------------------------------------------------------

% Approximate state population shares (2015 vintage, normalized)
pop_raw = [ ...
    4.858; 0.738 ; 6.828; 2.978; 39.145; ...  % AL AK AZ AR CA
    5.456; 3.591; 0.945; 0.672; 20.272; ...  % CO CT DE DC FL
   10.214; 1.431; 1.654; 12.860; 6.619; ...  % GA HI ID IL IN
    3.123; 2.911; 4.425; 4.671; 1.329; ...   % IA KS KY LA ME
    6.006; 6.795; 9.923; 5.490; 2.992; ...   % MD MA MI MN MS
    6.084; 1.033; 1.896; 2.890; 1.331; ...   % MO MT NE NV NH
    8.958; 2.085; 19.795; 10.042; 0.757; ... % NJ NM NY NC ND
   11.614; 3.911; 4.029; 12.803; 1.056; ...  % OH OK OR PA RI
    4.896; 0.858; 6.600; 27.470; 3.001; ...  % SC SD TN TX UT
    0.626; 8.383; 7.171; 1.844; 5.772; ...   % VT VA WA WV WI
    0.586];                                   % WY

P.pop_share = pop_raw / sum(pop_raw);   % [J x 1], sums to 1


%% -------------------------------------------------------------------------
%  3. GENERATE STATE-LEVEL WAGES  w_j
%     Paper: mean ~1, range [0.887, 1.209], log-std ~ 0.072 (Table 9)
% -------------------------------------------------------------------------

% Use actual values from Table 9 in the paper
w_table9 = [ ...
    0.914; 1.114; 0.977; 0.890; 1.094; ...  % AL AK AZ AR CA
    1.027; 1.125; 1.005; 1.209; 0.942; ...  % CO CT DE DC FL
    0.952; 1.047; 0.922; 1.021; 0.941; ...  % GA HI ID IL IN
    0.929; 0.931; 0.910; 0.935; 0.934; ...  % IA KS KY LA ME
    1.125; 1.122; 0.973; 1.051; 0.894; ...  % MD MA MI MN MS
    0.956; 0.902; 0.934; 1.006; 1.033; ...  % MO MT NE NV NH
    1.124; 0.949; 1.055; 0.947; 0.951; ...  % NJ NM NY NC ND
    0.957; 0.907; 1.022; 0.990; 1.037; ...  % OH OK OR PA RI
    0.915; 0.894; 0.896; 0.952; 0.974; ...  % SC SD TN TX UT
    0.947; 1.056; 1.106; 0.887; 0.970; ...  % VT VA WA WV WI
    0.964];                                  % WY

w_avg = sum(w_table9 .* P.pop_share);
P.w = w_table9 / w_avg;


%% -------------------------------------------------------------------------
%  4. GENERATE STATE-LEVEL PERSISTENT SHOCK VARIANCE  sigma2_eps_j
%     Paper: range [0.011, 0.023], corr with log(w) ~ 0.8 (Figure 5b)
%     These are the INNOVATION variances sigma2_epsilon
%     Stationary variance = sig2_eps / (1 - rho^2)
% -------------------------------------------------------------------------

sig2_eps_table9 = [ ...
    0.011; 0.014; 0.015; 0.013; 0.021; ...  % AL AK AZ AR CA
    0.018; 0.017; 0.013; 0.023; 0.015; ...  % CO CT DE DC FL
    0.015; 0.016; 0.017; 0.020; 0.015; ...  % GA HI ID IL IN
    0.015; 0.018; 0.015; 0.016; 0.017; ...  % IA KS KY LA ME
    0.020; 0.017; 0.014; 0.017; 0.014; ...  % MD MA MI MN MS
    0.016; 0.015; 0.015; 0.014; 0.019; ...  % MO MT NE NV NH
    0.023; 0.017; 0.019; 0.015; 0.013; ...  % NJ NM NY NC ND
    0.016; 0.014; 0.017; 0.018; 0.015; ...  % OH OK OR PA RI
    0.014; 0.013; 0.017; 0.018; 0.017; ...  % SC SD TN TX UT
    0.014; 0.021; 0.021; 0.013; 0.016; ...  % VT VA WA WV WI
    0.011];                                  % WY

P.sig2_eps = sig2_eps_table9;

%% -------------------------------------------------------------------------
%  5. GENERATE STATE-LEVEL TRANSITORY SHOCK VARIANCE  sigma2_eta_j
%     Paper: range [0.052, 0.112], corr with log(w) ~ 0.2 (Figure 6b)
% -------------------------------------------------------------------------

sig2_eta_table9 = [ ...
    0.098; 0.108; 0.092; 0.094; 0.100; ...  % AL AK AZ AR CA
    0.094; 0.093; 0.091; 0.089; 0.111; ...  % CO CT DE DC FL
    0.109; 0.109; 0.073; 0.090; 0.075; ...  % GA HI ID IL IN
    0.077; 0.076; 0.083; 0.100; 0.062; ...  % IA KS KY LA ME
    0.099; 0.107; 0.086; 0.066; 0.094; ...  % MD MA MI MN MS
    0.078; 0.080; 0.066; 0.112; 0.070; ...  % MO MT NE NV NH
    0.108; 0.102; 0.112; 0.106; 0.073; ...  % NJ NM NY NC ND
    0.064; 0.103; 0.075; 0.073; 0.077; ...  % OH OK OR PA RI
    0.087; 0.067; 0.083; 0.093; 0.079; ...  % SC SD TN TX UT
    0.075; 0.091; 0.070; 0.099; 0.052; ...  % VT VA WA WV WI
    0.075];                                  % WY

P.sig2_eta = sig2_eta_table9;

%% -------------------------------------------------------------------------
%  7. DERIVED WAGE OBJECTS
% -------------------------------------------------------------------------
P.w_tilde = P.w * (1 - P.tau);          % after-tax wage  [J x 1]
P.y_tilde = P.w_tilde - P.q;    

%% -------------------------------------------------------------------------
%  11. RESIDUAL LOG-WAGES  (synthetic individual-level moments by state)
%     Generate state-level moments of residual log wages that are
%     consistent with the model income process:
%       zeta_{i,j} = z_i + eta_i
%       Var(zeta|j) = Var(z|j) + sig2_eta_j
%       Cov(zeta_t, zeta_{t-1}|j, stayer) = rho * Var(z|j)
% -------------------------------------------------------------------------

% Stationary variance of persistent shock z in each state
P.sig2_z_stationary = P.sig2_eps / (1 - P.rho^2);   % [J x 1]

% State-level moments of residual log-earnings
%   Var(zeta|j) = sig2_z_stationary + sig2_eta  [total residual variance]
P.var_resid = P.sig2_z_stationary + P.sig2_eta;      % [J x 1]

%   Cov(zeta_t, zeta_{t-1}|j, stayer) = rho^2 * sig2_z_stat
%   (one-period-lag autocovariance for stayers)
P.cov_resid_lag1 = P.rho^2 * P.sig2_z_stationary;   % [J x 1]

%   Ratio used to identify rho: cov(t,t-2)/cov(t,t-1) = rho
P.cov_resid_lag2 = P.rho^3 * P.sig2_z_stationary;   % [J x 1]

% Mean residual log-wage by state (= E[z|j], from selection)
P.mean_resid = P.Ez;    % [J x 1]


fprintf('Synthetic parameter generation complete.\n');
end % function generate_synthetic_params
