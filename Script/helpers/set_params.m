%set the parameters of the model
function P = set_params()

P.state_names = { ...
        'Alabama','Alaska','Arizona','Arkansas','California', ...
        'Colorado','Connecticut','Delaware','District of Columbia','Florida', ...
        'Georgia','Hawaii','Idaho','Illinois','Indiana', ...
        'Iowa','Kansas','Kentucky','Louisiana','Maine', ...
        'Maryland','Massachusetts','Michigan','Minnesota','Mississippi', ...
        'Missouri','Montana','Nebraska','Nevada','New Hampshire', ...
        'New Jersey','New Mexico','New York','North Carolina','North Dakota', ...
        'Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island', ...
        'South Carolina','South Dakota','Tennessee','Texas','Utah', ...
        'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming'};

P.J       = 51;          % number of states
P.r       = 0.03;        % interest rate
P.beta    = 1/(1+P.r);   % discount factor
P.gamma   = 3;           % CARA risk aversion / prudence
P.rho     = 0.934;       % AR(1) persistence of z
P.nu      = 3.95;        % Weibull shape parameter
P.delta   = 0.025;       % death probability (2.5% annual mobility ≈ replacement)
P.lambda_bar = 0.2;    % Weibull scale (matched to national W/Y = 4.09)
P.a0      = 0;           % initial assets of newborns

P.w       = [];   % unit wages w_j
P.sig2_eps = [];  % variance of persistent shock innovation
P.sig2_eta = [];  % variance of transitory shock

P.A = zeros(P.J,1);   % amenities A_ell 

P.dist_mat = state_distances(P);
P.psi_d = 5;
P.psi = zeros(P.J, P.J);
    for j = 1:P.J
        for l = 1:P.J
            if j ~= l
                    P.psi(j, l) = P.psi_d * log(1 + P.dist_mat(j, l));
            end
        end
    end

% Derived
P.w_tilde = [];   % w_tilde = w
P.y_tilde = [];   % y_tilde = w - q

%state capacity and land 
[current_script_path, ~, ~] = fileparts(mfilename('fullpath'));
project_root = fullfile(current_script_path, '..', '..');
housing_file_path = fullfile(project_root, 'Data', 'ACSDT1Y2015.csv');
opts = detectImportOptions(housing_file_path);
opts.SelectedVariableNames = {'NAME', 'B25001_001E'};
opts.DataLines = [2, Inf];
housing_raw = readtable(housing_file_path, opts);
H = zeros(P.J, 1);
for i = 1:P.J
    idx = strcmp(housing_raw.NAME, P.state_names{i});
    H(i) = housing_raw.B25001_001E(idx);
end
P.state_capacity = H;
P.border_mat = readmatrix(fullfile(project_root, 'Data', 'border_matrix.csv'));

P.coastal_dummy = [
    0; % Alabama
    1; % Alaska
    0; % Arizona
    0; % Arkansas
    1; % California (West Coast)
    0; % Colorado
    1; % Connecticut
    1; % Delaware
    1; % District of Columbia
    1; % Florida
    1; % Georgia
    0; % Hawaii
    0; % Idaho
    0; % Illinois
    0; % Indiana
    0; % Iowa
    0; % Kansas
    0; % Kentucky
    1; % Louisiana (Gulf Coast)
    1; % Maine
    1; % Maryland
    1; % Massachusetts
    0; % Michigan
    0; % Minnesota
    1; % Mississippi (Gulf Coast)
    0; % Missouri
    0; % Montana
    0; % Nebraska
    0; % Nevada
    1; % New Hampshire
    1; % New Jersey
    0; % New Mexico
    1; % New York
    1; % North Carolina
    0; % North Dakota
    0; % Ohio
    0; % Oklahoma
    0; % Oregon
    1; % Pennsylvania
    1; % Rhode Island
    1; % South Carolina
    0; % South Dakota
    0; % Tennessee
    1; % Texas
    0; % Utah
    0; % Vermont (New England, but landlocked)
    1; % Virginia
    1; % Washington
    0; % West Virginia
    0; % Wisconsin
    0  % Wyoming
];


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
end


