# Mathematical Documentation: Dynamic Spatial Wealth and Climate Shock Framework

**Project**: SavingSIAM — Dynamic Quantitative Spatial Economics & Precautionary Savings Under Climate Change  
**Author / Documentation Specialist**: Technical Writer & Mathematical Documentation Specialist  
**Date**: 2026-07-25  
**Target Document**: `/Users/leon/Documents/SavingSIAM/Output/mathematical_documentation.md`  

---

## 1. Executive Summary & Model Overview

### 1.1 Theoretical Scope and Literature Position
This document delivers the comprehensive, publication-grade mathematical documentation for the **SavingSIAM** dynamic quantitative spatial economics framework. The model evaluates the interaction between climate change shocks, spatial reallocation, dynamic precautionary savings, and aggregate welfare across $J = 51$ U.S. states over the 2010–2100 transition period and into long-run steady states ($t = 500$).

The theoretical core builds upon and extends the benchmark dynamic quantitative spatial model of **Dvorkin & Greaney (2024)** ("*The Geography of Wealth: Shocks, Mobility, and Precautionary Savings*"). The key analytical innovation of this class of models is the adoption of **Constant Absolute Risk Aversion (CARA) preferences** combined with **Extreme Value Type I (Weibull) location choice shocks**. CARA utility eliminates wealth-income scale dependencies in spatial moving decisions, enabling a closed-form solution for optimal consumption-savings policy functions that are linear in individual wealth holdings $a$. 

This structural decoupling yields a tractable **three-stage sequential solution algorithm**:
1. **Backward Value Function Index Induction**: Solves for the location-productivity value index matrix $\boldsymbol{\Xi}_t \in \mathbb{R}^{J \times N_z}$ backward from terminal period $t=T$ to $t=1$, completely independently of individual asset distributions.
2. **Forward Population Propagation**: Propagates the joint population mass array $\boldsymbol{\Lambda}_t \in \mathbb{R}^{J \times N_z}$ forward from $t=1$ to $t=T$ using sparse time-varying transition matrices $\mathbf{T}_t$, independent of wealth holdings.
3. **Forward Linear Wealth Propagation**: Propagates per-capita asset holdings $\bar{a}_t(\ell, z)$ forward using affine savings rules and survivor wealth transfer operators.

### 1.2 Scenario Taxonomy Matrix
The codebase implements 9 primary structural scenario runs (Scenarios A through I), systematically exploring expectation regimes, structural saving constraints, climate uncertainty, and spatial mobility barriers.

| Scenario Identifier | Expectation Regime | Structural Variant | Migration Friction ($\boldsymbol{\psi}$) | Idiosyncratic Variance ($\sigma^2_\eta$) | Output Subdirectory |
|:---|:---|:---|:---|:---|:---|
| **Scenario A** | Perfect Foresight | Standard (Saving Allowed) | Baseline Calibrated $\boldsymbol{\psi}$ | Baseline ($\text{uncertainty} = 0$) | `Output/graph_perfect_foresight` |
| **Scenario B** | Static (Myopia) | Standard (Saving Allowed) | Baseline Calibrated $\boldsymbol{\psi}$ | Baseline ($\text{uncertainty} = 0$) | `Output/graph_myopia` |
| **Scenario C** | Perfect Foresight | Variance Shocks (Uncertainty) | Baseline Calibrated $\boldsymbol{\psi}$ | $\Delta \sigma^2_{\eta,\ell} = 50 \cdot \|\Delta w_\ell\|$ | `Output/graph_uncertainty_variance_shocks` |
| **Scenario D** | Perfect Foresight | Hand-to-Mouth (No Saving) | Baseline Calibrated $\boldsymbol{\psi}$ | Baseline ($\text{uncertainty} = 0$) | `Output/graph_no_saving_struct` |
| **Scenario E** | Perfect Foresight | Autarkic DiD (No Migration) | $\psi_{j,\ell} = 10^8$ for $j \neq \ell$ | Baseline ($\text{uncertainty} = 0$) | `Output/graph_no_migration_did` |
| **Scenario F** | Perfect Foresight | Small Migration DiD (+10% Cost) | $1.1 \times \boldsymbol{\psi}$ | Baseline ($\text{uncertainty} = 0$) | `Output/graph_small_migration_did` |
| **Scenario G** | Perfect Foresight | Big Migration DiD (-10% Cost) | $0.9 \times \boldsymbol{\psi}$ | Baseline ($\text{uncertainty} = 0$) | `Output/graph_big_migration_did` |
| **Scenario H** | Perfect Foresight | Variance Shocks + Autarkic DiD | $\psi_{j,\ell} = 10^8$ for $j \neq \ell$ | $\Delta \sigma^2_{\eta,\ell} = 50 \cdot \|\Delta w_\ell\|$ | `Output/graph_uncertainty_variance_shocks_no_migration_did` |
| **Scenario I** | Perfect Foresight | Damage Path Volatility | Baseline Calibrated $\boldsymbol{\psi}$ | Path Volatility ($\text{uncertainty} = 1$) | `Output/graph_uncertainty_damage_path` |

---

## 2. Theoretical Foundations & Micro-Foundations

### 2.1 Demographic Structure (Blanchard-Yaari Perpetual Youth)
The economy consists of a continuum of infinitely lived households facing an exogenous annual mortality probability $\delta > 0$. In each period $t$, a fraction $\delta$ of households dies and is replaced by an equal measure of newborn households, maintaining a constant total population normalized to $L_{\text{tot}} = 1$.

- **Mortality / Replacement Rate**: $\delta = 0.025$ (2.5% annual replacement rate, matching U.S. annual baseline gross migration/turnover in `set_params.m:22`).
- **Newborn Initial Wealth**: $a_0 = 0$ (`set_params.m:25`).
- **Newborn Initial Productivity Shock**: Newborns enter location $\ell$ at the median persistent productivity state $z_0 = \lceil N_z / 2 \rceil = 4$ for $N_z = 7$.

### 2.2 Income Processes and Idiosyncratic Risk
A worker residing in state $\ell \in \{1, \ldots, J\}$ with persistent productivity index $z_i \in \{z_1, \ldots, z_{N_z}\}$ and transitory shock $\eta_{i,t}$ earns gross labor income:
$$y_{i,\ell,t} = w_{\ell,t} \left(1 + z_{i,t} + \eta_{i,t}\right)$$

where $w_{\ell,t}$ is the base wage per unit of efficiency in state $\ell$. 

- **Taxation and Net Wage**: Disposable wage income is subject to proportional tax $\tau = 0.25$, yielding net wage $\tilde{w}_{\ell,t} = w_{\ell,t}(1 - \tau)$.
- **Housing Rent**: Each worker incurs local housing cost $q_{\ell,t} = \alpha_{\text{housing}} \cdot (\text{RPP}_{\ell} / 100)$, with $\alpha_{\text{housing}} = 0.30$ (`set_params.m:172-173`).
- **Net Base Disposable Income**:
  $$\tilde{y}_{\ell,t} = \tilde{w}_{\ell,t} - q_{\ell,t}$$

- **Persistent Shock Process**: Log efficiency $z$ follows a state-specific AR(1) process:
  $$z_{i,t+1} = \rho z_{i,t} + \varepsilon_{i,t+1}, \quad \varepsilon_{i,t+1} \sim \mathcal{N}\left(0, \sigma^2_{\varepsilon,\ell,t}\right)$$
  with persistence $\rho = 0.934$ (`set_params.m:20`). The process is discretized into $N_z = 7$ states using the Rouwenhorst method implemented in `discretize_normal.m`.
- **Transitory Shock Process**: $\eta_{i,t} \sim \mathcal{N}\left(0, \sigma^2_{\eta,\ell,t}\right)$, independent of $z$.

### 2.3 Closed-Form Policy Functions (Proposition 1)
Agents have Constant Absolute Risk Aversion (CARA) preferences over non-housing consumption $c$ and local amenity $A_\ell$:
$$u(c, A_\ell) = -\frac{1}{\gamma} \exp\left(-\gamma (c - A_\ell)\right)$$
where $\gamma = 3$ is the coefficient of absolute risk aversion / prudence (`set_params.m:19`).

The budget constraint for an agent in state $\ell$ holding wealth $a_t$ is:
$$a_{t+1} = (1 + r_t)\left( a_t + \tilde{y}_{\ell,t} + \tilde{w}_{\ell,t} z_i + \tilde{w}_{\ell,t} \eta_{i,t} - c_{t} \right)$$
where $r = 0.03$ is the real interest rate (`set_params.m:17`) and discount factor $\beta = \frac{1}{1+r} \approx 0.97087$ (`set_params.m:18`).

#### Proposition 1 (CARA Savings Policy Decoupling):
The value function factorizes additively into wealth holdings $a$ and non-wealth state index $\Xi_{\ell,t}(z)$:
$$V_{\ell,t}(a, z) = -\frac{1}{\gamma} \exp\left( -\gamma \kappa_{t+1} a - \gamma \Xi_{\ell,t}(z) \right)$$

where $\kappa_{t+1}$ is the marginal propensity to consume (MPC) out of wealth:
$$\kappa_{t+1} = \frac{r_t}{1 + r_t}$$
For stationary interest rate $r = 0.03$, $\kappa^* = \frac{0.03}{1.03} \approx 0.029126$.

The optimal consumption and next-period asset policy functions are:
$$c_{\ell,t}(a, z', \eta') = \Xi_{\ell,t}(z') + \kappa_{t+1} \left[ \tilde{y}_{\ell,t} + \tilde{w}_{\ell,t} \eta' + a \right]$$
$$a_{\ell,t+1}(a, z', \eta') = (1+r_t)(1-\kappa_{t+1}) \left[ \tilde{y}_{\ell,t} + \tilde{w}_{\ell,t} \eta' + a \right] + (1+r_t)\left[ \tilde{w}_{\ell,t} z' - \Xi_{\ell,t}(z') \right]$$

### 2.4 Value Function Index Recursion ($\Xi$)
The location-productivity index $\Xi_{\ell,t}(z)$ satisfies the non-linear Bellman recursion equation (equation (9) of Dvorkin & Greaney 2024, implemented in `compute_Xi.m:74-83` and `backward_step_Xi.m:50-65`):

$$\Xi_{\ell,t}(z) = -\frac{1}{\gamma(1+r_t)} \left[ \ln\left(\beta (1+r_t)\right) - \gamma r_t \tilde{w}_{\ell,t} z + A_{\ell,t} + \ln \bar{\lambda} - \frac{1}{\nu} \text{ls}_{\ell,t}(z) \right]$$

where $\text{ls}_{\ell,t}(z)$ is the log-sum-exp location value aggregator over all potential destination states $m \in \{1, \ldots, J\}$:
$$\text{ls}_{\ell,t}(z) = \ln \sum_{m=1}^J \exp\left( \nu B_{1,\ell m, t+1} - \nu \ln \mathbb{E}_{z'|z,m}\left[ \exp\left(-\gamma \Xi_{m,t+1}(z')\right) \right] \right)$$

and deterministic utility component $B_{1,\ell m, t+1}$ is:
$$B_{1,\ell m, t+1} = A_{m,t+1} - \psi_{\ell m} + \gamma \kappa_{t+1} \tilde{y}_{m,t+1} - \frac{1}{2}\left(\kappa_{t+1} \gamma \tilde{w}_{m,t+1}\right)^2 \sigma^2_{\eta,m,t+1}$$

### 2.5 Spatial Selection and Weibull Mobility Choice
At the end of period $t$, workers observe destination-specific i.i.d. preference shocks $\epsilon_m$ drawn from an Extreme Value Type I (Weibull) distribution with shape parameter $\nu = 3.65$ (`set_params.m:21`) and scale parameter $\bar{\lambda} = 0.5257$ (`set_params.m:24`).

Moving from origin $\ell$ to destination $m$ incurs bilateral migration cost $\psi_{\ell m}$:
$$\psi_{\ell m} = \psi_{\text{fixed}} + \psi_{\text{dist}} \ln(1 + d_{\ell m}) \quad \text{for } \ell \neq m, \quad \psi_{\ell \ell} = 0$$
where $d_{\ell m}$ is the Haversine great-circle distance between state centroids in miles (`helpers/state_distances.m`).

Conditional on persistent productivity $z$, the probability of a worker residing in origin state $j$ choosing destination state $\ell$ in period $t$ is given by the 3D tensor $\boldsymbol{\mu}_t \in \mathbb{R}^{J \times J \times N_z}$ (`helpers/compute_mu_all.m:35-50`):
$$\mu_{j \to \ell, t}(z) = \frac{\exp\left( \nu B_{1,j\ell,t} - \nu \ln \mathbb{E}_{z'|z,\ell}\left[ e^{-\gamma \Xi_{\ell,t+1}(z')} \right] \right)}{\sum_{m=1}^J \exp\left( \nu B_{1,jm,t} - \nu \ln \mathbb{E}_{z'|z,m}\left[ e^{-\gamma \Xi_{m,t+1}(z')} \right] \right)}$$

---

## 3. Static Model & Two-Layer Calibration

### 3.1 Structural Calibration (`calibrate_param.m`)
The baseline 2010 steady state is calibrated via a two-layer nested optimization framework matching macro targets and spatial population distributions:

```
+-----------------------------------------------------------------------------------+
| OUTER LOOP: MATLAB fsolve on global parameters theta = [psi_fixed, psi_dist, lambda_bar] |
| Matches: 1. National gross mobility rate m = 2.5%                                  |
|          2. Mean log migration distance                                           |
|          3. Aggregate wealth-to-earnings ratio W/Y = 4.0900                       |
+-----------------------------------------------------------------------------------+
                                         |
                                         v
+-----------------------------------------------------------------------------------+
| INNER LOOP: Berry-Levinsohn-Pakes (BLP 1995) contraction mapping on amenities A_ell|
| Formula: A_ell^(k+1) = A_ell^(k) + (1/nu) * [ ln(pop_target_ell) - ln(pop_model_ell) ]|
| Convergence criterion: max_ell |pop_model_ell - pop_target_ell| < 10^(-5)             |
+-----------------------------------------------------------------------------------+
```

### 3.2 Direct Matrix Inversion for Stationary Population Distribution
Let $S = J \times N_z = 51 \times 7 = 357$ be the total state space size. The composite state index $s = (\ell, z)$ tracks spatial location $\ell$ and productivity index $z$.

The full sparse transition matrix $\mathbf{T} \in \mathbb{R}^{S \times S}$ (`helpers/build_T_t.m:30-65`) governs worker transitions from origin state $(j, z)$ to destination state $(\ell, z')$:
$$T_{(j,z), (\ell,z')} = (1 - \delta) \mu_{j \to \ell}(z) \Pi^\ell_{z, z'} + \delta \mu_{j \to \ell}(z_0) \Pi^\ell_{z_0, z'}$$

The stationary population distribution vector $\boldsymbol{\Lambda} \in \mathbb{R}^{1 \times S}$ solves the fixed-point matrix equation:
$$\boldsymbol{\Lambda} = \boldsymbol{\Lambda} \mathbf{T} \iff \boldsymbol{\Lambda} (\mathbf{I} - \mathbf{T}) = \mathbf{0}$$

In computational MATLAB code (`static/compute_SS_moments.m:40-48`), $\boldsymbol{\Lambda}$ is extracted directly as the dominant left eigenvector associated with eigenvalue $\lambda_1 = 1$:
```matlab
[eig_vec, ~] = eigs(T_mat_full', 1, 'lr');
Lambda_vec = (eig_vec / sum(eig_vec))';
```

### 3.3 Neumann Linear Solve for Stationary Wealth Distribution
Let $\mathbf{A} \in \mathbb{R}^{1 \times S}$ be the vector of total asset mass held by workers in each composite state $(\ell, z')$. Combining the affine wealth policy function with Blanchard-Yaari mortality yields the matrix wealth recursion:
$$\mathbf{A} = \mathbf{A} \mathbf{T}_{\text{surv}} + \boldsymbol{b}_A$$

where the surviving wealth transfer matrix is:
$$\mathbf{T}_{\text{surv}} = (1+r)(1-\kappa)(1-\delta) \mathbf{T}$$

and the exogenous wealth injection vector $\boldsymbol{b}_A \in \mathbb{R}^{1 \times S}$ is (`static/compute_SS_moments.m:55-70`):
$$b_A(\ell, z') = \Lambda(\ell, z') C_{\text{income}}(\ell, z') + \text{ArrivingNB}(\ell, z') a_0$$
with income wealth increment:
$$C_{\text{income}}(\ell, z') = \tilde{y}_\ell + (1+r)\left[ \tilde{w}_\ell z' - \Xi_\ell(z') \right]$$

Since the spectral radius $\rho(\mathbf{T}_{\text{surv}}) < 1$, the exact stationary wealth vector is solved in a single step using the Neumann series matrix inversion:
$$\mathbf{A} = \boldsymbol{b}_A \left(\mathbf{I} - \mathbf{T}_{\text{surv}}\right)^{-1}$$

In MATLAB code, this is computed via efficient sparse left division:
```matlab
A_vec = B_A / (speye(N_states) - T_mat_surv);
a_bar = A_vec ./ Lambda_vec; % per-capita wealth holdings
```

### 3.4 Amenity OLS Regression Audit (`test_amenities_correlation.m`)
To validate the structural amenities $A_\ell$ extracted from the BLP contraction mapping, `test_amenities_correlation.m` estimates an OLS regression against empirical state quality-of-life indicators:
$$A_\ell = \alpha + \beta_1 \text{USDA}_\ell + \beta_2 \text{NAEP}_\ell + \beta_3 \text{EFNA}_\ell + \beta_4 \text{Crime}_\ell + \beta_5 \text{CBP}_\ell + \beta_6 \text{Density}_\ell + \varepsilon_\ell$$

---

## 4. Forward-Backward Dynamic Transition Algorithm

### 4.1 Parameter Path Construction (`build_paths.m`)
The climate change shock develops over a transition period $T_{\text{shock}} = 90$ years (2010 to 2100). To ensure complete numerical convergence to the 2100 steady state, the dynamic simulation horizon is extended to $T_{\text{sim}} = 500$ years with terminal parameter padding.

`build_paths.m` constructs smooth time-series parameter paths $P_t = (w_{\ell,t}, q_{\ell,t}, A_{\ell,t}, \sigma^2_{\varepsilon,\ell,t}, \sigma^2_{\eta,\ell,t})$ across $t \in \{1, \ldots, T_{\text{sim}}\}$:
$$x_{\ell,t} = \begin{cases}
x_{\ell, 2010} + \frac{t-1}{T_{\text{shock}}-1} \left( x_{\ell, 2100} - x_{\ell, 2010} \right) & \text{for } 1 \le t \le T_{\text{shock}} \\
x_{\ell, 2100} & \text{for } T_{\text{shock}} < t \le T_{\text{sim}}
\end{cases}$$

### 4.2 Stage 1: Dynamic Backward Induction (`solve_transition_backward.m`)
Under Perfect Foresight (Scenario A), agents know the full sequence of future spatial parameters $\{P_t\}_{t=1}^{T_{\text{sim}}}$. The value function index matrix array $\boldsymbol{\Xi} \in \mathbb{R}^{J \times N_z \times T_{\text{sim}}}$ is solved backward from $t = T_{\text{sim}}-1$ down to $t=1$:

1. **Terminal Condition**: Set $\boldsymbol{\Xi}_{T_{\text{sim}}} = \boldsymbol{\Xi}_{2100}$ from the terminal static fixed-point solution (`compute_Xi.m`).
2. **Backward Step** (`dynamics/backward_step_Xi.m`): For each $t = T_{\text{sim}}-1, T_{\text{sim}}-2, \ldots, 1$:
   $$\Xi_{\ell,t}(z) = -\frac{1}{\gamma(1+r_t)} \left[ \ln\left(\beta (1+r_t)\right) - \gamma r_t \tilde{w}_{\ell,t} z + A_{\ell,t} + \ln \bar{\lambda} - \frac{1}{\nu} \text{ls}_{\ell,t}(z) \right]$$
   where $\text{ls}_{\ell,t}(z)$ is evaluated using the known future index matrix $\boldsymbol{\Xi}_{t+1}$.

### 4.3 Stage 2: Forward Population Propagation (`simulate_Lambda_forward.m`)
Given the sequence of value indices $\{\boldsymbol{\Xi}_t\}_{t=1}^{T_{\text{sim}}}$, bilateral migration probability matrices $\boldsymbol{\mu}_t$ and sparse state transition matrices $\mathbf{T}_t$ are generated for each period $t$.

Starting from the baseline 2010 population distribution $\boldsymbol{\Lambda}_1 = \boldsymbol{\Lambda}_{2010}$, the joint population distribution mass array $\boldsymbol{\Lambda} \in \mathbb{R}^{J \times N_z \times T_{\text{sim}}}$ is propagated forward sequentially:
$$\boldsymbol{\Lambda}_{t+1} = \boldsymbol{\Lambda}_t \mathbf{T}_t \quad \text{for } t = 1, 2, \ldots, T_{\text{sim}}-1$$

### 4.4 Stage 3: Forward Linear Wealth Propagation (`compute_moments_dyn.m`)
Starting from baseline per-capita wealth holdings $\bar{a}_1(\ell, z) = \bar{a}_{2010}(\ell, z)$, total wealth mass array $\mathbf{A}_{t+1} \in \mathbb{R}^{J \times N_z}$ is updated forward:

$$A_{t+1}(\ell, z') = (1+r)(1-\kappa) \sum_{j=1}^J \sum_{z=1}^{N_z} A_t(j, z) (1-\delta) \mu_{j \to \ell, t}(z) \Pi^\ell_{z, z'} + \Lambda_{t+1}(\ell, z') C_{\text{income}, t}(\ell, z') + \text{ArrivingNB}_{t+1}(\ell, z') a_0$$

Per-capita wealth holdings for the next period are updated via pointwise array division:
$$\bar{a}_{t+1}(\ell, z') = \frac{A_{t+1}(\ell, z')}{\Lambda_{t+1}(\ell, z')}$$

---

## 5. Structural Counterfactuals & Scenario Formulations

### 5.1 Sectoral Damage Functions (Hsiang et al. 2017 Mapping)
Climate damage shocks for 2080–2100 are derived from county-level empirical projections in **Hsiang et al. (2017)** (*Science*), aggregated to state level in `data_cleaning.R`:

- **Wage Damage Shock**:
  $$\text{shock\_wage}_\ell = \frac{\text{labor\_tot}_\ell}{100} + \frac{\text{agri\_yield}_\ell}{100} \cdot \text{agri\_share}_\ell - \frac{\text{coastal}_\ell}{100} (1 - \alpha_{\text{coastal}})$$
- **Amenity Damage Shock**:
  $$\text{shock\_amenity}_\ell = \frac{\text{mortality}_\ell}{100} + \frac{\text{crime\_pct}_\ell}{100}$$
  where crime monetized cost per capita is $\text{cost\_crime\_pc}_\ell = \text{extra\_violent}_\ell \cdot 100{,}000 + \text{extra\_property}_\ell \cdot 2{,}000$, and $\text{crime\_pct}_\ell = (\text{cost\_crime\_pc}_\ell / \text{income\_pc}_\ell) \cdot 100$.
- **Housing Damage Shock**:
  $$\text{shock\_housing}_\ell = \frac{\text{energy\_exp}_\ell}{100} \cdot \text{energy\_share}_\ell - \frac{\text{coastal}_\ell}{100} \alpha_{\text{coastal}}$$

#### Parameters at Terminal Steady State ($t = 2100$):
$$w_{\ell, 2100} = w_{\ell, 2010} \left(1 + \text{shock\_wage}_\ell\right)$$
$$q_{\ell, 2100} = q_{\ell, 2010} \left(1 + \text{shock\_housing}_\ell\right)$$
$$A_{\ell, 2100} = A_{\ell, 2010} - \text{shock\_amenity}_\ell \cdot \tilde{w}_{\ell, 2010}$$

### 5.2 Scenario B: Static Expectations (Myopia)
In `solve_transition_myopia.m`, agents have myopic expectations. At each period $t$, agents assume current parameters $P_t$ will persist forever. Consequently, the value index $\boldsymbol{\Xi}_t^{\text{myopia}}$ is solved as a static fixed point (`compute_Xi.m`) using current parameter vector $P_t$, without backward induction. Forward population and wealth propagation use this myopic index sequence.

### 5.3 Scenario C & H: Climate Uncertainty (Variance Shocks)
Under climate uncertainty (`run_scenario('perfect_foresight', P_unc2)`), climate change increases transitory labor income volatility $\sigma^2_{\eta,\ell,t}$ proportionally to local wage damages:
$$\sigma^2_{\eta,\ell,t} = \sigma^2_{\eta,\ell,2010} \left(1 + 50 \cdot \|\text{shock\_wage}_\ell\| \cdot \frac{t-1}{T_{\text{shock}}-1}\right)$$

In `compute_Xi.m:46`, higher variance increases flow risk penalty $\hat{B}_{1,\ell,t} = \gamma \kappa \tilde{y}_{\ell,t} - \frac{1}{2}(\kappa \gamma \tilde{w}_{\ell,t})^2 \sigma^2_{\eta,\ell,t}$, inducing strong spatial precautionary savings and out-migration from high-variance states.

### 5.4 Scenario D: Hand-to-Mouth Structural Counterfactual (No Saving)
In the no-saving counterfactual (`no_saving/`), agents are constrained to consume their disposable income each period ($a_t \equiv 0, c_t = \tilde{y}_{\ell,t}(z, \eta)$). Precautionary wealth accumulation is completely shut off.

The value function index $\Xi_\ell^{\text{NS}}(z)$ satisfies the non-linear functional equation (`no_saving/compute_Xi_no_saving.m:60-80`):
$$\exp\left(-\gamma \Xi_\ell^{\text{NS}}(z)\right) = \exp\left(F_\ell(z)\right) + \beta \Omega_\ell(z)$$

where:
$$F_\ell(z) = -\gamma \tilde{y}_\ell - \gamma \tilde{w}_\ell z - A_\ell + \frac{1}{2}\left(\gamma \tilde{w}_\ell\right)^2 \sigma^2_{\eta,\ell}$$
$$\Omega_\ell(z) = \bar{\lambda} \left[ \frac{1}{J} \sum_{m=1}^J \exp\left(-\nu \psi_{\ell m} - \nu \ln \mathbb{E}_{z'|z,m}\left[ \exp\left(-\gamma \Xi_m^{\text{NS}}(z')\right) \right]\right) \right]^{-1/\nu}$$

Individual Certainty Equivalent welfare collapses to pure non-asset index $\text{CE}_{\ell,t}^{\text{NS}}(z) = \Xi_{\ell,t}^{\text{NS}}(z)$.

### 5.5 Scenarios E–H: Difference-in-Differences Autarkic Counterfactuals
`run_did.m` implements a Difference-in-Differences (DiD) structural identification framework comparing Treatment and Control arms:
- **Treatment Arm**: Exposes economy to climate shocks under migration friction $\boldsymbol{\psi}_{\text{variant}}$.
- **Control Arm**: Holds climate parameters constant at 2010 baseline under migration friction $\boldsymbol{\psi}_{\text{variant}}$.
- **Pure DiD Isolator**:
  $$\Delta^{\text{DiD}} X_t = X_t^{\text{Treat}} - X_t^{\text{Ctrl}}$$

- **Scenario E (No Migration Autarky)**: Prohibits inter-state moving by setting $\psi_{j\ell} = 10^8$ for all $j \neq \ell$.
- **Scenario F (Small Migration)**: Increases moving costs by 10% ($\boldsymbol{\psi}_{\text{new}} = 1.1 \times \boldsymbol{\psi}$).
- **Scenario G (Big Migration)**: Decreases moving costs by 10% ($\boldsymbol{\psi}_{\text{new}} = 0.9 \times \boldsymbol{\psi}$).
- **Scenario H (Variance Shocks + Autarky)**: Combines variance shocks with infinite moving friction ($\psi_{j\ell} = 10^8$).

---

## 6. Comprehensive Figure-to-Math & Code Mapping Matrix

The matrix below maps all 30 figure environments appearing in `Draft/leon_memoire_brainstorming.tex` to their generating MATLAB functions, script line numbers, output image file paths, and core mathematical formulas.

| LaTeX Figure Panel | Description | Output PNG File Path | Generating MATLAB Function & Lines | Scenario Runner Call | Core Mathematical Formula / Definition |
|---|---|---|---|---|---|
| **Fig 1(a)** (`fig:damage_maps`) | Wage Damage Map | `Output/graph_perfect_foresight/damage_maps/damage_maps_damage_wage.png` | `plot_damage_maps.m:25-50` | `run_scenario.m` | $\text{shock\_wage}_\ell = \frac{\bar L_\ell}{100} + \frac{\bar Y_\ell^{\text{agri}}}{100}s_\ell^{\text{agri}} - \frac{\bar C_\ell}{100}(1-\alpha_{\text{coastal}})$ |
| **Fig 1(b)** (`fig:damage_maps`) | Amenity Damage Map | `Output/graph_perfect_foresight/damage_maps/damage_maps_damage_amenity.png` | `plot_damage_maps.m:51-75` | `run_scenario.m` | $\text{shock\_amenity}_\ell = \frac{\text{mortality}_\ell}{100} + \frac{\text{cost\_crime\_pc}_\ell}{\text{income\_pc}_\ell}$ |
| **Fig 1(c)** (`fig:damage_maps`) | Housing Damage Map | `Output/graph_perfect_foresight/damage_maps/damage_maps_damage_housing.png` | `plot_damage_maps.m:76-100` | `run_scenario.m` | $\text{shock\_housing}_\ell = \frac{\bar E_\ell}{100} s_\ell^{\text{energy}} - \frac{\bar C_\ell}{100} \alpha_{\text{coastal}}$ |
| **Fig 2(a)** (`trans_dyn_pf`) | Earnings (PF) | `Output/graph_perfect_foresight/dynamic_moments/dynamic_moments_national_earnings.png` | `plot_dynamic_moments.m:30-55` | `run_scenario('perfect_foresight', ...)` | $Y_t = \sum_\ell \tilde{w}_{\ell,t} \left( L_{\ell,t} + \sum_z z \Lambda_t(\ell,z) \right)$ |
| **Fig 2(b)** (`trans_dyn_pf`) | Mobility (PF) | `Output/graph_perfect_foresight/dynamic_moments/dynamic_moments_national_mobility.png` | `plot_dynamic_moments.m:56-80` | `run_scenario('perfect_foresight', ...)` | $m_t = \sum_{\ell=1}^J \sum_z \Lambda_t(\ell, z) \left[ 1 - \mu_{\ell \to \ell, t}(z) \right]$ |
| **Fig 2(c)** (`trans_dyn_pf`) | Wealth (PF) | `Output/graph_perfect_foresight/dynamic_moments/dynamic_moments_national_wealth.png` | `plot_dynamic_moments.m:81-105` | `run_scenario('perfect_foresight', ...)` | $W_t = \sum_{\ell=1}^J \sum_{z=1}^{N_z} \bar a_t(\ell, z) \Lambda_t(\ell, z)$ |
| **Fig 2(d)** (`trans_dyn_pf`) | W/Y Ratio (PF) | `Output/graph_perfect_foresight/dynamic_moments/dynamic_moments_national_wy.png` | `plot_dynamic_moments.m:106-130` | `run_scenario('perfect_foresight', ...)` | $(W/Y)_t = W_t / Y_t$ |
| **Fig 3(a)** (`trans_map_pf`) | $\Delta$ Welfare (Level) | `Output/graph_perfect_foresight/transition_maps/transition_maps_delta_mean_welfare_ce.png` | `plot_transition_maps.m:35-60` | `run_scenario('perfect_foresight', ...)` | $\Delta \text{CE}_\ell = \overline{\text{CE}}_{\ell, 2100} - \overline{\text{CE}}_{\ell, 2010}, \quad \text{CE}_\ell(z) = \kappa \bar a_\ell(z) + \Xi_\ell(z)$ |
| **Fig 3(b)** (`trans_map_pf`) | $\Delta$ W/Y (Level) | `Output/graph_perfect_foresight/transition_maps/transition_maps_wealth_income_ratio_changes.png` | `plot_transition_maps.m:61-85` | `run_scenario('perfect_foresight', ...)` | $\Delta (W/Y)_\ell = (W/Y)_{\ell, 2100} - (W/Y)_{\ell, 2010}$ |
| **Fig 3(c)** (`trans_map_pf`) | Pop $\Delta$ (%) | `Output/graph_perfect_foresight/transition_maps/transition_maps_local_pop_changes_pct.png` | `plot_transition_maps.m:86-110` | `run_scenario('perfect_foresight', ...)` | $\Delta L_\ell\% = 100 \times (L_{\ell, 2100} - L_{\ell, 2010}) / L_{\ell, 2010}$ |
| **Fig 3(d)** (`trans_map_pf`) | Pop $\Delta$ (pp) | `Output/graph_perfect_foresight/transition_maps/transition_maps_pop_changes_pp.png` | `plot_transition_maps.m:111-135` | `run_scenario('perfect_foresight', ...)` | $\Delta p_\ell = 100 \times (p_{\ell, 2100} - p_{\ell, 2010}), \quad p_\ell = L_\ell / \sum_m L_m$ |
| **Fig 4(a)** (Ineq PF) | Welfare Ineq (Map) | `Output/graph_perfect_foresight/welfare_wealth_inequality_maps/..._welfare_inequality.png` | `plot_welfare_variance.m:25-50` | `run_scenario('perfect_foresight', ...)` | $\text{Var}_z(\text{CE}_\ell(z)) = \sum_z \Lambda(\ell,z) [\text{CE}_\ell(z) - \overline{\text{CE}}_\ell]^2 / L_\ell$ |
| **Fig 4(b)** (Ineq PF) | Wealth Ineq (Map) | `Output/graph_perfect_foresight/welfare_wealth_inequality_maps/..._wealth_inequality.png` | `plot_welfare_variance.m:51-75` | `run_scenario('perfect_foresight', ...)` | $\text{Var}_z(\bar a_\ell(z)) = \sum_z \Lambda(\ell,z) [\bar a_\ell(z) - \bar a_\ell]^2 / L_\ell$ |
| **Fig 4(c)** (Ineq PF) | Welfare Ineq (Trans) | `Output/graph_perfect_foresight/welfare_transition/welfare_transition_welfare_inequality.png` | `plot_welfare_transition.m:25-50` | `run_scenario('perfect_foresight', ...)` | $\text{Ineq}_t^{\text{CE}} = \sum_\ell p_{\ell,t} [\overline{\text{CE}}_{\ell,t} - \overline{\text{CE}}_t^{\text{agg}}]^2$ |
| **Fig 4(d)** (Ineq PF) | Wealth Ineq (Trans) | `Output/graph_perfect_foresight/welfare_transition/welfare_transition_wealth_inequality.png` | `plot_welfare_transition.m:51-75` | `run_scenario('perfect_foresight', ...)` | $\text{Ineq}_t^{\text{W}} = \sum_\ell p_{\ell,t} [\bar a_{\ell,t} - \bar a_t^{\text{agg}}]^2$ |
| **Fig 5(a-d)** (Myopia) | Dynamic Moments (Myopic) | `Output/graph_myopia/dynamic_moments/dynamic_moments_national_*.png` | `plot_dynamic_moments.m:30-130` | `run_scenario('myopia', ...)` | Same as Fig 2, but $\boldsymbol{\Xi}_t$ solved as static fixed point at each $t$ |
| **Fig 6(a-d)** (Myopia) | Transition Maps (Myopic) | `Output/graph_myopia/transition_maps/...` | `plot_transition_maps.m:35-135` | `run_scenario('myopia', ...)` | Same as Fig 3, evaluated under myopic expectation path |
| **Fig 7(a-d)** (Myopia) | Inequality Maps (Myopic) | `Output/graph_myopia/welfare_wealth_inequality_maps/...` | `plot_welfare_variance.m:25-75` | `run_scenario('myopia', ...)` | Same as Fig 4, evaluated under myopic expectation path |
| **Fig 8** (`variance_shock_map`) | Variance Shock Map | `Output/graph_uncertainty_variance_shocks/variance_map.png` | `plot_damage_maps.m:105-130` | `run_scenario('perfect_foresight', P_unc2)` | $\Delta \sigma^2_{\eta,\ell} = 50 \cdot \|\text{shock\_wage}_\ell\|$ |
| **Figs 9-11** (Uncertainty) | Moments/Maps (Uncertainty) | `Output/graph_uncertainty_variance_shocks/...` | `plot_*` utilities | `run_scenario('perfect_foresight', P_unc2)` | Evaluated with progressive variance shock path $\sigma^2_{\eta,\ell,t}$ |
| **Figs 12-14** (No Saving) | Moments/Maps (No Saving) | `Output/graph_no_saving_struct/...` | `plot_*` utilities | `run_scenario('perfect_foresight', P, 'no_saving_struct')` | Hand-to-mouth recursion: $e^{-\gamma \Xi_\ell^{\text{NS}}(z)} = e^{F_\ell(z)} + \beta \Omega_\ell(z)$ |
| **Fig 15(a-d)** (No Mig DiD) | Dynamic Diffs (No Mig) | `Output/graph_no_migration_did/did_dynamic_moments_diff_*.png` | `plot_did_dynamic_moments.m:30-100` | `run_did(P, 'no_migration')` | Pure DiD: $\Delta^{\text{DiD}} X_t = X_t^{\text{treat}} - X_t^{\text{ctrl}}$, with $\psi_{j\ell}^{\text{ban}} = 10^8$ |
| **Fig 16(a-d)** (No Mig DiD) | Transition Maps (No Mig) | `Output/graph_no_migration_did/transition_maps/...` | `plot_transition_maps.m:35-135` | `run_did(P, 'no_migration')` | Transition maps comparing autarkic treatment vs autarkic control |
| **Fig 17(a-b)** (No Mig DiD) | Inequality Maps (No Mig) | `Output/graph_no_migration_did/welfare_wealth_inequality_maps/...` | `plot_welfare_variance.m:25-75` | `run_did(P, 'no_migration')` | Inequality maps under autarkic Difference-in-Differences |
| **Figs 18-20** (Small Mig) | Moments/Maps (Small Mig) | `Output/graph_small_migration_did/...` | `plot_*` utilities | `run_did(P, 'small_migration')` | DiD evaluated with $\boldsymbol{\psi}_{\text{new}} = 1.1 \times \boldsymbol{\psi}$ |
| **Figs 21-23** (Big Mig) | Moments/Maps (Big Mig) | `Output/graph_big_migration_did/...` | `plot_*` utilities | `run_did(P, 'big_migration')` | DiD evaluated with $\boldsymbol{\psi}_{\text{new}} = 0.9 \times \boldsymbol{\psi}$ |
| **Figs 24-26** (Var No Mig) | Moments/Maps (Var No Mig) | `Output/graph_uncertainty_variance_shocks_no_migration_did/...` | `plot_*` utilities | `run_did(P_unc3, 'no_migration')` | DiD combining variance shocks and autarkic migration ban |
| **Fig 27** (Welfare Audit) | Aggregate CE Loss Across Scenarios | `Output/graph_welfare_audit/welfare_audit_all_scenarios/...` | `compare_all_welfare.m:40-90` | `main.m` orchestrator | $\text{CE}^{\text{agg}}_t = \sum_{\ell=1}^J \sum_{z=1}^{N_z} \Lambda_t(\ell,z) \left( \kappa \bar a_t(\ell,z) + \Xi_{t}(\ell,z) \right)$ |
| **Fig 28(a-b)** (Macro Audit) | Macro Diffs Across Scenarios | `Output/graph_macro_audit/macro_earning/...` & `macro_wealth/...` | `compare_all_macro.m:35-110` | `main.m` orchestrator | Cross-scenario time series overlay of earnings, wealth, W/Y, mobility |
| **Fig 29(a-c)** (U-Shape) | U-Shape Dynamics at $t=500$ | `Output/graph_ushape/ushape_t500_california.png` (CA, FL, WA) | `plot_ushape_t500_split.m:45-120` | `main.m` orchestrator | Trajectory of population share $p_{\ell,t}$ demonstrating initial dip and long-run recovery |
| **Fig 30** (Florida Comp) | Florida Population Trajectory | `Output/graph_perfect_foresight/florida_dynamics_comparison/...` | `plot_florida_dynamics.m:30-75` | `main.m` orchestrator | Florida population share comparison across foresight & expectations |

---

## 7. Mathematical Construction of Output Maps

### 7.1 Spatial Level and Difference Maps
`plot_transition_maps.m` computes 4 long-run steady-state transition difference choropleth maps comparing 2100 ($t=90$) against baseline 2010 ($t=1$):

1. **Welfare Change Map ($\Delta \text{CE}_\ell$)**:
   $$\Delta \text{CE}_\ell = \overline{\text{CE}}_{\ell, 2100} - \overline{\text{CE}}_{\ell, 2010}$$
   where state mean Certainty Equivalent is $\overline{\text{CE}}_\ell = \frac{1}{L_\ell} \sum_{z=1}^{N_z} \Lambda(\ell, z) \left[ \kappa \bar{a}_\ell(z) + \Xi_\ell(z) \right]$.
2. **Wealth-to-Earnings Ratio Change Map ($\Delta (W/Y)_\ell$)**:
   $$\Delta (W/Y)_\ell = \frac{W_{\ell, 2100}}{Y_{\ell, 2100}} - \frac{W_{\ell, 2010}}{Y_{\ell, 2010}}$$
3. **Relative Population Change Map ($\Delta L_\ell\%$)**:
   $$\Delta L_\ell\% = 100 \times \frac{L_{\ell, 2100} - L_{\ell, 2010}}{L_{\ell, 2010}}$$
4. **Absolute Population Share Change Map ($\Delta p_\ell$)**:
   $$\Delta p_\ell = 100 \times \left( p_{\ell, 2100} - p_{\ell, 2010} \right), \quad p_\ell = \frac{L_\ell}{\sum_{m=1}^J L_m}$$

### 7.2 Shift-Share Decomposition of Wealth-to-Earnings Ratio (`plot_spatial_decomposition.m`)
To isolate whether state wealth-to-earnings ratios $(W/Y)_\ell$ change due to individual precautionary savings or spatial population re-sorting, `plot_spatial_decomposition.m` executes a exact shift-share decomposition:

Let $\pi_\ell(z, t) = \frac{\Lambda_t(\ell, z)}{\sum_{z'} \Lambda_t(\ell, z')}$ be the fraction of state $\ell$'s population in productivity state $z$. The state wealth-to-earnings ratio is:
$$\left(\frac{W}{Y}\right)_{\ell,t} = \sum_{z=1}^{N_z} \pi_\ell(z, t) \left( \frac{\bar{a}_{\ell,t}(z)}{\bar{y}_{\ell,t}} \right)$$

The total change between 2010 and 2100 decomposes into two additive margins:
$$\Delta \left(\frac{W}{Y}\right)_\ell = \underbrace{\sum_{z=1}^{N_z} \pi_\ell(z, 2010) \cdot \Delta \left( \frac{\bar{a}_\ell(z)}{\bar{y}_\ell} \right)}_{\text{Intensive Precautionary Savings Margin}} + \underbrace{\sum_{z=1}^{N_z} \left( \frac{\bar{a}_{\ell,2100}(z)}{\bar{y}_{\ell,2100}} \right) \cdot \Delta \pi_\ell(z)}_{\text{Compositional Migration Sorting Margin}}$$

### 7.3 Spatial Inequality Metrics (`plot_welfare_variance.m`)
- **Within-State Worker Welfare Variance Map**:
  $$\text{Var}_z\left(\text{CE}_\ell(z)\right) = \frac{1}{L_\ell} \sum_{z=1}^{N_z} \Lambda(\ell,z) \left[ \text{CE}_\ell(z) - \overline{\text{CE}}_\ell \right]^2$$
- **National Cross-State Inequality Transition**:
  $$\text{Ineq}_t^{\text{CE}} = \sum_{\ell=1}^J p_{\ell,t} \left[ \overline{\text{CE}}_{\ell,t} - \overline{\text{CE}}_t^{\text{agg}} \right]^2$$

---

## 8. Mathematical Construction of Dynamic Moment Graphs

### 8.1 National Aggregate Macro Trajectories (`plot_dynamic_moments.m`)
`plot_dynamic_moments.m` plots 4 core national macroeconomic time series across $t \in \{1, \ldots, T_{\text{sim}}\}$:

1. **National Aggregate Annual Earnings ($Y_t$)**:
   $$Y_t = \sum_{\ell=1}^J \tilde{w}_{\ell,t} \left( L_{\ell,t} + \sum_{z=1}^{N_z} z \Lambda_t(\ell,z) \right)$$
2. **National Aggregate Wealth Mass ($W_t$)**:
   $$W_t = \sum_{\ell=1}^J \sum_{z=1}^{N_z} \bar{a}_t(\ell, z) \Lambda_t(\ell, z)$$
3. **National Aggregate Wealth-to-Earnings Ratio ($(W/Y)_t$)**:
   $$(W/Y)_t = \frac{W_t}{Y_t}$$
4. **National Gross Annual Mobility Rate ($m_t$)**:
   $$m_t = \frac{1}{\sum_{\ell=1}^J L_{\ell,t}} \sum_{\ell=1}^J \sum_{z=1}^{N_z} \Lambda_t(\ell, z) \left[ 1 - \mu_{\ell \to \ell, t}(z) \right]$$

### 8.2 Certainty Equivalent Welfare Metric
The aggregate Certainty Equivalent welfare level in period $t$ is:
$$\text{CE}_t^{\text{agg}} = \sum_{\ell=1}^J \sum_{z=1}^{N_z} \Lambda_t(\ell,z) \left[ \kappa \bar{a}_t(\ell,z) + \Xi_{\ell,t}(z) \right]$$

---

## 9. Difference-in-Differences & Counterfactual Analytics

### 9.1 DiD Identification Strategy (`run_did.m`)
The structural Difference-in-Differences counterfactual isolates climate damage losses under counterfactual spatial mobility regimes:
$$\Delta^{\text{DiD}} X_t = X_t^{\text{Treat}}\left(P_{\text{climate}}, \boldsymbol{\psi}_{\text{variant}}\right) - X_t^{\text{Ctrl}}\left(P_{\text{no\_climate}}, \boldsymbol{\psi}_{\text{variant}}\right)$$

### 9.2 Spatial Insurance Value Decomposition
By comparing scenario outcomes, we quantify two fundamental economic mechanisms:

1. **Value of Spatial Insurance (Geographic Mobility)**:
   $$\Delta \text{CE}_{\text{insurance}, t} = \text{CE}_t^{\text{agg}}(\text{Scenario A: Free Migration}) - \text{CE}_t^{\text{agg}}(\text{Scenario E: Autarky})$$
2. **Value of Precautionary Savings**:
   $$\Delta \text{CE}_{\text{savings}, t} = \text{CE}_t^{\text{agg}}(\text{Scenario A: Saving Allowed}) - \text{CE}_t^{\text{agg}}(\text{Scenario D: Hand-to-Mouth})$$

---

## 10. Long-Run Dynamics, Spatial Sorting ($t=500$), and Cross-Scenario Welfare Audit

### 10.1 Very Long-Run Dynamics ($t=500$) & U-Shape Trajectory (`plot_ushape_t500_split.m`)
When tracking state population shares $p_{\ell,t} = L_{\ell,t} / \sum_m L_{m,t}$ up to $t=500$, heavily exposed states (such as Florida, California, and Washington) exhibit a distinct non-monotonic **U-shaped population trajectory**:

```
Population Share p_{ell, t}
  ^
  |  * (2010 Baseline)
  |   \                                   /---* Long-run steady state (t=500)
  |    \                                 /
  |     \ (Direct Climate Damage Dip)   /  (Capital Re-accumulation & Rent Adjustment)
  |      \                             /
  |       \                           /
  |        \                         /
  |         +-----------------------+
  |                   t = 90 (2100 Peak Damage)
  +---------------------------------------------------------> Time t
```

- **Phase 1: Out-Migration Dip ($t \in [1, 90]$)**: Rapid climate wage degradation $w_{\ell,t}$ and amenity destruction $A_{\ell,t}$ induce out-migration.
- **Phase 2: Capital Re-accumulation & Spatial Re-equilibration ($t \in [90, 500]$)**: Declining local population lowers housing rent $q_{\ell,t}$, while long-run capital accumulation restores local per-capita wealth $\bar{a}_\ell(z)$. Higher effective returns attract population back, leading to a long-run recovery.

### 10.2 Cross-Scenario Welfare Audit (`compare_all_welfare.m`)
`compare_all_welfare.m` synthesizes aggregate Certainty Equivalent welfare losses across all 9 scenarios.

The welfare loss ranking demonstrates:
1. **Uncertainty + Autarky (Scenario H)** suffers the worst lifetime welfare degradation.
2. **Hand-to-Mouth (Scenario D)** reveals that shutting off savings amplifies climate damage by depriving households of self-insurance.
3. **Perfect Foresight (Scenario A)** yields lower welfare losses than **Myopia (Scenario B)**, confirming that forward-looking expectation allows households to pre-accumulate precautionary wealth ahead of climate shock arrival.

---

## 11. Technical Appendix: Numerical Algorithms

### 11.1 Rouwenhorst Discretization Formulas (`helpers/discretize_normal.m`)
For AR(1) process $z' = \rho z + \varepsilon, \varepsilon \sim \mathcal{N}(0, \sigma^2_\varepsilon)$, the $N_z$-state Markov chain transition matrix $\mathbf{\Pi}_{N_z}$ is constructed recursively:
$$\mathbf{\Pi}_2 = \begin{bmatrix} p & 1-p \\ 1-p & p \end{bmatrix}, \quad p = \frac{1 + \rho}{2}$$
$$\mathbf{\Pi}_N = p \begin{bmatrix} \mathbf{\Pi}_{N-1} & \mathbf{0} \\ \mathbf{0}^\top & 0 \end{bmatrix} + (1-p) \begin{bmatrix} \mathbf{0} & \mathbf{\Pi}_{N-1} \\ 0 & \mathbf{0}^\top \end{bmatrix} + (1-p) \begin{bmatrix} \mathbf{0}^\top & 0 \\ \mathbf{\Pi}_{N-1} & \mathbf{0} \end{bmatrix} + p \begin{bmatrix} 0 & \mathbf{0}^\top \\ \mathbf{0} & \mathbf{\Pi}_{N-1} \end{bmatrix}$$

Divisor normalization ensures row-stochasticity $\sum_{z'} \Pi(z, z') = 1$.

### 11.2 Log-Sum-Exp Numerical Stability Protocols
To prevent floating-point underflow/overflow when computing exponents in value recursions and Weibull probabilities, all sums are evaluated using the Log-Sum-Exp identity:
$$\ln \sum_{i=1}^M \exp(x_i) = x^* + \ln \sum_{i=1}^M \exp(x_i - x^*), \quad x^* = \max_{i} x_i$$
