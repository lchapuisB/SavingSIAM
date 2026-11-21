import re

tex_content = r"""\documentclass[aspectratio=169, 10pt]{beamer}
\usetheme{Berlin}
\useinnertheme{rounded}
\useoutertheme{infolines}
\usefonttheme{serif}

\definecolor{forestgreen}{RGB}{46, 90, 60}
\definecolor{lightleaf}{RGB}{240, 245, 238}
\definecolor{moss}{RGB}{138, 154, 91}

\setbeamercolor{frametitle}{bg=lightleaf, fg=forestgreen}
\setbeamercolor{structure}{fg=forestgreen}
\setbeamercolor{normal text}{fg=black, bg=white}
\setbeamercolor{block title}{bg=forestgreen!15, fg=forestgreen!70!black}
\setbeamercolor{block body}{bg=white}

\usepackage[most]{tcolorbox}
\usepackage{tikz}

\tcbset{
  leafbox/.style={
    enhanced,
    colback=lightleaf!80!white,
    colframe=forestgreen!60!black,
    coltitle=forestgreen!70!black,
    fonttitle=\bfseries,
    boxrule=0pt,
    frame hidden,
    borderline west={3pt}{0pt}{forestgreen!70!black},
    sharp corners,
    title style={colback=forestgreen!10!white},
  }
}

% Packages standards
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{subcaption} 
\usepackage{bbm} 
\newcommand{\ind}[1]{\mathbbm{1}_{\{#1\}}}

% Infos de la présentation
\title[Climate Adaptation: Saving \& Migration]{\textbf{Should I Save and Should I Go?}}
\subtitle{A Spatial Model of Climate Change Adaptation with Saving and Migration}
\author{Léon Chapuis}
\institute{Paris School of Economics / ENS-PSL}
\date{Winter 2026}

\begin{document}

\begin{frame}
  \titlepage
\end{frame}

\section{Introduction}

\begin{frame}{Climate change requires individual adaptation}
    \begin{itemize}
        \item Mitigation policies have repeatedly fallen short of international targets (e.g., missed Paris Agreement milestones).
        \item Extreme weather events, such as the intense forest fires this summer, are becoming more frequent.
        \item \textbf{Implication:} Households can no longer rely solely on global mitigation; \textit{individual adaptation} is increasingly critical.
    \end{itemize}
\end{frame}

\begin{frame}{Climate change impacts are spatially heterogeneous}
    \begin{itemize}
        \item Climate damages do not affect all regions equally (e.g., Burke et al., Hsiang et al. for the US).
        \item Some states suffer severe productivity drops and amenity losses, while others may experience milder effects or even gains.
        \item \textbf{Implication:} \textit{Migration} emerges as a key individual mechanism to escape localized climate shocks.
    \end{itemize}
\end{frame}

\begin{frame}{Saving is a crucial but often omitted adaptation mechanism}
    \begin{itemize}
        \item \textbf{The Dimensionality Curse:} Spatial macroeconomic models often omit wealth accumulation due to the intractable state space when combining forward-looking migration and saving.
        \item \textbf{Distributional Impacts:} Wealth inequality strongly shapes the ability to adapt.
        \item \textbf{Literature:}
        \begin{itemize}
            \item Migration models (Caliendo et al. 2019) often assume hand-to-mouth consumers.
            \item Recent efforts (Sun 2025, Giannone et al. 2023) tackle this but rely on heavy computational methods (e.g., neural networks, OLG).
        \end{itemize}
        \item \textbf{Implication:} We need a tractable structural framework to jointly evaluate saving and migration.
    \end{itemize}
\end{frame}

\begin{frame}{How do saving and migration jointly shape climate adaptation?}
    \begin{tcolorbox}[leafbox, title=Research Question]
        What are the aggregate and distributional impacts of climate change when individuals can both save and migrate, and how do these two adaptive margins interact?
    \end{tcolorbox}
    \vspace{0.3cm}
    \textbf{Approach:}
    \begin{itemize}
        \item Build on the macro-finance spatial framework of \textbf{Dvorkin \& Greaney (2025)}.
        \item \textbf{Extensions:} 
        \begin{enumerate}
            \item Introduce a deterministic climate change trajectory (2010--2100).
            \item Explicitly solve for the transition path.
            \item Disentangle migration costs from local amenities.
        \end{enumerate}
    \end{itemize}
\end{frame}

\section{Model}

\begin{frame}{Economic Environment: Geography, Endowments, and Risk}
    \begin{enumerate}
        \item \textbf{Spatial Structure:} 
        \begin{itemize}
            \item $J$ locations (US States) with migration frictions $\psi_{j,\ell}$.
            \item Idiosyncratic location preferences $\varepsilon$ drawn from a Weibull distribution.
        \end{itemize}
        \vspace{0.2cm}
        \item \textbf{Income Processes:} 
        \begin{itemize}
            \item Workers supply efficient units of labor subject to two uninsurable risks: a persistent shock $z'$ and a transitory shock $\eta'$.
            \item Wage income: $w_{j,t} (1 + z' + \eta')$
        \end{itemize}
        \vspace{0.2cm}
        \item \textbf{Incomplete Markets:} 
        \begin{itemize}
            \item No state-contingent claims. Workers self-insure via a risk-free bond (saving/borrowing).
        \end{itemize}
        \vspace{0.2cm}
        \item \textbf{Climate Change:} 
        \begin{itemize}
            \item Exogenous shocks affect local fundamentals: wages $D^w_{j,t}$, amenities $D^A_{j,t}$, and housing costs $D^q_{j,t}$.
        \end{itemize}
    \end{enumerate}
\end{frame}

\begin{frame}{Consumer Program \& Timing}
    \textbf{Intra-Period Timing:}
    \begin{enumerate}
        \item Enter period in location $j$ with assets $a$ and past productivity $z$.
        \item Draw preference shocks $\varepsilon$ and observe climate fundamentals.
        \item \textbf{Migration Choice:} Choose destination $\ell$.
        \item Relocate and draw new income shocks $z'$ and $\eta'$.
        \item \textbf{Financial Choice:} Choose consumption $c_{\ell,t}$ and next-period savings $a_{\ell,t+1}$.
    \end{enumerate}
    \vspace{0.3cm}
    \begin{tcolorbox}[leafbox, title=Household Value Function (CARA Utility)]
        $V_{j,t}(a, z, \varepsilon) = \max_{\ell} \left\{ \varepsilon_\ell e^{\psi_{j\ell}} \mathbb{E}_{\eta^\prime, z^\prime | z, \ell} \left[ \max_{c, a'} \left\{ -\frac{1}{\gamma} e^{-\gamma c_{\ell,t} - \tilde{A}_{\ell,t}} + \beta \mathbb{E}_{\varepsilon^\prime} \left[ V_{\ell,t+1}(a', z^\prime, \varepsilon^\prime) \right] \right\} \right] \right\}$
    \end{tcolorbox}
\end{frame}

\begin{frame}{Proposition 1: Tractable Closed-Form Solutions}
    \begin{tcolorbox}[leafbox, title=Optimal Asset and Consumption Choices]
        $c_{\ell,t}(a, z^\prime, \eta^\prime) = \Xi_{\ell,t}(z^\prime) + \frac{r}{1+r} \left[ \tilde{y}_{\ell,t} + \tilde{w}_{\ell,t}\eta^\prime + a \right]$ \\
        $a_{\ell,t+1}(a, z^\prime, \eta^\prime) = \left[ \tilde{y}_{\ell,t} + \tilde{w}_{\ell,t}\eta^\prime + a \right] + (1+r) \left( \tilde{w}_{\ell,t}z^\prime - \Xi_{\ell,t}(z^\prime) \right)$
    \end{tcolorbox}
    \begin{tcolorbox}[leafbox, title=Mobility Choice \& Continuation Value]
        $\mu_{j\ell,t}(z) = \frac{ \exp\left(\nu \mathbb{B}_{1,j\ell,t}\right) \left( \mathbb{E}_{z^\prime | z, \ell} \left[ e^{-\gamma \Xi_{\ell,t}(z^\prime)} \right] \right)^{-\nu} }{ \sum_{m} \exp\left(\nu \mathbb{B}_{1,jm,t}\right) \left( \mathbb{E}_{z^\prime | z, m} \left[ e^{-\gamma \Xi_{m,t}(z^\prime)} \right] \right)^{-\nu} }$
    \end{tcolorbox}
    
    \begin{itemize}
        \item $\Xi_{\ell,t}(z^\prime)$ captures all forward-looking components (future income, precautionary motives, spatial option value).
        \item \textbf{Superb tractability feature:} Due to CARA preferences, there are \textit{no wealth effects} on migration. The state space collapses from $(a, z, \ell)$ to just $(z, \ell)$.
    \end{itemize}
\end{frame}

\begin{frame}{Saving Motives \& Wealth Accumulation}
    \begin{equation}
        \Delta a_\ell^{\text{Saving Rule}}(z) \propto \sum_{s=0}^{\infty} \frac{1}{(1+r)^s} \left[ \underbrace{- z \rho^s (1-\rho) w \Delta^{SS} D^w_\ell}_{\text{Damage-Level Dampening}} + \underbrace{\frac{\gamma}{2} \Delta^{SS} \sigma_{C, \ell}^2}_{\Delta \text{Precautionary Risk}} + \underbrace{\frac{1}{\gamma \nu} \Delta^{SS} \ln \mu_{\ell \ell}(z)}_{\Delta \text{Spatial Option Value}} \right]
        \label{eq:intensive_margin_decomposition}
    \end{equation}
    
    \begin{equation}
        \Delta^{SS} \bar{a}_\ell = \sum_z \pi_\ell^{\text{init}}(z) \big[\Delta a_\ell^{\text{Saving Rule}}(z) + \Delta a_\ell^{\text{Transition}}(z) \big] + \underbrace{\sum_z \big[\pi_\ell^{\text{term}}(z) - \pi_\ell^{\text{init}}(z) \big] a_\ell^{\text{term}}(z)}_{\text{Compositional Sorting}}
        \label{eq:wealth_decomposition}
    \end{equation}
    
    \begin{itemize}
        \item \textbf{Permanent income} and \textbf{amenity smoothing} effects exist \textit{only} on the transition path ($\Delta a_\ell^{\text{Transition}}$).
        \item Long-run wealth changes ($\Delta^{SS} \bar{a}_\ell$) are driven by risk changes and compositional sorting of productivity types.
    \end{itemize}
\end{frame}

\section{Calibration \& Damages}

\begin{frame}{Calibration matches US spatial and economic moments}
    \begin{itemize}
        \item Migration costs are parameterized using a gravity structure: $\psi_{j,\ell} = \psi_0 + \psi_1 \ln(\text{distance}_{j,\ell})$.
        \item The variance of preference shocks $\nu$ governs migration elasticity (set to match inter-state flows).
    \end{itemize}
    \vspace{0.3cm}
    \begin{table}[h]
        \centering
        \small
        \begin{tabular}{clc}
            \toprule
            \textbf{Parameter} & \textbf{Description} & \textbf{Value} \\
            \midrule
            $r$ & Interest rate & 0.04 \\
            $\gamma$ & Absolute risk aversion & 2.0 \\
            $\delta$ & Death probability & 0.025 \\
            $\rho$ & Persistence of $z$ shock & 0.91 \\
            $\nu$ & Preference shock shape (migration elasticity) & 0.5 \\
            \bottomrule
        \end{tabular}
    \end{table}
\end{frame}

\begin{frame}{Climate damages affect wages, amenities, and housing costs}
    \begin{itemize}
        \item Damages mapped from Hsiang et al. (2017) across 6 channels: agriculture, crime, coastal storms, energy, mortality, and labor.
        \item Amenity damages enter additively into utility, reflecting non-market welfare losses (e.g., mortality risk, extreme heat).
    \end{itemize}
    \vspace{0.2cm}
    \begin{figure}
        \centering
        \includegraphics[width=0.32\textwidth]{../Output/big_nu/graph_perfect_foresight/damage_maps/damage_maps_damage_wage.png}\hfill
        \includegraphics[width=0.32\textwidth]{../Output/big_nu/graph_perfect_foresight/damage_maps/damage_maps_damage_housing.png}\hfill
        \includegraphics[width=0.32\textwidth]{../Output/big_nu/graph_perfect_foresight/damage_maps/damage_maps_damage_amenity.png}
    \end{figure}
\end{frame}

\section{Quantitative Results}

\begin{frame}{Anticipation smooths aggregate income and wealth over the transition}
    \begin{figure}
        \centering
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_perfect_foresight/dynamic_moments/dynamic_moments_national_earnings.png}\hfill
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_perfect_foresight/dynamic_moments/dynamic_moments_national_wealth.png}
    \end{figure}
\end{frame}

\begin{frame}{Climate change reshapes population and wealth distributions}
    \begin{figure}
        \centering
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_perfect_foresight/transition_maps/transition_maps_local_pop_changes_pct.png}\hfill
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_perfect_foresight/transition_maps/transition_maps_wealth_income_ratio_changes.png}
    \end{figure}
\end{frame}

\begin{frame}{Spatial reallocation offsets aggregate losses but creates inequality}
    \begin{figure}
        \centering
        \includegraphics[width=0.9\textwidth]{../Output/big_nu/graph_perfect_foresight/exact_percapita_decomposition_bar_no_counterfactual_no_counterfactual.png}
    \end{figure}
\end{frame}

\begin{frame}{Without foresight, workers suffer greater welfare losses}
    \begin{columns}
        \begin{column}{0.4\textwidth}
            \begin{itemize}
                \item When workers are myopic, they fail to accumulate precautionary savings ahead of the shocks.
                \item The transition path diverges significantly from the baseline.
                \item Convergence to the new steady state is delayed by $\sim100$ years.
            \end{itemize}
        \end{column}
        \begin{column}{0.6\textwidth}
            \begin{figure}
                \centering
                \includegraphics[width=0.8\textwidth]{../Output/big_nu/graph_myopia/dynamic_moments/dynamic_moments_national_wealth.png}\\
                \vspace{0.2cm}
                \includegraphics[width=0.8\textwidth]{../Output/big_nu/graph_myopia/dynamic_moments/dynamic_moments_national_wy.png}
            \end{figure}
        \end{column}
    \end{columns}
\end{frame}

\begin{frame}{Mobility and saving are substitute insurance mechanisms}
    \begin{columns}
        \begin{column}{0.4\textwidth}
            \begin{itemize}
                \item \textbf{No Saving:} Inter-state migration increases to compensate for lack of financial buffering.
                \item \textbf{No Migration:} Wealth accumulation surges as workers can no longer spatially self-insure.
                \item Shutting down migration homogenizes spatial wealth responses compared to the baseline.
            \end{itemize}
        \end{column}
        \begin{column}{0.6\textwidth}
            \begin{figure}
                \centering
                \includegraphics[width=0.8\textwidth]{../Output/big_nu/graph_no_saving_struct/transition_maps/transition_maps_local_pop_changes_pct.png}\\
                \vspace{0.2cm}
                \includegraphics[width=0.8\textwidth]{../Output/big_nu/graph_no_migration_did/did_dynamic_moments_diff_wy.png}
            \end{figure}
        \end{column}
    \end{columns}
\end{frame}

\begin{frame}{Comparison of Scenarios: Frictions amplify welfare losses}
    \begin{figure}
        \centering
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_welfare_audit/welfare_audit_all_scenarios_ce_pct_frictions.png}\hfill
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_macro_audit/macro_wy_frictions/macro_wy_frictions_macro_wy_delta_wy_pct.png}
    \end{figure}
\end{frame}

\begin{frame}{Conclusion}
    \begin{itemize}
        \item \textbf{Joint Adaptation:} Saving and migration are critical, interacting mechanisms for climate adaptation.
        \item \textbf{Distributional Shift:} While aggregate US losses appear moderate, climate change drastically reshapes the spatial distribution of wealth and population.
        \item \textbf{Foresight Matters:} Anticipating future damages allows for smoother transitions via precautionary saving, significantly reducing welfare drops.
        \item \textbf{Substitutability:} Restricting either spatial mobility or financial saving forces households to over-rely on the remaining mechanism, amplifying overall welfare losses.
    \end{itemize}
\end{frame}

\appendix

\begin{frame}{Appendix: Computational Algorithm}
    \begin{itemize}
        \item \textbf{Static Equilibrium:} Solved via an iterative fixed-point algorithm on the continuation value $\Xi_{\ell}(z)$.
        \item \textbf{Dynamic Transition:} 
        \begin{itemize}
            \item Guess a transition path of 500 years.
            \item Backward induction to find the sequence of value functions and policy rules.
            \item Forward simulation to track the distribution of population and wealth.
        \end{itemize}
    \end{itemize}
\end{frame}

\begin{frame}{Appendix: Additional Spatial Moments}
    \begin{figure}
        \centering
        \includegraphics[width=0.48\textwidth]{../Output/graph_comparison/out_sample_fit/out_sample_fit_model_mean_wealth.png}\hfill
        \includegraphics[width=0.48\textwidth]{../Output/graph_comparison/out_sample_fit/out_sample_fit_log_10_model_flow_share.png}
    \end{figure}
\end{frame}

\begin{frame}{Appendix: Uncertainty and Variance Shocks}
    \begin{figure}
        \centering
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_uncertainty_variance_shocks/variance_map.png}\hfill
        \includegraphics[width=0.48\textwidth]{../Output/big_nu/graph_uncertainty_variance_shocks/dynamic_moments/dynamic_moments_national_wealth.png}
    \end{figure}
\end{frame}

\end{document}
"""

with open('presentation.tex', 'w') as f:
    f.write(tex_content)
