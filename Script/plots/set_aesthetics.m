function set_aesthetics(mode)
    % Global settings for SavingSIAM MATLAB plots to make them more aesthetic,
    % use pastel colors, thinner lines/markers, and LaTeX style fonts.
    % mode: 'standalone' (default) for full-width figures,
    %       'subplot' for images assembled in LaTeX minipages (~0.4\textwidth).
    
    if nargin < 1, mode = 'standalone'; end
    
    % 1. Line Width, Marker Size, and Font Sizes (mode-dependent)
    switch mode
        case 'subplot'
            % Larger sizes to remain readable after LaTeX shrinks to ~40% width
            set(groot, 'defaultLineLineWidth', 4.0);
            set(groot, 'defaultLineMarkerSize', 10);
            set(groot, 'defaultAxesFontSize', 22);
            set(groot, 'defaultTextFontSize', 22);
            set(groot, 'defaultLegendFontSize', 22);
        otherwise % 'standalone'
            set(groot, 'defaultLineLineWidth', 3.0);
            set(groot, 'defaultLineMarkerSize', 8);
            set(groot, 'defaultAxesFontSize', 15);
            set(groot, 'defaultTextFontSize', 15);
            set(groot, 'defaultLegendFontSize', 15);
    end
    
    % 2. Font configuration (to match .tex style)
    set(groot, 'defaultTextInterpreter', 'tex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
    set(groot, 'defaultLegendInterpreter', 'tex');
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultTextFontName', 'Times New Roman');
    set(groot, 'defaultLegendFontName', 'Times New Roman');
    
    % 3. Grid configuration (small fine grey grid)
    set(groot, 'defaultAxesGridColor', [0.8 0.8 0.8]);
    set(groot, 'defaultAxesGridAlpha', 0.4);
    set(groot, 'defaultAxesMinorGridColor', [0.85 0.85 0.85]);
    set(groot, 'defaultAxesMinorGridAlpha', 0.3);
    
    % 4. Warm Academic Color Palette (Custom Color Order)
    pastel_colors = [
        0.33, 0.47, 0.65;  % Steel Blue
        0.82, 0.42, 0.35;  % Terracotta
        0.40, 0.65, 0.55;  % Sage Green
        0.85, 0.62, 0.35;  % Warm Amber
        0.58, 0.44, 0.65;  % Muted Lavender
        0.45, 0.68, 0.72;  % Teal Mist
        0.72, 0.65, 0.42   % Olive Gold
    ];
    set(groot, 'defaultAxesColorOrder', pastel_colors);
    
    set(groot, 'defaultAxesBox', 'on');
end
