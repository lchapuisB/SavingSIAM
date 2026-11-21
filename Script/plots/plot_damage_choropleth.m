function plot_choropleth(state_names, data_vec, ax, global_bounds)
    set_aesthetics();
    
    parent_fig = get(ax, 'Parent');
    set(parent_fig, 'CurrentAxes', ax, 'Visible', 'off');

    % =====================================================================
    % 1. BOUNDS DETERMINATION
    % =====================================================================
    valid_data = data_vec(~isnan(data_vec) & ~isinf(data_vec));
    if isempty(valid_data)
        local_min = -1e-5;
        local_max =  1e-5;
    else
        local_min = prctile(valid_data, 10); 
        local_max = prctile(valid_data, 90);
    end

    if nargin >= 4 && ~isempty(global_bounds)
        min_val = global_bounds(1);
        max_val = global_bounds(2);
    else
        min_val = local_min;
        max_val = local_max;
    end
    
    cmax = max(abs(min_val), abs(max_val));
    if cmax == 0
        cmax = 1e-5;
    end

    % =====================================================================
    % 2. SMOOTH-TRANSITION DIVERGING COLORMAP
    % =====================================================================
    num_colors = 256;
    x_orig = [-1, -0.6, -0.25, -0.05, 0, 0.05, 0.25, 0.6, 1]; 
    
    % Gamma compression is set to 1 (linear) to preserve visual differences 
    % at the extremes (e.g., between -5 and -20).
    gamma = 0.5; 
    x_anchors = sign(x_orig) .* (abs(x_orig) .^ gamma);
    
    color_anchors = [
        0.75, 0.30, 0.20;  % Punchy Terracotta
        0.85, 0.50, 0.35;  % Warm Orange
        0.92, 0.65, 0.45;  % Sand
        0.96, 0.82, 0.60;  % Pale Wheat
        0.98, 0.98, 0.90;  % Cream / Near White
        0.82, 0.88, 0.60;  % Pale Olive
        0.65, 0.75, 0.45;  % Light Moss Green
        0.45, 0.60, 0.30;  % Moss Green
        0.25, 0.45, 0.15   % Punchy Forest Green
    ];

    xi = linspace(-1, 1, num_colors);
    cmap = interp1(x_anchors, color_anchors, xi, 'pchip');
    cmap = max(min(cmap, 1), 0);
    
    % =====================================================================
    % 3. MAPPING LOGIC 
    % =====================================================================
    states = shaperead('usastatelo', 'UseGeoCoords', true);
    
    axesm('eqaconic', 'MapLatLimit', [24 50], 'MapLonLimit', [-125 -66], ...
          'MapParallels', [29.5 45.5]);

    framem off; gridm off; mlabel off; plabel off; tightmap; 
    
    us_poly = polyshape();
    for i = 1:length(states)
        s_name = states(i).Name;
        if ~strcmpi(s_name, 'Alaska') && ~strcmpi(s_name, 'Hawaii')
            [lat_dense, lon_dense] = interpm(states(i).Lat, states(i).Lon, 0.1);
            [lon_c, lat_c] = polysplit(lon_dense, lat_dense);
            for j = 1:length(lon_c)
                us_poly = union(us_poly, polyshape(lon_c{j}, lat_c{j}, 'Simplify', false));
            end
        end
    end
    
    % --- RENDER STATES  ---
    for i = 1:length(states)
        s_name = states(i).Name;
        idx = find(strcmpi(state_names, s_name));
        
        if ~isempty(idx) && ~strcmpi(s_name, 'Alaska') && ~strcmpi(s_name, 'Hawaii')
            val = data_vec(idx);
            if ~isnan(val) && ~isinf(val)
                val_capped = max(min(val, cmax), -cmax);
                c_idx = round(1 + (num_colors - 1) * (val_capped - (-cmax)) / (2 * cmax));
                c_idx = max(1, min(num_colors, c_idx)); 
                
                [lat_dense, lon_dense] = interpm(states(i).Lat, states(i).Lon, 0.1);
                geoshow(lat_dense, lon_dense, 'DisplayType', 'polygon', ...
                    'FaceColor', cmap(c_idx, :), 'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 0.1);
            end
        end
    end
    
    % --- RENDER NATIONAL BOUNDARY ---
    us_poly = rmholes(us_poly);
    [lon_bnd, lat_bnd] = boundary(us_poly);
    [lat_bnd_dense, lon_bnd_dense] = interpm(lat_bnd, lon_bnd, 0.1);
    geoshow(lat_bnd_dense, lon_bnd_dense, 'DisplayType', 'polygon', 'FaceColor', 'none', ...
        'EdgeColor', 'k', 'LineWidth', 0.3);
    
    colormap(ax, cmap);
    set(ax, 'XColor', 'none', 'YColor', 'none');
    axis(ax, 'off');
    
    try
        clim([-cmax, cmax]); 
    catch
        caxis([-cmax, cmax]); 
    end

    pos = get(ax, 'Position'); 
    
    % =====================================================================
    cb = colorbar(ax, 'Location', 'eastoutside');
    if local_max > local_min
        cb.Limits = [local_min, local_max];
    else
        cb.Limits = [local_min - 1e-5, local_min + 1e-5];
    end
    set(cb, 'Position', [pos(1) + pos(3) * 0.98, pos(2) + pos(4) * 0.20, 0.025, pos(4) * 0.6]);
    cb_font = round(get(groot, 'defaultAxesFontSize') );
    set(cb, 'FontSize', cb_font, 'LineWidth', 0.5);
    

    % =====================================================================
    % 5. INSETS (ALASKA & HAWAÏ)
    % =====================================================================
    % --- ALASKA ---
    ak_pos = [pos(1) + pos(3) * 0.05, pos(2) + pos(4) * 0.06, pos(3) * 0.20, pos(4) * 0.20];
    ax_ak = axes('Parent', parent_fig, 'Position', ak_pos, 'Visible', 'off');
    axesm('eqaconic', 'MapLatLimit', [51 72], 'MapLonLimit', [-180 -129], ...
          'MapParallels', [55 65]);
    framem off; gridm off; mlabel off; plabel off; tightmap;
    set(ax_ak, 'XColor', 'none', 'YColor', 'none');
    axis(ax_ak, 'off');
    
    for i = 1:length(states)
        if strcmpi(states(i).Name, 'Alaska')
            idx = find(strcmpi(state_names, 'Alaska'));
            if ~isempty(idx)
                val = data_vec(idx);
                if ~isnan(val) && ~isinf(val)
                    val_capped = max(min(val, cmax), -cmax);
                    c_idx = round(1 + (num_colors - 1) * (val_capped - (-cmax)) / (2 * cmax));
                    c_idx = max(1, min(num_colors, c_idx)); 
                    
                    [lat_dense, lon_dense] = interpm(states(i).Lat, states(i).Lon, 0.1);
                    geoshow(lat_dense, lon_dense, 'DisplayType', 'polygon', ...
                        'FaceColor', cmap(c_idx, :), 'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 0.1);
                    geoshow(lat_dense, lon_dense, 'DisplayType', 'polygon', ...
                        'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.3);
                end
            end
        end
    end
    
    % --- HAWAÏ ---
    hi_pos = [pos(1) + pos(3) * 0.26, pos(2) + pos(4) * 0.06, pos(3) * 0.11, pos(4) * 0.11];
    ax_hi = axes('Parent', parent_fig, 'Position', hi_pos, 'Visible', 'off');
    axesm('eqaconic', 'MapLatLimit', [18 23], 'MapLonLimit', [-161 -154], ...
          'MapParallels', [19 22]);
    framem off; gridm off; mlabel off; plabel off; tightmap;
    set(ax_hi, 'XColor', 'none', 'YColor', 'none');
    axis(ax_hi, 'off');
    
    for i = 1:length(states)
        if strcmpi(states(i).Name, 'Hawaii')
            idx = find(strcmpi(state_names, 'Hawaii'));
            if ~isempty(idx)
                val = data_vec(idx);
                if ~isnan(val) && ~isinf(val)
                    val_capped = max(min(val, cmax), -cmax);
                    c_idx = round(1 + (num_colors - 1) * (val_capped - (-cmax)) / (2 * cmax));
                    c_idx = max(1, min(num_colors, c_idx)); 
                    
                    [lat_dense, lon_dense] = interpm(states(i).Lat, states(i).Lon, 0.1);
                    geoshow(lat_dense, lon_dense, 'DisplayType', 'polygon', ...
                        'FaceColor', cmap(c_idx, :), 'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 0.1);
                    geoshow(lat_dense, lon_dense, 'DisplayType', 'polygon', ...
                        'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.3);
                end
            end
        end
    end
    
    set(parent_fig, 'Visible', 'off');
end