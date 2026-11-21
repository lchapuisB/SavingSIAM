function ax = subplot(varargin)
    % Overloaded subplot to create independent figures instead of a grid.
    % This satisfies the requirement to save each plot separately into a subfolder.
    global my_sub_figs;
    
    fig = figure('Visible', 'off', 'Position', [100 100 1200 900]);
    
    if isempty(my_sub_figs)
        my_sub_figs = fig;
    else
        my_sub_figs(end+1) = fig;
    end
    
    % Return the axes of this new figure
    ax = gca;
end
