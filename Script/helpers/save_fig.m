function save_fig(fig, output_dir, filename)
    global my_sub_figs; %#ok<GVMIS>
    
    if ~isempty(my_sub_figs)
        [~, name, ext] = fileparts(filename);
        if isempty(ext), ext = '.png'; end
        
        sub_dir = fullfile(output_dir, name);
        if ~exist(sub_dir, 'dir')
            mkdir(sub_dir);
        end
        
        for i = 1:length(my_sub_figs)
            subfig = my_sub_figs(i);
            if ishandle(subfig)
                sub_tag = '';
                if isprop(subfig, 'Tag') && ~isempty(subfig.Tag)
                    sub_tag = subfig.Tag;
                else
                    all_axes = findall(subfig, 'Type', 'axes');
                    for a = 1:length(all_axes)
                        tmp_ax = all_axes(a);
                        if isprop(tmp_ax, 'Tag') && ~isempty(tmp_ax.Tag)
                            sub_tag = tmp_ax.Tag;
                            break;
                        elseif ~isempty(tmp_ax.YLabel.String)
                            sub_tag = char(tmp_ax.YLabel.String);
                            break;
                        elseif ~isempty(tmp_ax.XLabel.String)
                            sub_tag = char(tmp_ax.XLabel.String);
                            break;
                        end
                    end
                end
                if isempty(sub_tag)
                    txt_objs = findall(subfig, 'Type', 'Text');
                    for t = 1:length(txt_objs)
                        str = get(txt_objs(t), 'String');
                        if ~isempty(str)
                            if iscell(str), str = str{1}; end
                            sub_tag = str;
                            break;
                        end
                    end
                end
                
                if isempty(sub_tag)
                    sub_tag = num2str(i);
                else
                    sub_tag = regexprep(sub_tag, '[^a-zA-Z0-9]', '_');
                    sub_tag = regexprep(sub_tag, '_+', '_');
                    sub_tag = regexprep(sub_tag, '^_+|_+$', '');
                    sub_tag = lower(sub_tag);
                    if length(sub_tag) > 30
                        sub_tag = sub_tag(1:30);
                    end
                end
                
                sub_name = sprintf('%s_%s%s', name, sub_tag, ext);
                
                % Force la finalisation de la mise en page hors-écran
                drawnow;
                
                try
                    exportgraphics(subfig, fullfile(sub_dir, sub_name), 'Resolution', 150, 'ContentType', 'image');
                catch
                    try
                        print(subfig, fullfile(sub_dir, sub_name), '-dpng', '-r150');
                    catch
                    end
                end
                delete(subfig);
            end
        end
        
        my_sub_figs = [];
        drawnow limitrate;
        
        if ishandle(fig)
            delete(fig);
        end
    else
        if ishandle(fig)
            drawnow;
            
            try
                exportgraphics(fig, fullfile(output_dir, filename), 'Resolution', 150, 'ContentType', 'image');
            catch
                try
                    print(fig, fullfile(output_dir, filename), '-dpng', '-r150');
                catch
                end
            end
            delete(fig);
            drawnow limitrate;
        end
    end
end