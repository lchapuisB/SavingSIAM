function [mdl, final_data] = test_amenities_correlation(P)
% TEST_AMENITIES_COMPONENTS Projette l'amenité structurelle A_l sur un 
% vecteur de caractéristiques géographiques

    data_dir = '/Users/leon/Documents/SavingSIAM/Output/';
    usda_tb  = readtable(fullfile(data_dir, 'usda_amenities.csv'));
    edu_tb   = readtable(fullfile(data_dir, 'education_amenities.csv'));
    efna_tb  = readtable(fullfile(data_dir, 'freedom_amenities.csv'));
    crime_tb = readtable(fullfile(data_dir, 'amenities_crime.csv'));
    urban_tb = readtable(fullfile(data_dir, 'amenities_variety.csv'));
    density_tb = readtable(fullfile(data_dir, 'amenities_density.csv'));

 
    A_calibrated = P.A(:);
    state_names = P.state_names(:);
    model_results = table(state_names, A_calibrated, ...
        'VariableNames', {'state_name', 'calibrated_amenity'});

   
    final_data = innerjoin(model_results, usda_tb, 'Keys', 'state_name');
    final_data = innerjoin(final_data, edu_tb, 'Keys', 'state_name');
    final_data = innerjoin(final_data, efna_tb, 'Keys', 'state_name');
    final_data = innerjoin(final_data, crime_tb, 'Keys', 'state_name');
    final_data = innerjoin(final_data, urban_tb, 'Keys', 'state_name');
    final_data = innerjoin(final_data, density_tb, 'Keys', 'state_name');

    formula = 'calibrated_amenity ~ pop_weighted_amenity + edu_amenity_index + efna_z + crime_z + urban_z + urban_density_z';
    mdl = fitlm(final_data, formula, 'RobustOpts', 'on');
    disp(mdl);

    figure('Name', 'Decomposition of the Structural Amenities', 'Color', 'w', ...
           'Position', [100, 100, 1000, 600]);

    vars_to_plot = {'pop_weighted_amenity', 'edu_amenity_index', 'efna_z','crime_z', 'urban_z' };
    titles = {'USDA (Natural Amenities)', 'Education (NAEP)', 'Economic Freedom (EFNA)', 'Crime', 'Urban Amenities (CBP)'};

    for i = 1:5
        subplot(3, 2,  i);
        plotAdded(mdl, vars_to_plot{i});
        title(sprintf('Partial Correlation : %s', titles{i}));
        grid on;
    end
end