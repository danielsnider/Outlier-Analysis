function color = experiment_to_color(experiment_name)
    color_map = containers.Map;
    color_map('Ctrl 1') = 'b';
    color_map('Ctrl 2') = 'c';
    color_map('Low dox 1') = 'k';
    color_map('Low dox 2') = 'g';
    color_map('High dox 1') = 'r';
    color_map('High dox 2') = 'm';

    color = color_map(experiment_name);
end

