

colour_scheme;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};

figure('Position', [680 558 680 420]);

% Formatting
axes_properties;

for iRat = 1:length(rat_names)
    
    rat = rat_names{iRat};
    
    % Load behavioural results
    load(['./data/behavioural_data/behaviour_' rat '.mat'])
    
    % Get the arm identity chosen on each trial
    evalc(['actions = cellfun(@(x, y) arrayfun(@(z) x(z), y), ' rat '.reward_probs, ' rat '.actions, ''UniformOutput'', false)']);
    
    % Plot
    subplot(2, 3, iRat); hold on
    line([0 23], [1/3 1/3], 'Color', colourscheme.chance, 'LineStyle', '--', 'HandleVisibility', 'off')
    plot(cellfun(@(x) mean(strcmp(x, 'high')), actions), 'Color', colourscheme.high, 'LineWidth', 2)
    plot(cellfun(@(x) mean(strcmp(x, 'medium')), actions), 'Color', colourscheme.medium, 'LineWidth', 2)
    plot(cellfun(@(x) mean(strcmp(x, 'low')), actions), 'Color', colourscheme.low, 'LineWidth', 2)
    
    % Title
    text(2, 0.65, ['Rat ' rat(1)], 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Formatting
    box off;
    xlim([0.5 22.5])
    ylim([0 0.7])
    
    % Axis labels
    set(gca, 'XTick', [1 22])
    if iRat <= 3
        set(gca, 'XTickLabel', [])
    end
    if iRat==5 xlabel('Training session'); end
    set(gca, 'YTick', [0 1/3, 0.7]);
    if iRat == 1 || iRat == 4
        set(gca, 'YTickLabel', {'0%', '33%', '70%'});
    else
        set(gca, 'YTickLabel', [])
    end
    
end

% Axis label
text(-70, 1, 'Percentage of arm choices', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 16, 'FontName', 'Arial')

% Legend
l = legend('High prob. arm', 'Mid prob. arm', 'Low prob. arm');
l.Box = 'off';
l.FontSize = 14;
l.Position = [0.7 0.45 0.2 0.15];